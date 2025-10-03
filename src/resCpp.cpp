#include <Rcpp.h>
#include <unordered_map>
#include <sstream>
#include <iomanip>

using namespace Rcpp;

// Normal pdf / cdf
inline double phi(double x) { return R::dnorm4(x, 0.0, 1.0, false); }
inline double Phi(double x) { return R::pnorm5(x, 0.0, 1.0, true, false); }
inline double logPhiNeg(double x) {
  return R::pnorm5(-x, 0.0, 1.0, true, true);
}

// Build stable (time,name) key
inline std::string make_key(double t, const Rcpp::String& n) {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss << std::setprecision(17) << t << "|" << std::string(n);
  return oss.str();
}

template <typename VecT>
inline VecT get_df_col(const DataFrame& df, const char* nm) {
  if (!df.containsElementNamed(nm))
    stop("Column '%s' not found.", nm);
  return df[nm];
}

// Match numeric values (equivalent to match.num in R)
inline IntegerVector match_num(NumericVector x, NumericVector table) {
  int n = x.size();
  int m = table.size();
  IntegerVector result(n);
  
  for (int i = 0; i < n; ++i) {
    bool found = false;
    for (int j = 0; j < m; ++j) {
      if (std::abs(x[i] - table[j]) < 1e-12) {
        result[i] = j + 1; // R uses 1-based indexing
        found = true;
        break;
      }
    }
    if (!found) result[i] = NA_INTEGER;
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List resCpp(Rcpp::DataFrame data, SEXP out, Rcpp::Nullable<Rcpp::DataFrame> err = R_NilValue, std::string optBLOQ = "none") {
  
  // --- 1) Extract required columns from data ----
  NumericVector time_data  = get_df_col<NumericVector>(data, "time");
  CharacterVector name_data = get_df_col<CharacterVector>(data, "name");
  NumericVector value = clone(get_df_col<NumericVector>(data, "value"));
  NumericVector sigma = clone(get_df_col<NumericVector>(data, "sigma"));
  int n = data.nrows();
  
  // Check for lloq column
  NumericVector lloq;
  bool has_lloq = data.containsElementNamed("lloq");
  if (has_lloq) {
    lloq = data["lloq"];
  } else {
    lloq = NumericVector(n, R_NegInf);
  }
  
  // --- 2) Build prediction vector by matching (time,name) like res() does ----
  NumericVector prediction(n, NA_REAL);
  
  // Check if out is a data.frame or matrix
  if (Rf_isFrame(out)) {
    // out is a data.frame with time and named columns (like ODE output)
    DataFrame outdf = as<DataFrame>(out);
    
    if (outdf.containsElementNamed("time") && outdf.containsElementNamed("name") && 
        (outdf.containsElementNamed("prediction") || outdf.containsElementNamed("value"))) {
      // Long format: has time, name, prediction columns
      NumericVector ot = get_df_col<NumericVector>(outdf, "time");
      CharacterVector on = get_df_col<CharacterVector>(outdf, "name");
      
      NumericVector oy;
      if (outdf.containsElementNamed("prediction")) {
        oy = outdf["prediction"];
      } else {
        oy = outdf["value"];
      }
      
      const int m = outdf.nrows();
      std::unordered_map<std::string, double> pred_map;
      pred_map.reserve(static_cast<size_t>(m) * 2u);
      
      for (int k = 0; k < m; ++k) {
        pred_map.emplace(make_key(ot[k], on[k]), oy[k]);
      }
      
      for (int i = 0; i < n; ++i) {
        const std::string key = make_key(time_data[i], name_data[i]);
        auto it = pred_map.find(key);
        if (it == pred_map.end()) {
          stop("No prediction found in `out` for (time=%g, name='%s').",
               static_cast<double>(time_data[i]),
               std::string(name_data[i]).c_str());
        }
        prediction[i] = it->second;
      }
      
    } else {
      // Wide format: time in first column, observables in other columns (ODE output)
      NumericVector out_times = outdf[0];
      CharacterVector out_names = outdf.names();
      
      // Get unique times and names from data
      NumericVector unique_times = sort_unique(time_data);
      CharacterVector unique_names = sort_unique(name_data);
      
      // Match data times/names in unique times/names
      IntegerVector data_time_idx = match(time_data, unique_times);
      IntegerVector data_name_idx = match(name_data, unique_names);
      
      // Match unique times in out times
      IntegerVector out_time_idx = match_num(unique_times, out_times);
      
      // Match unique names in out column names (skip first column = time)
      IntegerVector out_name_idx(unique_names.size());
      for (int i = 0; i < unique_names.size(); ++i) {
        bool found = false;
        for (int j = 1; j < out_names.size(); ++j) { // skip first column
          if (std::string(unique_names[i]) == std::string(out_names[j])) {
            out_name_idx[i] = j + 1; // R 1-based
            found = true;
            break;
          }
        }
        if (!found) {
          stop("The following observable in data does not have a prediction: %s", 
               std::string(unique_names[i]).c_str());
        }
      }
      
      // Extract predictions
      for (int i = 0; i < n; ++i) {
        int time_idx = out_time_idx[data_time_idx[i] - 1] - 1; // to 0-based
        int name_idx = out_name_idx[data_name_idx[i] - 1] - 1; // to 0-based
        
        if (IntegerVector::is_na(out_time_idx[data_time_idx[i] - 1])) {
          stop("Time %g from data not found in out", static_cast<double>(time_data[i]));
        }
        
        NumericVector col = outdf[name_idx];
        prediction[i] = col[time_idx];
      }
    }
    
  } else if (Rf_isMatrix(out) && TYPEOF(out) == REALSXP) {
    // out is a numeric matrix: wide format with time in column 1
    NumericMatrix om = as<NumericMatrix>(out);
    
    if (om.ncol() < 2) {
      stop("Matrix `out` must have at least 2 columns (time + observables).");
    }
    
    NumericVector out_times = om(_, 0);
    CharacterVector col_names = colnames(om);
    
    if (col_names.size() != om.ncol()) {
      stop("Matrix `out` must have column names for observables.");
    }
    
    // Get unique times and names from data
    NumericVector unique_times = sort_unique(time_data);
    CharacterVector unique_names = sort_unique(name_data);
    
    // Match data times/names in unique times/names
    IntegerVector data_time_idx = match(time_data, unique_times);
    IntegerVector data_name_idx = match(name_data, unique_names);
    
    // Match unique times in out times
    IntegerVector out_time_idx = match_num(unique_times, out_times);
    
    // Match unique names in out column names (skip first column = time)
    IntegerVector out_name_idx(unique_names.size());
    for (int i = 0; i < unique_names.size(); ++i) {
      bool found = false;
      for (int j = 1; j < col_names.size(); ++j) { // skip first column
        if (std::string(unique_names[i]) == std::string(col_names[j])) {
          out_name_idx[i] = j; // 0-based for matrix indexing
          found = true;
          break;
        }
      }
      if (!found) {
        stop("The following observable in data does not have a prediction: %s", 
             std::string(unique_names[i]).c_str());
      }
    }
    
    // Extract predictions
    for (int i = 0; i < n; ++i) {
      int time_idx = out_time_idx[data_time_idx[i] - 1] - 1; // to 0-based
      int name_idx = out_name_idx[data_name_idx[i] - 1];     // already 0-based
      
      if (IntegerVector::is_na(out_time_idx[data_time_idx[i] - 1])) {
        stop("Time %g from data not found in out", static_cast<double>(time_data[i]));
      }
      
      prediction[i] = om(time_idx, name_idx);
    }
    
  } else {
    stop("`out` must be a data.frame or a numeric matrix.");
  }
  
  // --- 3) Handle error model if provided ----
  LogicalVector sNAIndex(n);
  for (int i = 0; i < n; ++i) {
    sNAIndex[i] = NumericVector::is_na(sigma[i]);
  }
  
  bool any_na_sigma = is_true(any(sNAIndex));
  
  if (any_na_sigma && err.isNull()) {
    stop("In data, some sigmas are NA and no errmodel exists. Please fix data$sigma or supply errmodel.");
  }
  
  if (!err.isNull() && any_na_sigma) {
    
    if (Rf_isFrame(err)) {
      DataFrame err_df = as<DataFrame>(err);
      
      if (err_df.containsElementNamed("time") && err_df.containsElementNamed("name") && 
          (err_df.containsElementNamed("prediction") || err_df.containsElementNamed("value"))) {
        // Long format error model
        NumericVector et = get_df_col<NumericVector>(err_df, "time");
        CharacterVector en = get_df_col<CharacterVector>(err_df, "name");
        NumericVector eval;
        
        if (err_df.containsElementNamed("prediction")) {
          eval = err_df["prediction"];
        } else {
          eval = err_df["value"];
        }
        
        std::unordered_map<std::string, double> err_map;
        for (int k = 0; k < err_df.nrows(); ++k) {
          err_map.emplace(make_key(et[k], en[k]), eval[k]);
        }
        
        for (int i = 0; i < n; ++i) {
          if (sNAIndex[i]) {
            const std::string key = make_key(time_data[i], name_data[i]);
            auto it = err_map.find(key);
            if (it == err_map.end() || NumericVector::is_na(it->second)) {
              stop("errmodel predicts NA for observable (time=%g, name='%s').",
                   static_cast<double>(time_data[i]),
                   std::string(name_data[i]).c_str());
            }
            sigma[i] = it->second;
          }
        }
      } else {
        // Wide format error model (time in first column, observables in other columns)
        NumericVector err_times = err_df[0];
        CharacterVector err_names = err_df.names();
        
        // Get unique times and names from data
        NumericVector unique_times = sort_unique(time_data);
        CharacterVector unique_names = sort_unique(name_data);
        
        // Match data times/names in unique times/names
        IntegerVector data_time_idx = match(time_data, unique_times);
        IntegerVector data_name_idx = match(name_data, unique_names);
        
        // Match unique times in err times
        IntegerVector err_time_idx = match_num(unique_times, err_times);
        
        // Match unique names in err column names (skip first column = time)
        IntegerVector err_name_idx(unique_names.size());
        for (int i = 0; i < unique_names.size(); ++i) {
          bool found = false;
          for (int j = 1; j < err_names.size(); ++j) { // skip first column
            if (std::string(unique_names[i]) == std::string(err_names[j])) {
              err_name_idx[i] = j; // 0-based for column indexing
              found = true;
              break;
            }
          }
          if (!found) {
            // If error model doesn't have this observable, that's OK - just skip it
            err_name_idx[i] = NA_INTEGER;
          }
        }
        
        // Extract error predictions for NA sigmas
        for (int i = 0; i < n; ++i) {
          if (sNAIndex[i]) {
            int time_idx_pos = data_time_idx[i] - 1; // to 0-based
            int name_idx_pos = data_name_idx[i] - 1; // to 0-based
            
            if (IntegerVector::is_na(err_name_idx[name_idx_pos])) {
              stop("errmodel does not contain observable '%s'", 
                   std::string(name_data[i]).c_str());
            }
            
            int err_t_idx = err_time_idx[time_idx_pos] - 1; // to 0-based
            int err_n_idx = err_name_idx[name_idx_pos];     // already 0-based
            
            if (IntegerVector::is_na(err_time_idx[time_idx_pos])) {
              stop("Time %g from data not found in errmodel", 
                   static_cast<double>(time_data[i]));
            }
            
            NumericVector err_col = err_df[err_n_idx];
            double err_val = err_col[err_t_idx];
            
            if (NumericVector::is_na(err_val)) {
              stop("errmodel predicts NA for observable (time=%g, name='%s').",
                   static_cast<double>(time_data[i]),
                   std::string(name_data[i]).c_str());
            }
            
            sigma[i] = err_val;
          }
        }
      }
    } else if (Rf_isMatrix(err) && TYPEOF(err) == REALSXP) {
      // Wide format error matrix
      NumericMatrix em = as<NumericMatrix>(err);
      
      if (em.ncol() < 2) {
        stop("Matrix `err` must have at least 2 columns (time + observables).");
      }
      
      NumericVector err_times = em(_, 0);
      CharacterVector col_names = colnames(em);
      
      if (col_names.size() != em.ncol()) {
        stop("Matrix `err` must have column names for observables.");
      }
      
      // Get unique times and names from data
      NumericVector unique_times = sort_unique(time_data);
      CharacterVector unique_names = sort_unique(name_data);
      
      // Match data times/names in unique times/names
      IntegerVector data_time_idx = match(time_data, unique_times);
      IntegerVector data_name_idx = match(name_data, unique_names);
      
      // Match unique times in err times
      IntegerVector err_time_idx = match_num(unique_times, err_times);
      
      // Match unique names in err column names (skip first column = time)
      IntegerVector err_name_idx(unique_names.size());
      for (int i = 0; i < unique_names.size(); ++i) {
        bool found = false;
        for (int j = 1; j < col_names.size(); ++j) { // skip first column
          if (std::string(unique_names[i]) == std::string(col_names[j])) {
            err_name_idx[i] = j; // 0-based for matrix indexing
            found = true;
            break;
          }
        }
        if (!found) {
          err_name_idx[i] = NA_INTEGER;
        }
      }
      
      // Extract error predictions
      for (int i = 0; i < n; ++i) {
        if (sNAIndex[i]) {
          int time_idx_pos = data_time_idx[i] - 1; // to 0-based
          int name_idx_pos = data_name_idx[i] - 1; // to 0-based
          
          if (IntegerVector::is_na(err_name_idx[name_idx_pos])) {
            stop("errmodel does not contain observable '%s'", 
                 std::string(name_data[i]).c_str());
          }
          
          int err_t_idx = err_time_idx[time_idx_pos] - 1; // to 0-based
          int err_n_idx = err_name_idx[name_idx_pos];     // already 0-based
          
          if (IntegerVector::is_na(err_time_idx[time_idx_pos])) {
            stop("Time %g from data not found in errmodel", 
                 static_cast<double>(time_data[i]));
          }
          
          double err_val = em(err_t_idx, err_n_idx);
          
          if (NumericVector::is_na(err_val)) {
            stop("errmodel predicts NA for observable (time=%g, name='%s').",
                 static_cast<double>(time_data[i]),
                 std::string(name_data[i]).c_str());
          }
          
          sigma[i] = err_val;
        }
      }
    } else {
      stop("`err` must be a data.frame or a numeric matrix.");
    }
  }
  
  // --- 4) Set value to lloq if below lloq (like res() does) ----
  LogicalVector bloq(n);
  for (int i = 0; i < n; ++i) {
    if (value[i] <= lloq[i]) {
      bloq[i] = true;
      value[i] = lloq[i];  // This is what res() does!
    } else {
      bloq[i] = false;
    }
  }
  
  // --- 5) Base residuals (like res()) ----
  NumericVector residual(n);
  NumericVector weighted_residual(n);
  NumericVector weighted0(n);
  
  for (int i = 0; i < n; ++i) {
    residual[i] = prediction[i] - value[i];
    weighted_residual[i] = residual[i] / sigma[i];
    weighted0[i] = prediction[i] / sigma[i];
  }
  
  // --- 6) Optional BLOQ handling (only if optBLOQ != "none") ----
  if (optBLOQ != "none") {
    for (int i = 0; i < n; ++i) {
      if (!bloq[i]) continue;
      
      const double wr = weighted_residual[i];
      const double w0 = weighted0[i];
      
      if (optBLOQ == "M1") {
        weighted_residual[i] = 0.0;
        
      } else if (optBLOQ == "M3") {
        double val = -2.0 * logPhiNeg(wr);
        if (!R_finite(val) || val < 0.0) val = 0.0;
        weighted_residual[i] = std::sqrt(val);
        
      } else if (optBLOQ == "M4NM") {
        if (w0 < 0.0) {
          stop("M4NM requires LLOQ >= 0 and prediction/sigma >= 0; consider M3 for negative LLOQ.");
        }
        const double Phi_wr = Phi(wr);
        const double Phi_w0 = Phi(w0);
        double frac = Phi_wr / Phi_w0;
        if (frac >= 1.0) frac = 1.0 - 1e-12;
        double val = -2.0 * std::log(1.0 - frac);
        if (!R_finite(val) || val < 0.0) val = 0.0;
        weighted_residual[i] = std::sqrt(val);
        
      } else {
        stop("Unknown optBLOQ '%s' (use 'none','M1','M3','M4NM').", optBLOQ.c_str());
      }
    }
  }
  
  // --- 7) Assemble objframe-like data.frame ----
  DataFrame outDf = DataFrame::create(
    Named("time") = time_data,
    Named("name") = name_data,
    Named("value") = value,
    Named("prediction") = prediction,
    Named("sigma") = sigma,
    Named("residual") = residual,
    Named("weighted.residual") = weighted_residual,
    Named("weighted.0") = weighted0,
    Named("bloq") = bloq
  );
  
  outDf.attr("class") = CharacterVector::create("objframe", "data.frame");
  return outDf;
}