#include <Rcpp.h>
#include <map>
#include <cmath>
#include <string>

using namespace Rcpp;

// =======================
// Helper: numeric matching
// =======================

// Match a numeric value x in vector vec with tolerance tol.
// Returns 0-based index, or -1 if no match.
int match_num(double x, const NumericVector& vec, double tol = 1e-8) {
  int n = vec.size();
  for (int i = 0; i < n; ++i) {
    double v = vec[i];
    if (NumericVector::is_na(v)) continue;
    if (std::fabs(x - v) < tol) return i;
  }
  return -1;
}

// Same, but in a given matrix column (0-based).
int match_num_in_col(double x, const NumericMatrix& mat, int col, double tol = 1e-8) {
  int n = mat.nrow();
  for (int i = 0; i < n; ++i) {
    double v = mat(i, col);
    if (NumericVector::is_na(v)) continue;
    if (std::fabs(x - v) < tol) return i;
  }
  return -1;
}

// ===========================
// Helper: name → index mapping
// ===========================

// Build a std::map from string name to 0-based index.
std::map<std::string, int> build_name_index(const CharacterVector& names) {
  std::map<std::string, int> m;
  int n = names.size();
  for (int i = 0; i < n; ++i) {
    // Convert Rcpp::String to std::string once here
    std::string s = as<std::string>(names[i]);
    m[s] = i;
  }
  return m;
}

// Lookup function with nice error if not found.
int get_index_or_stop(const std::map<std::string,int>& m,
                      const std::string& key,
                      const std::string& context) {
  std::map<std::string,int>::const_iterator it = m.find(key);
  if (it == m.end()) {
    stop(("Name '" + key + "' not found in " + context).c_str());
  }
  return it->second;
}


// =======================
// Main function: res_cpp
// =======================

//' Compute residuals (Rcpp implementation, optimized)
//'
//' @param data DataFrame with columns: time, name, value, sigma, lloq
//' @param out NumericMatrix with predictions (first column = time)
//' @param err Nullable NumericMatrix with error model predictions
//' @return List with:
//'   - data: DataFrame with residual columns
//'   - deriv:  NumericMatrix n_data × n_pars (first derivatives of out)
//'   - deriv2: NumericArray  n_data × n_pars × n_pars (second derivatives of out)
//'   - deriv.err:  NumericMatrix n_data × n_pars (first derivatives of err sigma-part)
//'   - deriv2.err: NumericArray  n_data × n_pars × n_pars (second derivatives of err sigma-part)
//' @export
// [[Rcpp::export]]
List res_cpp(DataFrame data, 
            NumericMatrix out,
            Nullable<NumericMatrix> err = R_NilValue) {
 
 // -------------------------
 // 1. Extract data columns
 // -------------------------
 NumericVector   data_time  = data["time"];
 CharacterVector data_name  = data["name"];
 NumericVector   data_value = data["value"];
 NumericVector   data_sigma = data["sigma"];
 NumericVector   data_lloq  = data["lloq"];
 int n_data = data_time.size();
 
 // Basic checks
 if (data_name.size() != n_data ||
     data_value.size() != n_data ||
     data_sigma.size() != n_data ||
     data_lloq.size()  != n_data) {
   stop("All data columns must have the same length.");
 }
 
 // -------------------------
 // 2. 'out' structure & maps
 // -------------------------
 int n_out_rows = out.nrow();
 int n_out_cols = out.ncol();
 
 if (n_out_cols < 2) {
   stop("'out' must have at least one time column and one observable column.");
 }
 
 CharacterVector out_colnames = colnames(out);
 if (out_colnames.size() != n_out_cols) {
   stop("'out' must have column names for all columns.");
 }
 
 // Time vector from 'out'
 NumericVector out_times(n_out_rows);
 for (int i = 0; i < n_out_rows; ++i) {
   out_times[i] = out(i, 0);
 }
 
 // Name → index map for 'out' columns
 std::map<std::string,int> out_name_map = build_name_index(out_colnames);
 
 // -------------------------
 // 3. Map each data row to 'out' row/col
 // -------------------------
 IntegerVector idx_time_out(n_data); // row index in out (and in deriv / deriv2)
 IntegerVector idx_name_out(n_data); // column index in out
 
 // Cache for time matching (data time → out index) to avoid repeated linear scans
 std::map<double,int> time_out_cache;
 
 for (int i = 0; i < n_data; ++i) {
   double t = data_time[i];
   
   // Time index: use cache if possible
   int ti;
   std::map<double,int>::const_iterator hit = time_out_cache.find(t);
   if (hit != time_out_cache.end()) {
     ti = hit->second;
   } else {
     ti = match_num(t, out_times);
     if (ti < 0) {
       stop(("Time point from data (time = " + std::to_string(t) +
         ") not found in 'out' time column.").c_str());
     }
     time_out_cache[t] = ti;
   }
   idx_time_out[i] = ti;
   
   // Name index: map via out_name_map
   std::string nm = as<std::string>(data_name[i]);
   int ni = get_index_or_stop(out_name_map, nm,
                              "column names of 'out'");
   idx_name_out[i] = ni;
 }
 
 // -------------------------
 // 4. Extract predictions
 // -------------------------
 NumericVector prediction(n_data);
 for (int i = 0; i < n_data; ++i) {
   prediction[i] = out(idx_time_out[i], idx_name_out[i]);
 }
 
 // -------------------------
 // 5. First derivatives of 'out' -> [n_data, n_pars]
 // -------------------------
 SEXP deriv_out = R_NilValue;
 
 if (out.hasAttribute("deriv")) {
   NumericVector deriv_array = out.attr("deriv");
   IntegerVector deriv_dims  = deriv_array.attr("dim");
   
   if (deriv_dims.size() != 3) {
     stop("'out$deriv' must be a 3D array [time, observable, parameter].");
   }
   
   int n_times_d = deriv_dims[0];
   int n_vars    = deriv_dims[1];
   int n_pars    = deriv_dims[2];
   
   if (n_times_d != n_out_rows) {
     stop("First dimension of 'out$deriv' must match number of rows in 'out'.");
   }
   
   // Parameter names
   CharacterVector pars;
   if (deriv_array.hasAttribute("dimnames")) {
     List dn = deriv_array.attr("dimnames");
     if (dn.size() >= 3 && !Rf_isNull(dn[2])) {
       pars = as<CharacterVector>(dn[2]);
     }
   }
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; ++p) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   // Observable names for second dimension
   CharacterVector deriv_varnames;
   if (deriv_array.hasAttribute("dimnames")) {
     List dn = deriv_array.attr("dimnames");
     if (dn.size() >= 2 && !Rf_isNull(dn[1])) {
       deriv_varnames = as<CharacterVector>(dn[1]);
     }
   }
   if (deriv_varnames.size() != n_vars) {
     stop("Second dimension of 'out$deriv' must have dimnames for observables.");
   }
   
   // Map observable name → index in deriv
   std::map<std::string,int> deriv_name_map = build_name_index(deriv_varnames);
   
   // Result matrix [n_data, n_pars]
   NumericMatrix deriv_result(n_data, n_pars);
   
   // Pointer to underlying raw data for speed
   double* dptr = REAL(deriv_array);
   
   // For each data row: pick correct [time, observable, :]
   for (int i = 0; i < n_data; ++i) {
     int time_idx = idx_time_out[i];
     std::string var = as<std::string>(data_name[i]);
     int var_idx = get_index_or_stop(deriv_name_map, var,
                                     "dimnames[[2]] of 'out$deriv'");
     
     // Base offset for [time_idx, var_idx, 0]
     // Linear index: t + T * (v + V * p)
     int base_tv = time_idx + n_times_d * var_idx;
     
     for (int p = 0; p < n_pars; ++p) {
       int src_idx = base_tv + n_times_d * n_vars * p;
       deriv_result(i, p) = dptr[src_idx];
     }
   }
   
   colnames(deriv_result) = pars;
   deriv_out = deriv_result;
 }
 
 // -------------------------
 // 6. Second derivatives of 'out' -> [n_data, n_pars, n_pars]
 // -------------------------
 SEXP deriv2_out = R_NilValue;
 
 if (out.hasAttribute("deriv2")) {
   NumericVector deriv2_array = out.attr("deriv2");
   IntegerVector deriv2_dims  = deriv2_array.attr("dim");
   
   if (deriv2_dims.size() != 4) {
     stop("'out$deriv2' must be a 4D array [time, observable, par, par].");
   }
   
   int n_times_d = deriv2_dims[0];
   int n_vars    = deriv2_dims[1];
   int n_pars    = deriv2_dims[2];
   int n_pars2   = deriv2_dims[3];
   
   if (n_times_d != n_out_rows) {
     stop("First dimension of 'out$deriv2' must match number of rows in 'out'.");
   }
   
   // Parameter names
   CharacterVector pars;
   if (deriv2_array.hasAttribute("dimnames")) {
     List dn = deriv2_array.attr("dimnames");
     if (dn.size() >= 3 && !Rf_isNull(dn[2])) {
       pars = as<CharacterVector>(dn[2]);
     }
   }
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; ++p) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   // Observable names
   CharacterVector deriv2_varnames;
   if (deriv2_array.hasAttribute("dimnames")) {
     List dn = deriv2_array.attr("dimnames");
     if (dn.size() >= 2 && !Rf_isNull(dn[1])) {
       deriv2_varnames = as<CharacterVector>(dn[1]);
     }
   }
   if (deriv2_varnames.size() != n_vars) {
     stop("Second dimension of 'out$deriv2' must have dimnames for observables.");
   }
   
   // Name map for observables
   std::map<std::string,int> deriv2_name_map = build_name_index(deriv2_varnames);
   
   // Result array [n_data, n_pars, n_pars2]
   NumericVector deriv2_result(n_data * n_pars * n_pars2);
   IntegerVector result_dims = IntegerVector::create(n_data, n_pars, n_pars2);
   deriv2_result.attr("dim") = result_dims;
   
   double* d2ptr = REAL(deriv2_array);
   double* r2ptr = REAL(deriv2_result);
   
   for (int i = 0; i < n_data; ++i) {
     int time_idx = idx_time_out[i];
     std::string var = as<std::string>(data_name[i]);
     int var_idx = get_index_or_stop(deriv2_name_map, var,
                                     "dimnames[[2]] of 'out$deriv2'");
     
     // Offset base for (time, var)
     int base_tv = time_idx + n_times_d * var_idx;
     
     for (int p1 = 0; p1 < n_pars; ++p1) {
       for (int p2 = 0; p2 < n_pars2; ++p2) {
         // src index: t + T*(v + V*p1 + V*P*p2)
         int src_idx = base_tv
         + n_times_d * n_vars * p1
         + n_times_d * n_vars * n_pars * p2;
         
         // dest index: i + n_data*(p1 + n_pars*p2)
         int dest_idx = i
         + n_data * p1
         + n_data * n_pars * p2;
         
         r2ptr[dest_idx] = d2ptr[src_idx];
       }
     }
   }
   
   // dimnames: list(NULL, pars, pars)
   List result_dimnames = List::create(R_NilValue, pars, pars);
   deriv2_result.attr("dimnames") = result_dimnames;
   
   deriv2_out = deriv2_result;
 }
 
 // -------------------------
 // 7. Error model handling
 // -------------------------
 LogicalVector sNAIndex = is_na(data_sigma);
 bool has_err = err.isNotNull();
 
 if (any(sNAIndex).is_true() && !has_err) {
   stop("Some sigmas are NA and no errmodel exists. Please fix data$sigma or supply errmodel.");
 }
 
 NumericMatrix err_mat;
 CharacterVector err_colnames;
 IntegerVector idx_time_err(n_data, -1);
 IntegerVector idx_name_err(n_data, -1);
 std::map<double,int> time_err_cache;
 std::map<std::string,int> err_name_map;
 
 if (has_err) {
   err_mat = as<NumericMatrix>(err);
   int n_err_rows = err_mat.nrow();
   int n_err_cols = err_mat.ncol();
   
   if (n_err_cols < 2) {
     stop("'err' must have at least one time column and one observable column.");
   }
   
   err_colnames = colnames(err_mat);
   if (err_colnames.size() != n_err_cols) {
     stop("'err' must have column names for all columns.");
   }
   
   // Time vector for err
   NumericVector err_times(n_err_rows);
   for (int i = 0; i < n_err_rows; ++i) {
     err_times[i] = err_mat(i, 0);
   }
   
   // Name map for err columns
   err_name_map = build_name_index(err_colnames);
   
   // Map each data row to err row/col (if needed)
   for (int i = 0; i < n_data; ++i) {
     // Time
     double t = data_time[i];
     int ti;
     std::map<double,int>::const_iterator hit = time_err_cache.find(t);
     if (hit != time_err_cache.end()) {
       ti = hit->second;
     } else {
       ti = match_num(t, err_times);
       if (ti < 0) {
         stop(("Time point from data (time = " + std::to_string(t) +
           ") not found in 'err' time column.").c_str());
       }
       time_err_cache[t] = ti;
     }
     idx_time_err[i] = ti;
     
     // Name
     std::string nm = as<std::string>(data_name[i]);
     int ni = get_index_or_stop(err_name_map, nm,
                                "column names of 'err'");
     idx_name_err[i] = ni;
   }
   
   // Fill missing sigmas from err-model
   for (int i = 0; i < n_data; ++i) {
     if (sNAIndex[i]) {
       int ti = idx_time_err[i];
       int ni = idx_name_err[i];
       double err_pred = err_mat(ti, ni);
       
       if (NumericVector::is_na(err_pred)) {
         stop("errmodel predicts NA for some observables with is.na(data$sigma).");
       }
       data_sigma[i] = err_pred;
     }
   }
 }
 
 // -------------------------
 // 8. First derivatives of error model (deriv.err)
 // -------------------------
 SEXP deriv_err_out = R_NilValue;
 
 if (has_err && err_mat.hasAttribute("deriv")) {
   NumericVector deriv_err_array = err_mat.attr("deriv");
   IntegerVector deriv_err_dims  = deriv_err_array.attr("dim");
   
   if (deriv_err_dims.size() != 3) {
     stop("'err$deriv' must be a 3D array [time, observable, parameter].");
   }
   
   int n_times_d = deriv_err_dims[0];
   int n_vars    = deriv_err_dims[1];
   int n_pars    = deriv_err_dims[2];
   
   if (n_times_d != err_mat.nrow()) {
     stop("First dimension of 'err$deriv' must match number of rows in 'err'.");
   }
   
   CharacterVector pars;
   if (deriv_err_array.hasAttribute("dimnames")) {
     List dn = deriv_err_array.attr("dimnames");
     if (dn.size() >= 3 && !Rf_isNull(dn[2])) {
       pars = as<CharacterVector>(dn[2]);
     }
   }
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; ++p) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   CharacterVector deriv_err_varnames;
   if (deriv_err_array.hasAttribute("dimnames")) {
     List dn = deriv_err_array.attr("dimnames");
     if (dn.size() >= 2 && !Rf_isNull(dn[1])) {
       deriv_err_varnames = as<CharacterVector>(dn[1]);
     }
   }
   if (deriv_err_varnames.size() != n_vars) {
     stop("Second dimension of 'err$deriv' must have dimnames for observables.");
   }
   
   // Name map for observables in err deriv
   std::map<std::string,int> deriv_err_name_map = build_name_index(deriv_err_varnames);
   
   NumericMatrix deriv_err_result(n_data, n_pars);
   double* derr_ptr = REAL(deriv_err_array);
   
   for (int i = 0; i < n_data; ++i) {
     int time_idx = idx_time_err[i];
     std::string var = as<std::string>(data_name[i]);
     int var_idx = get_index_or_stop(deriv_err_name_map, var,
                                     "dimnames[[2]] of 'err$deriv'");
     
     int base_tv = time_idx + n_times_d * var_idx;
     
     for (int p = 0; p < n_pars; ++p) {
       int src_idx = base_tv + n_times_d * n_vars * p;
       double val = derr_ptr[src_idx];
       
       // As in R version:
       // - derivative is 0 if original sigma was not NA
       // - or if derivative is NA
       if (NumericVector::is_na(val) || !sNAIndex[i]) {
         deriv_err_result(i, p) = 0.0;
       } else {
         deriv_err_result(i, p) = val;
       }
     }
   }
   
   colnames(deriv_err_result) = pars;
   deriv_err_out = deriv_err_result;
 }
 
 // -------------------------
 // 9. Second derivatives of error model (deriv2.err)
 // -------------------------
 SEXP deriv2_err_out = R_NilValue;
 
 if (has_err && err_mat.hasAttribute("deriv2")) {
   NumericVector deriv2_err_array = err_mat.attr("deriv2");
   IntegerVector deriv2_err_dims  = deriv2_err_array.attr("dim");
   
   if (deriv2_err_dims.size() != 4) {
     stop("'err$deriv2' must be a 4D array [time, observable, par, par].");
   }
   
   int n_times_d = deriv2_err_dims[0];
   int n_vars    = deriv2_err_dims[1];
   int n_pars    = deriv2_err_dims[2];
   int n_pars2   = deriv2_err_dims[3];
   
   if (n_times_d != err_mat.nrow()) {
     stop("First dimension of 'err$deriv2' must match number of rows in 'err'.");
   }
   
   CharacterVector pars;
   if (deriv2_err_array.hasAttribute("dimnames")) {
     List dn = deriv2_err_array.attr("dimnames");
     if (dn.size() >= 3 && !Rf_isNull(dn[2])) {
       pars = as<CharacterVector>(dn[2]);
     }
   }
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; ++p) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   CharacterVector deriv2_err_varnames;
   if (deriv2_err_array.hasAttribute("dimnames")) {
     List dn = deriv2_err_array.attr("dimnames");
     if (dn.size() >= 2 && !Rf_isNull(dn[1])) {
       deriv2_err_varnames = as<CharacterVector>(dn[1]);
     }
   }
   if (deriv2_err_varnames.size() != n_vars) {
     stop("Second dimension of 'err$deriv2' must have dimnames for observables.");
   }
   
   std::map<std::string,int> deriv2_err_name_map = build_name_index(deriv2_err_varnames);
   
   NumericVector deriv2_err_result(n_data * n_pars * n_pars2);
   IntegerVector result_dims = IntegerVector::create(n_data, n_pars, n_pars2);
   deriv2_err_result.attr("dim") = result_dims;
   
   double* d2e_ptr = REAL(deriv2_err_array);
   double* r2e_ptr = REAL(deriv2_err_result);
   
   for (int i = 0; i < n_data; ++i) {
     int time_idx = idx_time_err[i];
     std::string var = as<std::string>(data_name[i]);
     int var_idx = get_index_or_stop(deriv2_err_name_map, var,
                                     "dimnames[[2]] of 'err$deriv2'");
     
     int base_tv = time_idx + n_times_d * var_idx;
     
     for (int p1 = 0; p1 < n_pars; ++p1) {
       for (int p2 = 0; p2 < n_pars2; ++p2) {
         int src_idx = base_tv
         + n_times_d * n_vars * p1
         + n_times_d * n_vars * n_pars * p2;
         
         int dest_idx = i
         + n_data * p1
         + n_data * n_pars * p2;
         
         double val = d2e_ptr[src_idx];
         
         if (NumericVector::is_na(val) || !sNAIndex[i]) {
           r2e_ptr[dest_idx] = 0.0;
         } else {
           r2e_ptr[dest_idx] = val;
         }
       }
     }
   }
   
   List result_dimnames = List::create(R_NilValue, pars, pars);
   deriv2_err_result.attr("dimnames") = result_dimnames;
   
   deriv2_err_out = deriv2_err_result;
 }
 
 // -------------------------
 // 10. Apply LLOQ constraint
 // -------------------------
 for (int i = 0; i < n_data; ++i) {
   if (data_value[i] < data_lloq[i]) {
     data_value[i] = data_lloq[i];
   }
 }
 
 // -------------------------
 // 11. Compute residuals
 // -------------------------
 LogicalVector is_bloq(n_data);
 NumericVector residual(n_data);
 NumericVector weighted_residual(n_data);
 NumericVector weighted_0(n_data);
 
 for (int i = 0; i < n_data; ++i) {
   is_bloq[i]          = (data_value[i] <= data_lloq[i]);
   residual[i]         = prediction[i] - data_value[i];
   weighted_residual[i]= residual[i] / data_sigma[i];
   weighted_0[i]       = prediction[i] / data_sigma[i];
 }
 
 // -------------------------
 // 12. Build output DataFrame
 // -------------------------
 DataFrame result = DataFrame::create(
   _["time"]              = data_time,
   _["name"]              = data_name,
   _["value"]             = data_value,
   _["sigma"]             = data_sigma,
   _["lloq"]              = data_lloq,
   _["prediction"]        = prediction,
   _["residual"]          = residual,
   _["weighted.residual"] = weighted_residual,
   _["weighted.0"]        = weighted_0,
   _["bloq"]              = is_bloq
 );
 
 // -------------------------
 // 13. Final list
 // -------------------------
 List output = List::create(
   _["data"]       = result,
   _["deriv"]      = deriv_out,
   _["deriv2"]     = deriv2_out,
   _["deriv.err"]  = deriv_err_out,
   _["deriv2.err"] = deriv2_err_out
 );
 
 return output;
}