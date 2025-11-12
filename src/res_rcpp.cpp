#include <Rcpp.h>
using namespace Rcpp;

//' Compute residuals (Rcpp implementation)
//' 
//' @param data DataFrame with columns: time, name, value, sigma, lloq
//' @param out NumericMatrix with predictions (first column = time)
//' @param err Nullable NumericMatrix with error model predictions
//' @return List with data and derivatives (3D/4D arrays)
//' @export
// [[Rcpp::export]]
List res_cpp(DataFrame data, 
             NumericMatrix out,
             Nullable<NumericMatrix> err = R_NilValue) {
 
  // Extract data columns
  NumericVector data_time = data["time"];
  CharacterVector data_name = data["name"];
  NumericVector data_value = data["value"];
  NumericVector data_sigma = data["sigma"];
  NumericVector data_lloq = data["lloq"];
  int n_data = data_time.size();
  
  // Get unique times and names
  NumericVector times = clone(data_time);
  times = sort_unique(times);
  CharacterVector names = unique(data_name);
  
  // Get out column names
  CharacterVector out_colnames = colnames(out);
  
  // Match data to unique times/names
  IntegerVector data_time_idx(n_data);
  IntegerVector data_name_idx(n_data);
  
  for (int i = 0; i < n_data; i++) {
   // Find time index (exact match)
   for (int t = 0; t < times.size(); t++) {
     if (data_time[i] == times[t]) {
       data_time_idx[i] = t;
       break;
     }
   }
   // Find name index (exact match)
   for (int j = 0; j < names.size(); j++) {
     if (data_name[i] == names[j]) {
       data_name_idx[i] = j;
       break;
     }
   }
  }
  
  // Match unique times in out (exact match)
  IntegerVector out_time_idx(times.size());
  for (int t = 0; t < times.size(); t++) {
   for (int i = 0; i < out.nrow(); i++) {
     if (times[t] == out(i, 0)) {
       out_time_idx[t] = i;
       break;
     }
   }
  }
  
  // Match names in out columns
  IntegerVector out_name_idx(names.size());
  for (int j = 0; j < names.size(); j++) {
   bool found = false;
   for (int k = 0; k < out_colnames.size(); k++) {
     if (names[j] == out_colnames[k]) {
       out_name_idx[j] = k;
       found = true;
       break;
     }
   }
   if (!found) {
     stop("Observable in data does not have a prediction: " + 
       std::string(names[j]));
   }
  }
  
  // Extract predictions
  NumericVector prediction(n_data);
  for (int i = 0; i < n_data; i++) {
   int time_idx = out_time_idx[data_time_idx[i]];
   int name_idx = out_name_idx[data_name_idx[i]];
   prediction[i] = out(time_idx, name_idx);
  }
  
  // Handle first derivatives if available (3D array format)
  SEXP deriv_out = R_NilValue;
  
  if (out.hasAttribute("deriv")) {
   NumericVector deriv_array = out.attr("deriv");
   IntegerVector deriv_dims = deriv_array.attr("dim");
   
   // deriv_dims: [n_times, n_vars, n_pars]
   int n_times_d = deriv_dims[0];
   int n_vars = deriv_dims[1];
   int n_pars = deriv_dims[2];
   
   // Get parameter names from dimnames[[3]]
   CharacterVector pars;
   List original_dimnames = R_NilValue;
   if (deriv_array.hasAttribute("dimnames")) {
     original_dimnames = deriv_array.attr("dimnames");
     if (original_dimnames.size() >= 3 && !Rf_isNull(original_dimnames[2])) {
       pars = as<CharacterVector>(original_dimnames[2]);
     }
   }
   
   // If no parameter names, create generic ones
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; p++) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   // Create 3D array: [n_data, n_vars_in_data, n_pars]
   // Get unique variable names in data
   CharacterVector data_vars = unique(data_name);
   int n_data_vars = data_vars.size();
   
   NumericVector deriv_result(n_data * n_data_vars * n_pars);
   IntegerVector result_dims = IntegerVector::create(n_data, n_data_vars, n_pars);
   
   // Fill the 3D array
   for (int i = 0; i < n_data; i++) {
     int time_idx = out_time_idx[data_time_idx[i]];
     int var_idx = data_name_idx[i];
     
     for (int p = 0; p < n_pars; p++) {
       // Access source 3D array: [time_idx, var_idx, p]
       int src_idx = time_idx + n_times_d * var_idx + n_times_d * n_vars * p;
       // Store in result 3D array: [i, var_idx, p]
       int dest_idx = i + n_data * var_idx + n_data * n_data_vars * p;
       deriv_result[dest_idx] = deriv_array[src_idx];
     }
   }
   
   deriv_result.attr("dim") = result_dims;
   
   // Set dimnames: list(NULL, data_vars, pars)
   List result_dimnames = List::create(R_NilValue, data_vars, pars);
   deriv_result.attr("dimnames") = result_dimnames;
   
   deriv_out = deriv_result;
  }
  
  // Handle second derivatives if available (4D array format)
  SEXP deriv2_out = R_NilValue;
  
  if (out.hasAttribute("deriv2")) {
   NumericVector deriv2_array = out.attr("deriv2");
   IntegerVector deriv2_dims = deriv2_array.attr("dim");
   
   // deriv2_dims: [n_times, n_vars, n_pars, n_pars]
   int n_times_d = deriv2_dims[0];
   int n_vars = deriv2_dims[1];
   int n_pars = deriv2_dims[2];
   int n_pars2 = deriv2_dims[3];
   
   // Get parameter names
   CharacterVector pars;
   if (deriv2_array.hasAttribute("dimnames")) {
     List original_dimnames = deriv2_array.attr("dimnames");
     if (original_dimnames.size() >= 3 && !Rf_isNull(original_dimnames[2])) {
       pars = as<CharacterVector>(original_dimnames[2]);
     }
   }
   
   if (pars.size() == 0) {
     pars = CharacterVector(n_pars);
     for (int p = 0; p < n_pars; p++) {
       pars[p] = "p" + std::to_string(p + 1);
     }
   }
   
   // Create 4D array: [n_data, n_vars_in_data, n_pars, n_pars]
   CharacterVector data_vars = unique(data_name);
   int n_data_vars = data_vars.size();
   
   NumericVector deriv2_result(n_data * n_data_vars * n_pars * n_pars2);
   IntegerVector result_dims = IntegerVector::create(n_data, n_data_vars, n_pars, n_pars2);
   
   // Fill the 4D array
   for (int i = 0; i < n_data; i++) {
     int time_idx = out_time_idx[data_time_idx[i]];
     int var_idx = data_name_idx[i];
     
     for (int p1 = 0; p1 < n_pars; p1++) {
       for (int p2 = 0; p2 < n_pars2; p2++) {
         // Access source 4D array: [time_idx, var_idx, p1, p2]
         int src_idx = time_idx + n_times_d * var_idx + 
           n_times_d * n_vars * p1 + 
           n_times_d * n_vars * n_pars * p2;
         // Store in result 4D array: [i, var_idx, p1, p2]
         int dest_idx = i + n_data * var_idx + 
           n_data * n_data_vars * p1 + 
           n_data * n_data_vars * n_pars * p2;
         deriv2_result[dest_idx] = deriv2_array[src_idx];
       }
     }
   }
   
   deriv2_result.attr("dim") = result_dims;
   
   // Set dimnames: list(NULL, data_vars, pars, pars)
   List result_dimnames = List::create(R_NilValue, data_vars, pars, pars);
   deriv2_result.attr("dimnames") = result_dimnames;
   
   deriv2_out = deriv2_result;
  }
  
  // Handle error model
  LogicalVector sNAIndex = is_na(data_sigma);
  bool has_err = err.isNotNull();
  
  if (any(sNAIndex).is_true() && !has_err) {
   stop("Some sigmas are NA and no errmodel exists. Please fix data$sigma or supply errmodel.");
  }
  
  if (has_err && any(sNAIndex).is_true()) {
   NumericMatrix err_mat = as<NumericMatrix>(err);
   CharacterVector err_colnames = colnames(err_mat);
   
   // Match times and names in err (exact match)
   IntegerVector err_time_idx(times.size());
   for (int t = 0; t < times.size(); t++) {
     for (int i = 0; i < err_mat.nrow(); i++) {
       if (times[t] == err_mat(i, 0)) {
         err_time_idx[t] = i;
         break;
       }
     }
   }
   
   IntegerVector err_name_idx(names.size());
   for (int j = 0; j < names.size(); j++) {
     for (int k = 0; k < err_colnames.size(); k++) {
       if (names[j] == err_colnames[k]) {
         err_name_idx[j] = k;
         break;
       }
     }
   }
   
   // Extract error predictions and fill NA sigmas
   for (int i = 0; i < n_data; i++) {
     if (sNAIndex[i]) {
       int time_idx = err_time_idx[data_time_idx[i]];
       int name_idx = err_name_idx[data_name_idx[i]];
       double err_pred = err_mat(time_idx, name_idx);
       
       if (NumericVector::is_na(err_pred)) {
         stop("errmodel predicts NA for some observables with is.na(data$sigma).");
       }
       data_sigma[i] = err_pred;
     }
   }
  }
  
  // Handle error model first derivatives (3D array format)
  SEXP deriv_err_out = R_NilValue;
  
  if (has_err) {
   NumericMatrix err_mat = as<NumericMatrix>(err);
   
   if (err_mat.hasAttribute("deriv")) {
     NumericVector deriv_err_array = err_mat.attr("deriv");
     IntegerVector deriv_err_dims = deriv_err_array.attr("dim");
     
     int n_times_d = deriv_err_dims[0];
     int n_vars = deriv_err_dims[1];
     int n_pars = deriv_err_dims[2];
     
     CharacterVector pars;
     if (deriv_err_array.hasAttribute("dimnames")) {
       List deriv_err_dimnames = deriv_err_array.attr("dimnames");
       if (deriv_err_dimnames.size() >= 3 && !Rf_isNull(deriv_err_dimnames[2])) {
         pars = as<CharacterVector>(deriv_err_dimnames[2]);
       }
     }
     
     if (pars.size() == 0) {
       pars = CharacterVector(n_pars);
       for (int p = 0; p < n_pars; p++) {
         pars[p] = "p" + std::to_string(p + 1);
       }
     }
     
     // Create 3D array for error derivatives
     CharacterVector data_vars = unique(data_name);
     int n_data_vars = data_vars.size();
     
     NumericVector deriv_err_result(n_data * n_data_vars * n_pars);
     IntegerVector result_dims = IntegerVector::create(n_data, n_data_vars, n_pars);
     
     // Match indices
     CharacterVector err_colnames = colnames(err_mat);
     IntegerVector err_time_idx(times.size());
     for (int t = 0; t < times.size(); t++) {
       for (int i = 0; i < err_mat.nrow(); i++) {
         if (times[t] == err_mat(i, 0)) {
           err_time_idx[t] = i;
           break;
         }
       }
     }
     
     for (int i = 0; i < n_data; i++) {
       int time_idx = err_time_idx[data_time_idx[i]];
       int var_idx = data_name_idx[i];
       
       for (int p = 0; p < n_pars; p++) {
         int src_idx = time_idx + n_times_d * var_idx + n_times_d * n_vars * p;
         int dest_idx = i + n_data * var_idx + n_data * n_data_vars * p;
         
         double val = deriv_err_array[src_idx];
         
         // Set to 0 if NA or if sigma was not NA
         if (NumericVector::is_na(val) || !sNAIndex[i]) {
           deriv_err_result[dest_idx] = 0.0;
         } else {
           deriv_err_result[dest_idx] = val;
         }
       }
     }
     
     deriv_err_result.attr("dim") = result_dims;
     List result_dimnames = List::create(R_NilValue, data_vars, pars);
     deriv_err_result.attr("dimnames") = result_dimnames;
     
     deriv_err_out = deriv_err_result;
   }
  }
  
  // Handle error model second derivatives (4D array format)
  SEXP deriv2_err_out = R_NilValue;
  
  if (has_err) {
   NumericMatrix err_mat = as<NumericMatrix>(err);
   
   if (err_mat.hasAttribute("deriv2")) {
     NumericVector deriv2_err_array = err_mat.attr("deriv2");
     IntegerVector deriv2_err_dims = deriv2_err_array.attr("dim");
     
     int n_times_d = deriv2_err_dims[0];
     int n_vars = deriv2_err_dims[1];
     int n_pars = deriv2_err_dims[2];
     int n_pars2 = deriv2_err_dims[3];
     
     CharacterVector pars;
     if (deriv2_err_array.hasAttribute("dimnames")) {
       List deriv2_err_dimnames = deriv2_err_array.attr("dimnames");
       if (deriv2_err_dimnames.size() >= 3 && !Rf_isNull(deriv2_err_dimnames[2])) {
         pars = as<CharacterVector>(deriv2_err_dimnames[2]);
       }
     }
     
     if (pars.size() == 0) {
       pars = CharacterVector(n_pars);
       for (int p = 0; p < n_pars; p++) {
         pars[p] = "p" + std::to_string(p + 1);
       }
     }
     
     // Create 4D array for error second derivatives
     CharacterVector data_vars = unique(data_name);
     int n_data_vars = data_vars.size();
     
     NumericVector deriv2_err_result(n_data * n_data_vars * n_pars * n_pars2);
     IntegerVector result_dims = IntegerVector::create(n_data, n_data_vars, n_pars, n_pars2);
     
     // Match indices
     CharacterVector err_colnames = colnames(err_mat);
     IntegerVector err_time_idx(times.size());
     for (int t = 0; t < times.size(); t++) {
       for (int i = 0; i < err_mat.nrow(); i++) {
         if (times[t] == err_mat(i, 0)) {
           err_time_idx[t] = i;
           break;
         }
       }
     }
     
     for (int i = 0; i < n_data; i++) {
       int time_idx = err_time_idx[data_time_idx[i]];
       int var_idx = data_name_idx[i];
       
       for (int p1 = 0; p1 < n_pars; p1++) {
         for (int p2 = 0; p2 < n_pars2; p2++) {
           int src_idx = time_idx + n_times_d * var_idx + 
             n_times_d * n_vars * p1 + 
             n_times_d * n_vars * n_pars * p2;
           int dest_idx = i + n_data * var_idx + 
             n_data * n_data_vars * p1 + 
             n_data * n_data_vars * n_pars * p2;
           
           double val = deriv2_err_array[src_idx];
           
           // Set to 0 if NA or if sigma was not NA
           if (NumericVector::is_na(val) || !sNAIndex[i]) {
             deriv2_err_result[dest_idx] = 0.0;
           } else {
             deriv2_err_result[dest_idx] = val;
           }
         }
       }
     }
     
     deriv2_err_result.attr("dim") = result_dims;
     List result_dimnames = List::create(R_NilValue, data_vars, pars, pars);
     deriv2_err_result.attr("dimnames") = result_dimnames;
     
     deriv2_err_out = deriv2_err_result;
   }
  }
  
  // Apply LLOQ constraint
  for (int i = 0; i < n_data; i++) {
   if (data_value[i] < data_lloq[i]) {
     data_value[i] = data_lloq[i];
   }
  }
  
  // Compute residuals
  LogicalVector is_bloq(n_data);
  NumericVector residual(n_data);
  NumericVector weighted_residual(n_data);
  NumericVector weighted_0(n_data);
  
  for (int i = 0; i < n_data; i++) {
   is_bloq[i] = (data_value[i] <= data_lloq[i]);
   residual[i] = prediction[i] - data_value[i];
   weighted_residual[i] = residual[i] / data_sigma[i];
   weighted_0[i] = prediction[i] / data_sigma[i];
  }
  
  // Create output data frame
  DataFrame result = DataFrame::create(
   _["time"] = data_time,
   _["name"] = data_name,
   _["value"] = data_value,
   _["sigma"] = data_sigma,
   _["lloq"] = data_lloq,
   _["prediction"] = prediction,
   _["residual"] = residual,
   _["weighted.residual"] = weighted_residual,
   _["weighted.0"] = weighted_0,
   _["bloq"] = is_bloq
  );
  
  // Return as list with 3D/4D array attributes
  List output = List::create(
   _["data"] = result,
   _["deriv"] = deriv_out,
   _["deriv2"] = deriv2_out,
   _["deriv.err"] = deriv_err_out,
   _["deriv2.err"] = deriv2_err_out
  );
  
  return output;
}