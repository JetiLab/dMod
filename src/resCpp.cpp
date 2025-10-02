#include <Rcpp.h>
using namespace Rcpp;

//' Residual calculation with optional BLOQ handling
 //'
 //' @description C++ implementation of `res()` with support for BLOQ methods (M3, M4NM, M4BEAL, M1).
 //' Combines data, model predictions, residuals, weights, and optional derivatives.
 //'
 //' @param data A data.frame with columns: time, name, value, prediction, sigma, bloq.
 //' @param out Numeric matrix of model predictions.
 //' @param deriv (optional) Matrix of derivatives of predictions wrt. structural parameters.
 //' @param deriv_err (optional) Matrix of derivatives wrt. error model parameters.
 //' @param optBLOQ Character string, one of "none", "M3", "M4NM", "M4BEAL", "M1".
 //'
 //' @return A list with elements:
 //' \itemize{
 //'   \item `data` = data.frame with residuals and weighted residuals
 //'   \item `deriv` = combined Jacobian (if provided)
 //' }
 //'
 //' @export
 // [[Rcpp::export]]
 List resCpp(DataFrame data,
             NumericMatrix out,
             Nullable<NumericMatrix> deriv = R_NilValue,
             Nullable<NumericMatrix> deriv_err = R_NilValue,
             std::string optBLOQ = "none") {
   
   // Extract inputs
   NumericVector value = data["value"];
   NumericVector prediction = data["prediction"];
   NumericVector sigma = data["sigma"];
   LogicalVector bloq = data.containsElementNamed("bloq") ? 
   as<LogicalVector>(data["bloq"]) : LogicalVector(value.size(), false);
   
   int n = value.size();
   NumericVector residual(n);
   NumericVector weighted_residual(n);
   
   // BLOQ handling
   for (int i = 0; i < n; i++) {
     double v = value[i];
     double p = prediction[i];
     double s = sigma[i];
     
     if (!bloq[i]) {
       residual[i] = v - p;
       weighted_residual[i] = residual[i] / s;
     } else {
       if (optBLOQ == "M3") {
         residual[i] = R::pnorm(v, p, s, 1, 1);
         weighted_residual[i] = residual[i];
       } else if (optBLOQ == "M4NM") {
         residual[i] = 0.0;
         weighted_residual[i] = 0.0;
       } else if (optBLOQ == "M4BEAL") {
         residual[i] = v - p;
         weighted_residual[i] = residual[i] / s;
       } else if (optBLOQ == "M1") {
         residual[i] = v / 2.0 - p;
         weighted_residual[i] = residual[i] / s;
       } else {
         residual[i] = v - p;
         weighted_residual[i] = residual[i] / s;
       }
     }
   }
   
   // Add residuals to data
   data["residual"] = residual;
   data["weighted.residual"] = weighted_residual;
   
   // Derivatives
   List result;
   result["data"] = data;
   
   if (deriv.isNotNull() || deriv_err.isNotNull()) {
     NumericMatrix J_pred, J_err;
     
     if (deriv.isNotNull()) J_pred = as<NumericMatrix>(deriv);
     if (deriv_err.isNotNull()) J_err = as<NumericMatrix>(deriv_err);
     
     // Combine Jacobians: here just row-bind or cbind depending on convention
     if (deriv.isNotNull() && deriv_err.isNotNull()) {
       int nr = J_pred.nrow();
       int nc1 = J_pred.ncol();
       int nc2 = J_err.ncol();
       NumericMatrix J(nr, nc1 + nc2);
       
       for (int i = 0; i < nr; i++) {
         for (int j = 0; j < nc1; j++) J(i, j) = J_pred(i, j);
         for (int j = 0; j < nc2; j++) J(i, nc1 + j) = J_err(i, j);
       }
       result["deriv"] = J;
     } else if (deriv.isNotNull()) {
       result["deriv"] = J_pred;
     } else if (deriv_err.isNotNull()) {
       result["deriv"] = J_err;
     }
   }
   
   return result;
 }
 