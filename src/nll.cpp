#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// =======================
// Helper: extract named logical flag
// =======================

// Safely extract a named logical flag from a LogicalVector.
// If not found, returns default_value.
bool get_flag(const LogicalVector& lv,
              const std::string& name,
              bool default_value = true) {
  if (lv.size() == 0) return default_value;
  CharacterVector nms = lv.names();
  if (nms.size() != lv.size()) return default_value;
  
  for (int i = 0; i < lv.size(); ++i) {
    if (TYPEOF(nms[i]) != STRSXP) continue;
    std::string s = as<std::string>(nms[i]);
    if (s == name) {
      if (LogicalVector::is_na(lv[i])) return false;
      return lv[i];
    }
  }
  return default_value;
}

// =======================
// Helper: G_by_Phi(w) = φ(w)/Φ(w)
// =======================

inline double G_by_Phi_single(double w) {
  double log_phi = R::dnorm4(w, 0.0, 1.0, 1);    // log density
  double log_Phi = R::pnorm5(w, 0.0, 1.0, 1, 1); // log CDF
  return std::exp(log_phi - log_Phi);
}

// =======================
// nll_ALOQ_cpp
// =======================

//' Non-linear log likelihood for the ALOQ part (C++ implementation)
//'
//' @param nout DataFrame from res/res_cpp, must contain:
//'   - weighted.residual (numeric) : wr (unscaled)
//'   - weighted.0       (numeric) : w0 (unscaled)
//'   - sigma            (numeric) : s
//' @param derivs NumericMatrix of first derivatives of prediction (NOT wr),
//'   with columns: time, name, then parameter columns.
//' @param derivs_err NumericMatrix of first derivatives of sigma,
//'   same structure as derivs (time, name, parameters) or NULL.
//' @param opt_BLOQ Character scalar, e.g. "M3", "M4NM", "M4BEAL", "M1".
//'   Only "M4BEAL" changes the ALOQ contribution (extra Phi(w0) term).
//' @param opt_hessian Named LogicalVector controlling Hessian components
//'   ("ALOQ_part1", "ALOQ_part2", "ALOQ_part3") when 2nd derivatives are NOT supplied.
//' @param bessel_correction Numeric scalar. If != 1, weighted residuals
//'   and w0 are multiplied by this factor (Bessel correction) for the
//'   optimization objective. The TRUE ML (-2logL without Bessel) is stored
//'   in the "neg2ll" attribute.
//' @param deriv2 Optional 3D array (NumericVector with dim attribute)
//'   of second derivatives of prediction: dim = c(n_data, n_pars, n_pars).
//' @param deriv2_err Optional 3D array of second derivatives of sigma,
//'   same dims as deriv2.
//'
//' @details
//' If deriv2 is NULL, the Hessian is computed using the same first-order
//' approximation as in the original R implementation, controlled by opt_hessian,
//' including the approximate M4BEAL Hessian block.
//'
//' If deriv2 is NOT NULL, a full exact Hessian for the Gaussian ALOQ part
//' and for the M4BEAL ALOQ contribution is computed using the 2nd derivatives.
//' In that case, opt_hessian is ignored.
//'
//' @return List with elements:
//'   - value    : corrected negative log-likelihood (for optimization)
//'   - gradient : numeric vector (or NULL) w.r.t parameters
//'   - hessian  : numeric matrix (or NULL)
//'   Attributes:
//'   - "chisquare"       : TRUE chisquare (no Bessel correction)
//'   - "neg2ll"          : TRUE -2*log(L) (no Bessel, but with M4BEAL term if used)
//'   - "besselcorrected" : logical
//'
//' @export
// [[Rcpp::export]]
Rcpp::List nll_ALOQ_cpp(const Rcpp::DataFrame& nout,
                       Rcpp::Nullable<Rcpp::NumericMatrix> derivs       = R_NilValue,
                       Rcpp::Nullable<Rcpp::NumericMatrix> derivs_err   = R_NilValue,
                       std::string opt_BLOQ                             = "M3",
                       Rcpp::LogicalVector opt_hessian                  = Rcpp::LogicalVector(),
                       double bessel_correction                         = 1.0,
                       Rcpp::Nullable<Rcpp::NumericVector> deriv2       = R_NilValue,
                       Rcpp::Nullable<Rcpp::NumericVector> deriv2_err   = R_NilValue) {
 
 // -------------------------
 // 1. Extract basic quantities from nout
 // -------------------------
 NumericVector wr = nout["weighted.residual"]; // unscaled
 NumericVector w0 = nout["weighted.0"];        // unscaled
 NumericVector s  = nout["sigma"];
 
 int n = wr.size();
 if (w0.size() != n || s.size() != n) {
   stop("nout columns 'weighted.residual', 'weighted.0', 'sigma' must have same length.");
 }
 
 const double two_pi = 2.0 * M_PI;
 
 // -------------------------
 // 2. TRUE ML quantities (no Bessel correction)
 // -------------------------
 double chisquare_ml = 0.0;
 double log_term_sum = 0.0;
 
 for (int i = 0; i < n; ++i) {
   double wi = wr[i];
   double si = s[i];
   chisquare_ml += wi * wi;
   log_term_sum += std::log(two_pi * si * si);
 }
 
 double neg2ll_ml = chisquare_ml + log_term_sum;
 
 // -------------------------
 // 3. Bessel correction: used only for optimization objective
 // -------------------------
 bool use_bessel = (std::fabs(bessel_correction - 1.0) > 1e-14);
 double b = bessel_correction;
 
 if (use_bessel) {
   for (int i = 0; i < n; ++i) {
     wr[i] *= b;   // wr' = b * wr0
     w0[i] *= b;   // w0' = b * w0_unscaled
   }
 }
 
 // -------------------------
 // 4. Corrected objective value (for optimization)
 // -------------------------
 double chisquare = 0.0;
 for (int i = 0; i < n; ++i) {
   chisquare += wr[i] * wr[i];
 }
 double obj = chisquare + log_term_sum;
 
 bool use_M4BEAL = (opt_BLOQ == "M4BEAL");
 
 if (use_M4BEAL) {
   double bloq_term = 0.0;
   for (int i = 0; i < n; ++i) {
     // 2 * log( Phi( w0_i ) ), w0_i already Bessel-scaled if use_bessel
     bloq_term += 2.0 * R::pnorm5(w0[i], 0.0, 1.0, 1, 1);
   }
   obj       += bloq_term;
   neg2ll_ml += bloq_term; // TRUE ML includes this extra term too
 }
 
 // -------------------------
 // 5. If no first derivatives are provided -> no grad/hessian
 // -------------------------
 NumericVector grad;
 NumericMatrix hessian;
 
 if (derivs.isNull()) {
   Rcpp::List out = Rcpp::List::create(
     _["value"]    = obj,
     _["gradient"] = R_NilValue,
     _["hessian"]  = R_NilValue
   );
   out.attr("chisquare")       = chisquare_ml;
   out.attr("neg2ll")          = neg2ll_ml;
   out.attr("besselcorrected") = use_bessel;
   return out;
 }
 
 // -------------------------
 // 6. Extract first derivatives (dxdp, dsdp)
 // -------------------------
 NumericMatrix deriv_mat(derivs.get());
 int n_rows = deriv_mat.nrow();
 int n_cols = deriv_mat.ncol();
 
 if (n_rows != n) {
   stop("Number of rows in 'derivs' must match number of residuals.");
 }
 if (n_cols <= 2) {
   stop("'derivs' must have at least 3 columns: time, name, parameter(s).");
 }
 
 int n_pars = n_cols - 2; // drop time and name
 CharacterVector deriv_colnames = colnames(deriv_mat);
 CharacterVector par_names(n_pars);
 for (int j = 0; j < n_pars; ++j) {
   par_names[j] = deriv_colnames[j + 2];
 }
 
 // dxdp: [n, n_pars]
 NumericMatrix dxdp(n, n_pars);
 for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n_pars; ++j) {
     dxdp(i, j) = deriv_mat(i, j + 2);
   }
 }
 
 // dsdp: [n, n_pars]
 NumericMatrix dsdp(n, n_pars);
 if (!derivs_err.isNull()) {
   NumericMatrix deriv_err_mat(derivs_err.get());
   int ne_rows = deriv_err_mat.nrow();
   int ne_cols = deriv_err_mat.ncol();
   if (ne_rows != n || ne_cols != n_cols) {
     stop("'derivs_err' must have the same dimensions as 'derivs'.");
   }
   for (int i = 0; i < n; ++i) {
     for (int j = 0; j < n_pars; ++j) {
       dsdp(i, j) = deriv_err_mat(i, j + 2);
     }
   }
 } else {
   std::fill(dsdp.begin(), dsdp.end(), 0.0);
 }
 
 // Check if there is any non-zero dsdp at all
 bool has_dsdp = false;
 for (int i = 0; i < n && !has_dsdp; ++i) {
   for (int j = 0; j < n_pars; ++j) {
     if (dsdp(i,j) != 0.0) { has_dsdp = true; break; }
   }
 }
 
 // -------------------------
 // 7. Compute first-order derivative matrices: dwrdp, dw0dp, dlogsdp
 // -------------------------
 NumericMatrix dwrdp(n, n_pars);
 NumericMatrix dw0dp(n, n_pars);
 NumericMatrix dlogsdp(n, n_pars);
 
 for (int i = 0; i < n; ++i) {
   double si    = s[i];
   double inv_s = 1.0 / si;
   double wr_i  = wr[i];
   double w0_i  = w0[i];
   
   for (int j = 0; j < n_pars; ++j) {
     double dx = dxdp(i, j);
     double ds = dsdp(i, j);
     
     // NOTE: We treat wr and w0 as already Bessel-scaled values.
     // These formulas match your R implementation structure.
     double tmp_wr  = inv_s * dx - (wr_i * inv_s) * ds;  // = b/s*dxdp - wr'/s*dsdp
     double tmp_w0  = inv_s * dx - (w0_i * inv_s) * ds;  // = b/s*dxdp - w0'/s*dsdp
     double tmp_log = inv_s * ds;                        // d(log s) / dθ
     
     dwrdp(i, j)  = tmp_wr;
     dw0dp(i, j)  = tmp_w0;
     dlogsdp(i,j) = tmp_log;
   }
 }
 
 // Precompute G(w0) if M4BEAL is active
 std::vector<double> G_w0(n, 0.0);
 if (use_M4BEAL) {
   for (int i = 0; i < n; ++i) {
     G_w0[i] = G_by_Phi_single(w0[i]);  // w0 already Bessel-scaled if needed
   }
 }
 
 // -------------------------
 // 8. Gradient (same form in both approx and exact Hessian modes)
 // -------------------------
 grad = NumericVector(n_pars);
 grad.names() = par_names;
 
 // Gaussian part: 2 * wr^T * dwrdp
 for (int j = 0; j < n_pars; ++j) {
   double sum_wr = 0.0;
   for (int i = 0; i < n; ++i) {
     sum_wr += wr[i] * dwrdp(i, j);
   }
   grad[j] = 2.0 * sum_wr;
 }
 
 // Sigma/log part: 2 * sum(dlogsdp_i)
 for (int j = 0; j < n_pars; ++j) {
   double sum_log = 0.0;
   for (int i = 0; i < n; ++i) {
     sum_log += dlogsdp(i, j);
   }
   grad[j] += 2.0 * sum_log;
 }
 
 // M4BEAL gradient: 2 * G(w0)^T * dw0dp
 if (use_M4BEAL) {
   for (int j = 0; j < n_pars; ++j) {
     double sum_g = 0.0;
     for (int i = 0; i < n; ++i) {
       sum_g += G_w0[i] * dw0dp(i, j);
     }
     grad[j] += 2.0 * sum_g;
   }
 }
 
 // -------------------------
 // 9. Hessian
 // -------------------------
 hessian = NumericMatrix(n_pars, n_pars);
 rownames(hessian) = par_names;
 colnames(hessian) = par_names;
 
 // ===========================================================
 // CASE A: NO second derivatives (deriv2 == NULL)
 // → use approximate Hessian exactly as in R implementation
 // ===========================================================
 if (deriv2.isNull()) {
   bool ALOQ_part1 = get_flag(opt_hessian, "ALOQ_part1", true);
   bool ALOQ_part2 = get_flag(opt_hessian, "ALOQ_part2", true);
   bool ALOQ_part3 = get_flag(opt_hessian, "ALOQ_part3", true);
   
   // Base Gaussian term: 2 * t(dwrdp) %*% dwrdp
   for (int j = 0; j < n_pars; ++j) {
     for (int k = 0; k < n_pars; ++k) {
       double sum_dw = 0.0;
       for (int i = 0; i < n; ++i) {
         sum_dw += dwrdp(i, j) * dwrdp(i, k);
       }
       hessian(j, k) = 2.0 * sum_dw;
     }
   }
   
   if (has_dsdp) {
     // Precompute useful scalars
     std::vector<double> inv_s2(n), a(n), bvec(n);
     for (int i = 0; i < n; ++i) {
       double si   = s[i];
       double si2  = si * si;
       inv_s2[i]   = 1.0 / si2;
       a[i]        = -wr[i] * inv_s2[i];
       bvec[i]     =  2.0 * wr[i] * wr[i] * inv_s2[i];
     }
     
     // ALOQ_part1: 2 * (t(-wr/s^2 * dxdp) %*% dsdp + t(-wr/s^2 * dsdp) %*% dxdp)
     if (ALOQ_part1) {
       for (int j = 0; j < n_pars; ++j) {
         for (int k = 0; k < n_pars; ++k) {
           double sum1 = 0.0, sum2 = 0.0;
           for (int i = 0; i < n; ++i) {
             double ai  = a[i];
             double dxj = dxdp(i,j);
             double dxk = dxdp(i,k);
             double dsj = dsdp(i,j);
             double dsk = dsdp(i,k);
             sum1 += ai * dxj * dsk;
             sum2 += ai * dsj * dxk;
           }
           hessian(j,k) += 2.0 * (sum1 + sum2);
         }
       }
     }
     
     // ALOQ_part2: 2 * t(2 * wr^2 / s^2 * dsdp) %*% dsdp
     if (ALOQ_part2) {
       for (int j = 0; j < n_pars; ++j) {
         for (int k = 0; k < n_pars; ++k) {
           double sum_b = 0.0;
           for (int i = 0; i < n; ++i) {
             sum_b += bvec[i] * dsdp(i,j) * dsdp(i,k);
           }
           hessian(j,k) += 2.0 * sum_b;
         }
       }
     }
     
     // ALOQ_part3: - 2 * t(dlogsdp) %*% dlogsdp
     if (ALOQ_part3) {
       for (int j = 0; j < n_pars; ++j) {
         for (int k = 0; k < n_pars; ++k) {
           double sum_log = 0.0;
           for (int i = 0; i < n; ++i) {
             sum_log += dlogsdp(i,j) * dlogsdp(i,k);
           }
           hessian(j,k) -= 2.0 * sum_log;
         }
       }
     }
     
     // ---- Approximate extra M4BEAL Hessian (same algebra as R) ----
     if (use_M4BEAL) {
       // Only needed once
       bool ALOQ_part1_flag = ALOQ_part1; // reuse
       
       // M4BEAL Term 1:
       //   2 * t( (-w0 * G - G^2) * dw0dp ) %*% dw0dp
       for (int j = 0; j < n_pars; ++j) {
         for (int k = 0; k < n_pars; ++k) {
           double sum1 = 0.0;
           for (int i = 0; i < n; ++i) {
             double Gi    = G_w0[i];
             double w0i   = w0[i];
             double coef  = (-w0i * Gi - Gi * Gi);
             sum1 += coef * dw0dp(i,j) * dw0dp(i,k);
           }
           hessian(j,k) += 2.0 * sum1;
         }
       }
       
       // M4BEAL Term 2:
       //   2 * t( G(w0)*(-1)/s^2 * dxdp ) %*% dsdp
       // + 2 * t( G(w0)*(-1)/s^2 * dsdp ) %*% dxdp
       std::vector<double> inv_s2_2(n);
       for (int i = 0; i < n; ++i)
         inv_s2_2[i] = 1.0 / (s[i] * s[i]);
       
       for (int j = 0; j < n_pars; ++j) {
         for (int k = 0; k < n_pars; ++k) {
           double sum_dxds = 0.0;
           double sum_dsdx = 0.0;
           for (int i = 0; i < n; ++i) {
             double Gi    = G_w0[i];
             double coef  = -Gi * inv_s2_2[i];
             sum_dxds += coef * dxdp(i,j) * dsdp(i,k);
             sum_dsdx += coef * dsdp(i,j) * dxdp(i,k);
           }
           hessian(j,k) += 2.0 * (sum_dxds + sum_dsdx);
         }
       }
       
       // M4BEAL Term 3 (only if ALOQ_part1 = TRUE):
       //   2 * t( 2 * G(w0) * w0 / s^2 * dsdp ) %*% dsdp
       if (ALOQ_part1_flag) {
         for (int j = 0; j < n_pars; ++j) {
           for (int k = 0; k < n_pars; ++k) {
             double sum3 = 0.0;
             for (int i = 0; i < n; ++i) {
               double Gi   = G_w0[i];
               double coef = 2.0 * Gi * w0[i] * inv_s2_2[i];
               sum3 += coef * dsdp(i,j) * dsdp(i,k);
             }
             hessian(j,k) += 2.0 * sum3;
           }
         }
       }
     } // use_M4BEAL
   }   // has_dsdp
   
   // ===========================================================
   // CASE B: 2nd derivatives provided (deriv2 != NULL)
   // → compute exact Hessian using d2x/dp2 and d2s/dp2
   //    (for both Gaussian ALOQ part and M4BEAL term)
   // ===========================================================
 } else {
   NumericVector d2x_vec(deriv2.get());
   IntegerVector d2x_dim = d2x_vec.attr("dim");
   if (d2x_dim.size() != 3) {
     stop("'deriv2' must be a 3D array with dim = c(n_data, n_pars, n_pars).");
   }
   int n_d2  = d2x_dim[0];
   int np_d2 = d2x_dim[1];
   int np_d2b= d2x_dim[2];
   if (n_d2 != n || np_d2 != n_pars || np_d2b != n_pars) {
     stop("'deriv2' dimensions must be (n_data, n_pars, n_pars).");
   }
   
   NumericVector d2s_vec;
   bool has_d2s = false;
   if (!deriv2_err.isNull()) {
     d2s_vec = NumericVector(deriv2_err.get());
     IntegerVector d2s_dim = d2s_vec.attr("dim");
     if (d2s_dim.size() != 3 ||
         d2s_dim[0] != n || d2s_dim[1] != n_pars || d2s_dim[2] != n_pars) {
       stop("'deriv2_err' must be a 3D array with dim = c(n_data, n_pars, n_pars).");
     }
     has_d2s = true;
   }
   
   // Accessors for 3D arrays: index (i,j,k) in column-major
   auto d2x = [&](int i, int j, int k) -> double {
     return d2x_vec[i + n * (j + n_pars * k)];
   };
   auto d2s = [&](int i, int j, int k) -> double {
     if (!has_d2s) return 0.0;
     return d2s_vec[i + n * (j + n_pars * k)];
   };
   
   // Exact Hessian: sum over data points of per-observation Hessians
   for (int j = 0; j < n_pars; ++j) {
     for (int k = 0; k < n_pars; ++k) {
       double Hjk = 0.0;
       
       for (int i = 0; i < n; ++i) {
         double si    = s[i];
         double si2   = si * si;
         double inv_s = 1.0 / si;
         double inv_s2= 1.0 / si2;
         double wr_i  = wr[i];
         double w0_i  = w0[i];
         
         // first-order pieces at (i,j,k)
         double dx_j = dxdp(i,j), dx_k = dxdp(i,k);
         double ds_j = dsdp(i,j), ds_k = dsdp(i,k);
         
         double dwr_j = dwrdp(i,j), dwr_k = dwrdp(i,k);
         
         double d2x_jk = d2x(i,j,k);
         double d2s_jk = d2s(i,j,k);
         
         // Second derivative of wr' (Bessel-scaled residual):
         // wr' = b * (x/s - y/s)
         // d2wr' = (b/s)*d2x - (b/s^2)*(dx_j*ds_k + dx_k*ds_j)
         //         + 2*wr_i/s^2 * ds_j*ds_k - (wr_i/s)*d2s_jk
         double d2wr_jk =
           (b * inv_s)     * d2x_jk
         - (b * inv_s2)    * (dx_j * ds_k + dx_k * ds_j)
           + (2.0 * wr_i * inv_s2) * ds_j * ds_k
           - (wr_i * inv_s) * d2s_jk;
           
           // Second derivative of log s:
           // d2log(s)/dθj dθk = (1/s)*d2s_jk - (1/s^2)*ds_j*ds_k
           double d2logs_jk = inv_s * d2s_jk - inv_s2 * ds_j * ds_k;
           
           // Gaussian ALOQ contribution:
           // H_i(j,k) = 2( dwr_j*dwr_k + wr_i*d2wr_jk ) + 2*d2logs_jk
           double H_i = 2.0 * (dwr_j * dwr_k + wr_i * d2wr_jk) + 2.0 * d2logs_jk;
           
           // M4BEAL extra ALOQ contribution:
           // term per i, j, k:
           //   2 [ (-w0_i G_i - G_i^2) * dw0_j * dw0_k + G_i * d2w0_jk ]
           if (use_M4BEAL) {
             double G_i = G_w0[i];
             
             double dw0_j = dw0dp(i,j);
             double dw0_k = dw0dp(i,k);
             
             // Second derivative of w0' (Bessel-scaled prediction/sigma):
             // w0' = b * x/s
             // d2w0' = (b/s)*d2x - (b/s^2)*(dx_j*ds_k + dx_k*ds_j)
             //         + 2*w0_i/s^2 * ds_j*ds_k - (w0_i/s)*d2s_jk
             double d2w0_jk =
               (b * inv_s)     * d2x_jk
             - (b * inv_s2)    * (dx_j * ds_k + dx_k * ds_j)
               + (2.0 * w0_i * inv_s2) * ds_j * ds_k
               - (w0_i * inv_s) * d2s_jk;
               
               double term_M4 =
               2.0 * ( (-w0_i * G_i - G_i * G_i) * dw0_j * dw0_k
                         + G_i * d2w0_jk );
               
               H_i += term_M4;
           }
           
           Hjk += H_i;
       } // i
       
       hessian(j, k) = Hjk;
     } // k
   }   // j
 }     // CASE B exact
 
 // -------------------------
 // 10. Build output list
 // -------------------------
 Rcpp::List out = Rcpp::List::create(
   _["value"] = obj
 );
 
 if (grad.size() > 0)
   out["gradient"] = grad;
 else
   out["gradient"] = R_NilValue;
 
 if (hessian.nrow() > 0)
   out["hessian"] = hessian;
 else
   out["hessian"] = R_NilValue;
 
 
 out.attr("chisquare")       = chisquare_ml;  // TRUE chi-square (no Bessel)
 out.attr("neg2ll")          = neg2ll_ml;    // TRUE -2logL (no Bessel, but incl. M4BEAL if used)
 out.attr("besselcorrected") = use_bessel;
 
 return out;
}
