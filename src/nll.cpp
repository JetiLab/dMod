#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
#include "nll_helpers.hpp"
using namespace Rcpp;

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
//'   with columns corresponding to parameters only (no time/name columns).
//' @param derivs_err NumericMatrix of first derivatives of sigma,
//'   same dimension as derivs (n_data × n_pars) or NULL.
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
 if (n_cols <= 0) {
   stop("'derivs' must have at least one parameter column.");
 }
 
 int n_pars = n_cols; // ALL columns are parameter derivatives
 CharacterVector par_names = colnames(deriv_mat);
 
 // dxdp: [n, n_pars]
 NumericMatrix dxdp(n, n_pars);
 for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n_pars; ++j) {
     dxdp(i, j) = deriv_mat(i, j);
   }
 }
 
 // dsdp: [n, n_pars]
 NumericMatrix dsdp(n, n_pars);
 if (!derivs_err.isNull()) {
   NumericMatrix deriv_err_mat(derivs_err.get());
   int ne_rows = deriv_err_mat.nrow();
   int ne_cols = deriv_err_mat.ncol();
   if (ne_rows != n || ne_cols != n_pars) {
     stop("'derivs_err' must have dimensions [n_data × n_pars].");
   }
   for (int i = 0; i < n; ++i) {
     for (int j = 0; j < n_pars; ++j) {
       dsdp(i, j) = deriv_err_mat(i, j);
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
     double tmp_wr  = inv_s * dx - (wr_i * inv_s) * ds;
     double tmp_w0  = inv_s * dx - (w0_i * inv_s) * ds;
     double tmp_log = inv_s * ds;    // d(log s) / dθ
     
     dwrdp(i, j)  = tmp_wr;
     dw0dp(i, j)  = tmp_w0;
     dlogsdp(i,j) = tmp_log;
   }
 }
 
 // Precompute G(w0) if M4BEAL is active
 std::vector<double> G_w0(n, 0.0);
 if (use_M4BEAL) {
   for (int i = 0; i < n; ++i) {
     G_w0[i] = G_single(w0[i]);
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
       bool ALOQ_part1_flag = ALOQ_part1;
       
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


// -----------------------------------------------------------
// nll_BLOQ_cpp
// -----------------------------------------------------------

//' Non-linear log likelihood for the BLOQ part (C++ implementation)
//'
//' @param nout_bloq DataFrame with BLOQ rows from res()/res_cpp, must contain:
//'   - value             (numeric) : original DV (for LLOQ check in M4)
//'   - weighted.residual (numeric) : wr
//'   - weighted.0        (numeric) : w0
//'   - sigma             (numeric) : s
//' @param derivs_bloq NumericMatrix of first derivatives of prediction,
//'   dimension n_data × n_pars, only parameter columns (no time/name).
//' @param derivs_err_bloq NumericMatrix of first derivatives of sigma,
//'   same dimension as derivs_bloq, or NULL.
//' @param opt_BLOQ Character scalar, one of "M3","M4NM","M4BEAL","M1".
//' @param opt_hessian Named LogicalVector controlling Hessian pieces
//'   ("BLOQ_part1","BLOQ_part2","BLOQ_part3").
//'
//' @return An object of class "objlist" (list with value, gradient, hessian).
//' @export
// [[Rcpp::export]]
Rcpp::List nll_BLOQ_cpp(const Rcpp::DataFrame& nout_bloq,
                        Rcpp::Nullable<Rcpp::NumericMatrix> derivs_bloq = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericMatrix> derivs_err_bloq = R_NilValue,
                        std::string opt_BLOQ = "M3",
                        Rcpp::LogicalVector opt_hessian = Rcpp::LogicalVector()) {
  
  // -----------------------------------------------
  // Extract required vectors
  // -----------------------------------------------
  NumericVector value = nout_bloq["value"];
  NumericVector wr    = nout_bloq["weighted.residual"];
  NumericVector w0    = nout_bloq["weighted.0"];
  NumericVector s     = nout_bloq["sigma"];
  
  int n = wr.size();
  
  // -----------------------------------------------
  // Objective value
  // -----------------------------------------------
  NumericVector objvals(n);
  
  if (opt_BLOQ == "M3") {
    for (int i = 0; i < n; ++i)
      objvals[i] = -2.0 * R::pnorm5(-wr[i], 0,1,1,1);
  }
  
  else if (opt_BLOQ == "M4NM" || opt_BLOQ == "M4BEAL") {
    for (int i = 0; i < n; ++i) {
      double Phi_wr = R::pnorm5(wr[i],0,1,1,0);
      double Phi_w0 = R::pnorm5(w0[i],0,1,1,0);
      double r      = Phi_wr / Phi_w0;
      objvals[i] = (r>=1) ? R_PosInf : -2*std::log(1-r);
    }
    
    // numeric fix
    for (int i=0;i<n;++i) {
      if (!R_finite(objvals[i])) {
        double d= w0[i]-wr[i];
        double ld=std::log(d);
        double intercept = (ld>0?1.8:-1.9*ld+0.9);
        double lin       = (ld>0?0.9:0.5);
        objvals[i] = intercept + lin*w0[i] + 0.95*w0[i]*w0[i];
      }
    }
  }
  
  double obj = 0;
  for (int i=0;i<n;++i) obj+=objvals[i];
  
  // ------------------------------------------------
  // No derivatives → return objlist(value,NULL,NULL)
  // ------------------------------------------------
  if (derivs_bloq.isNull()) {
    Rcpp::List out = Rcpp::List::create(
      _["value"]=obj,
      _["gradient"]=R_NilValue,
      _["hessian"]=R_NilValue
    );
    out.attr("class") = CharacterVector::create("objlist","list");
    return out;
  }
  
  // ------------------------------------------------
  // Extract dxdp, dsdp
  // ------------------------------------------------
  NumericMatrix dxdp(derivs_bloq.get());
  int p = dxdp.ncol();
  
  NumericMatrix dsdp(n,p);
  if (!derivs_err_bloq.isNull()) {
    NumericMatrix D(derivs_err_bloq.get());
    if (D.nrow()!=n || D.ncol()!=p)
      stop("derivs.err.bloq wrong dims");
    for (int i=0;i<n;++i)
      for (int j=0;j<p;++j)
        dsdp(i,j) = D(i,j);
  } else {
    std::fill(dsdp.begin(), dsdp.end(), 0.0);
  }
  
  CharacterVector parnames = colnames(dxdp);
  
  // ------------------------------------------------
  // First-order derivatives dwrdp, dw0dp, dlogsdp
  // ------------------------------------------------
  NumericMatrix dwrdp(n,p), dw0dp(n,p), dlogsdp(n,p);
  
  for (int i=0;i<n;++i) {
    double invs=1.0/s[i];
    for (int j=0;j<p;++j) {
      double dx=dxdp(i,j);
      double ds=dsdp(i,j);
      dwrdp(i,j)= invs*dx - (wr[i]*invs)*ds;
      dw0dp(i,j)= invs*dx - (w0[i]*invs)*ds;
      dlogsdp(i,j)= invs*ds;
    }
  }
  
  // ------------------------------------------------
  // Gradient
  // ------------------------------------------------
  NumericVector grad(p); grad.names()=parnames;
  
  if (opt_BLOQ=="M3") {
    NumericVector Gm(n);
    for (int i=0;i<n;++i) Gm[i]=G_single(-wr[i],-wr[i]);
    
    for (int j=0;j<p;++j) {
      double s_=0;
      for (int i=0;i<n;++i)
        s_+= Gm[i]*dwrdp(i,j);
      grad[j]=2*s_;
    }
  }
  
  else if (opt_BLOQ=="M4NM" || opt_BLOQ=="M4BEAL") {
    
    NumericVector G_wr_w0 = G_vec(wr, w0);
    NumericVector G_wr_wr = G_vec(wr, wr);
    NumericVector G_w0_w0 = G_vec(w0, w0);
    NumericVector G_w0_wr = G_vec(w0, wr);
    NumericVector G0      = G_vec(w0, w0);
    
    NumericVector c1(n),c2(n),c3(n);
    for (int i=0;i<n;++i) {
      c1[i]=2.0/(1.0/G_wr_w0[i] - 1.0/G_wr_wr[i]);
      c2[i]=2.0/(1.0/G_w0_w0[i] - 1.0/G_w0_wr[i]);
      c3[i]=2.0*G0[i];
    }
    
    for (int j=0;j<p;++j) {
      double S1=0,S2=0,S3=0;
      for (int i=0;i<n;++i) {
        S1+=c1[i]*dwrdp(i,j);
        S2+=c2[i]*dw0dp(i,j);
        S3+=c3[i]*dw0dp(i,j);
      }
      grad[j] = S1 - S2 + S3;
    }
  }
  
  // ------------------------------------------------
  // Hessian
  // ------------------------------------------------
  NumericMatrix H(p,p); rownames(H)=parnames; colnames(H)=parnames;
  
  bool P1 = get_flag(opt_hessian,"BLOQ_part1",true);
  bool P2 = get_flag(opt_hessian,"BLOQ_part2",true);
  bool P3 = get_flag(opt_hessian,"BLOQ_part3",true);
  
  if (opt_BLOQ=="M3") {
    
    NumericVector Gm(n);
    for (int i=0;i<n;++i) Gm[i]=G_single(-wr[i],-wr[i]);
    
    // part1
    if (P1) {
      for (int j=0;j<p;++j)
        for (int k=0;k<p;++k) {
          
          double sum_=0;
          for (int i=0;i<n;++i) {
            
            double coef = (-wr[i]*Gm[i] + Gm[i]*Gm[i]);
            sum_ += coef * dwrdp(i,j)*dwrdp(i,k);
          }
          H(j,k) += 2*sum_;
        }
    }
    
    // part2
    if (P2) {
      for (int j=0;j<p;++j)
        for (int k=0;k<p;++k) {
          
          double s1=0,s2=0;
          for (int i=0;i<n;++i) {
            double invs2=1.0/(s[i]*s[i]);
            double coef = Gm[i]*invs2;
            s1+=coef*dxdp(i,j)*dsdp(i,k);
            s2+=coef*dsdp(i,j)*dxdp(i,k);
          }
          H(j,k) -= 2*(s1+s2);
        }
    }
    
    // part3
    if (P3) {
      for (int j=0;j<p;++j)
        for (int k=0;k<p;++k) {
          
          double sum3=0;
          for (int i=0;i<n;++i) {
            double invs2=1.0/(s[i]*s[i]);
            double coef = Gm[i]*(2*(-wr[i]))*invs2;
            sum3 += coef * dsdp(i,j) * dsdp(i,k);
          }
          H(j,k) -= 2*sum3;
        }
    }
  }
  
  // ------------------------------------------------
  // M4 Hessian
  // ------------------------------------------------
  else if (opt_BLOQ == "M4NM" || opt_BLOQ == "M4BEAL") {
    
    // ------------------------------------------------------------------
    // Precompute stable factors
    // ------------------------------------------------------------------
    NumericVector stable_wr = stable_fun(wr, w0, wr);
    NumericVector stable_w0 = stable_fun(w0, w0, wr);
    NumericVector G0        = G_vec(w0, w0);
    
    // Number of parameters
    int p = dxdp.ncol();
    
    // ------------------------------------------------------------------
    // Build coefficient vectors A1..A6
    // ------------------------------------------------------------------
    NumericVector A1(n), A2(n), A3(n), A4(n), A5(n), A6(n);
    
    for (int i = 0; i < n; ++i) {
      A1[i] = -wr[i] * stable_wr[i];
      A2[i] =  stable_wr[i];
      
      A3[i] = -w0[i] * stable_w0[i];
      A4[i] =  stable_w0[i];
      
      A5[i] = -w0[i] * G0[i] - G0[i] * G0[i];
      A6[i] =  G0[i];
    }
    
    // ------------------------------------------------------------------
    // part1 = 2 * ( d(A1) + d2(A2) − d(A3) − d2(A4) )
    // ------------------------------------------------------------------
    NumericMatrix tmp11 = d_dp_sq_cpp(A1, wr, s, dxdp, dsdp);
    NumericMatrix tmp12 = d2_dp2_cpp(A2, wr, s, dxdp, dsdp);
    NumericMatrix tmp21 = d_dp_sq_cpp(A3, w0, s, dxdp, dsdp);
    NumericMatrix tmp22 = d2_dp2_cpp(A4, w0, s, dxdp, dsdp);
    
    NumericMatrix part1(p, p);
    
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) {
        part1(j, k) = 2.0 * (
          tmp11(j, k) +
            tmp12(j, k) -
            tmp21(j, k) -
            tmp22(j, k)
        );
      }
    }
    
    // ------------------------------------------------------------------
    // part2 = -2 * Σ_i (a_ij * a_ik)
    //         with a_i = stable_wr * dwrdp − stable_w0 * dw0dp
    // ------------------------------------------------------------------
    NumericMatrix part2(p, p);
    
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) {
        
        double sum = 0.0;
        
        for (int i = 0; i < n; ++i) {
          double aij = stable_wr[i] * dwrdp(i, j)
          - stable_w0[i] * dw0dp(i, j);
          double aik = stable_wr[i] * dwrdp(i, k)
            - stable_w0[i] * dw0dp(i, k);
          
          sum += aij * aik;
        }
        
        part2(j, k) = -2.0 * sum;
      }
    }
    
    // ------------------------------------------------------------------
    // part3 = 2 * ( d(A5) + d2(A6) )
    // ------------------------------------------------------------------
    NumericMatrix tmp31 = d_dp_sq_cpp(A5, w0, s, dxdp, dsdp);
    NumericMatrix tmp32 = d2_dp2_cpp(A6, w0, s, dxdp, dsdp);
    
    NumericMatrix part3(p, p);
    
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) {
        part3(j, k) = 2.0 * (
          tmp31(j, k) +
            tmp32(j, k)
        );
      }
    }
    
    // ------------------------------------------------------------------
    // Accumulate contributions into Hessian H
    // ------------------------------------------------------------------
    if (P1) {
      for (int j = 0; j < p; ++j)
        for (int k = 0; k < p; ++k)
          H(j, k) += part1(j, k);
    }
    
    if (P2) {
      for (int j = 0; j < p; ++j)
        for (int k = 0; k < p; ++k)
          H(j, k) += part2(j, k);
    }
    
    if (P3) {
      for (int j = 0; j < p; ++j)
        for (int k = 0; k < p; ++k)
          H(j, k) += part3(j, k);
    }
  }
  
  // ------------------------------------------------
  // Build objlist
  // ------------------------------------------------
  Rcpp::List out = Rcpp::List::create(
    _["value"]    = obj,
    _["gradient"] = grad,
    _["hessian"]  = H
  );
  out.attr("class") = CharacterVector::create("objlist","list");
  return out;
}
