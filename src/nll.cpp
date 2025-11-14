#include "nll.h"
#include "nll_helpers.hpp"
#include <vector>
#include <string>
#include <cmath>

// =============================================================================
// Helper functions to extract data from R objects
// =============================================================================

static double* get_numeric_column(SEXP df, const char* name, int* n_out) {
  SEXP names = Rf_getAttrib(df, R_NamesSymbol);
  int ncols = Rf_length(df);
  
  for (int i = 0; i < ncols; ++i) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      SEXP col = VECTOR_ELT(df, i);
      *n_out = Rf_length(col);
      return REAL(col);
    }
  }
  
  Rf_error("Column '%s' not found in data.frame", name);
  return nullptr;
}

static bool get_named_flag(SEXP opt_hessian, const char* name, bool default_value) {
  if (Rf_isNull(opt_hessian) || Rf_length(opt_hessian) == 0)
    return default_value;
  
  SEXP names = Rf_getAttrib(opt_hessian, R_NamesSymbol);
  if (Rf_isNull(names))
    return default_value;
  
  int n = Rf_length(opt_hessian);
  int* lgl = LOGICAL(opt_hessian);
  
  for (int i = 0; i < n; ++i) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      if (lgl[i] == NA_LOGICAL)
        return false;
      return lgl[i];
    }
  }
  
  return default_value;
}

// =============================================================================
// C_nll_ALOQ
// =============================================================================

extern "C" SEXP C_nll_ALOQ(SEXP nout, SEXP derivs, SEXP derivs_err,
                          SEXP opt_BLOQ, SEXP opt_hessian, SEXP bessel_correction,
                          SEXP deriv2, SEXP deriv2_err) {
  
  // -------------------------
  // 1. Extract basic quantities from nout
  // -------------------------
  int n;
  double* wr = get_numeric_column(nout, "weighted.residual", &n);
  
  int n_check;
  double* w0 = get_numeric_column(nout, "weighted.0", &n_check);
  if (n_check != n) Rf_error("Column lengths mismatch");
  
  double* s = get_numeric_column(nout, "sigma", &n_check);
  if (n_check != n) Rf_error("Column lengths mismatch");
  
  const double two_pi = 2.0 * M_PI;
  
  // Make copies for potential Bessel correction
  std::vector<double> wr_vec(wr, wr + n);
  std::vector<double> w0_vec(w0, w0 + n);
  
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
  // 3. Bessel correction
  // -------------------------
  double b = Rf_asReal(bessel_correction);
  bool use_bessel = (std::fabs(b - 1.0) > 1e-14);
  
  if (use_bessel) {
    for (int i = 0; i < n; ++i) {
      wr_vec[i] *= b;
      w0_vec[i] *= b;
    }
  }
  
  // -------------------------
  // 4. Corrected objective value
  // -------------------------
  double chisquare = 0.0;
  for (int i = 0; i < n; ++i) {
    chisquare += wr_vec[i] * wr_vec[i];
  }
  double obj = chisquare + log_term_sum;
  
  const char* opt_BLOQ_str = CHAR(STRING_ELT(opt_BLOQ, 0));
  bool use_M4BEAL = (strcmp(opt_BLOQ_str, "M4BEAL") == 0);
  
  if (use_M4BEAL) {
    double bloq_term = 0.0;
    for (int i = 0; i < n; ++i) {
      bloq_term += 2.0 * R::pnorm5(w0_vec[i], 0.0, 1.0, 1, 1);
    }
    obj += bloq_term;
    neg2ll_ml += bloq_term;
  }
  
  // -------------------------
  // 5. If no derivatives -> return early
  // -------------------------
  if (Rf_isNull(derivs)) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    
    SET_STRING_ELT(names, 0, Rf_mkChar("value"));
    SET_STRING_ELT(names, 1, Rf_mkChar("gradient"));
    SET_STRING_ELT(names, 2, Rf_mkChar("hessian"));
    
    SEXP value = PROTECT(Rf_ScalarReal(obj));
    SET_VECTOR_ELT(out, 0, value);
    SET_VECTOR_ELT(out, 1, R_NilValue);
    SET_VECTOR_ELT(out, 2, R_NilValue);
    
    Rf_setAttrib(out, R_NamesSymbol, names);
    
    // Attributes
    SEXP attr_chisq = PROTECT(Rf_ScalarReal(chisquare_ml));
    Rf_setAttrib(out, Rf_install("chisquare"), attr_chisq);
    
    SEXP attr_neg2ll = PROTECT(Rf_ScalarReal(neg2ll_ml));
    Rf_setAttrib(out, Rf_install("neg2ll"), attr_neg2ll);
    
    SEXP attr_bessel = PROTECT(Rf_ScalarLogical(use_bessel));
    Rf_setAttrib(out, Rf_install("besselcorrected"), attr_bessel);
    
    UNPROTECT(6);
    return out;
  }
  
  // -------------------------
  // 6. Extract first derivatives
  // -------------------------
  int n_rows = Rf_nrows(derivs);
  int n_pars = Rf_ncols(derivs);
  
  if (n_rows != n) Rf_error("Number of rows in 'derivs' must match number of residuals");
  if (n_pars <= 0) Rf_error("'derivs' must have at least one parameter column");
  
  double* dxdp_data = REAL(derivs);
  
  // Get parameter names
  SEXP dimnames = Rf_getAttrib(derivs, R_DimNamesSymbol);
  SEXP par_names = R_NilValue;
  if (!Rf_isNull(dimnames) && Rf_length(dimnames) >= 2) {
    par_names = VECTOR_ELT(dimnames, 1);
  }
  
  // dsdp
  std::vector<double> dsdp_data(n * n_pars, 0.0);
  bool has_dsdp = false;
  
  if (!Rf_isNull(derivs_err)) {
    int ne_rows = Rf_nrows(derivs_err);
    int ne_cols = Rf_ncols(derivs_err);
    if (ne_rows != n || ne_cols != n_pars)
      Rf_error("'derivs_err' must have dimensions [n_data × n_pars]");
    
    double* derivs_err_data = REAL(derivs_err);
    std::copy(derivs_err_data, derivs_err_data + n * n_pars, dsdp_data.begin());
    
    for (int i = 0; i < n * n_pars; ++i) {
      if (dsdp_data[i] != 0.0) {
        has_dsdp = true;
        break;
      }
    }
  }
  
  // -------------------------
  // 7. Compute first-order derivatives: dwrdp, dw0dp, dlogsdp
  // -------------------------
  std::vector<double> dwrdp(n * n_pars);
  std::vector<double> dw0dp(n * n_pars);
  std::vector<double> dlogsdp(n * n_pars);
  
  for (int i = 0; i < n; ++i) {
    double si = s[i];
    double inv_s = 1.0 / si;
    double wr_i = wr_vec[i];
    double w0_i = w0_vec[i];
    
    for (int j = 0; j < n_pars; ++j) {
      int idx = i + j * n;  // Column-major indexing
      
      double dx = dxdp_data[idx];
      double ds = dsdp_data[idx];
      
      dwrdp[idx] = inv_s * dx - (wr_i * inv_s) * ds;
      dw0dp[idx] = inv_s * dx - (w0_i * inv_s) * ds;
      dlogsdp[idx] = inv_s * ds;
    }
  }
  
  // Precompute G(w0) if M4BEAL
  std::vector<double> G_w0(n, 0.0);
  if (use_M4BEAL) {
    for (int i = 0; i < n; ++i) {
      G_w0[i] = G_single(w0_vec[i]);
    }
  }
  
  // -------------------------
  // 8. Gradient
  // -------------------------
  SEXP grad = PROTECT(Rf_allocVector(REALSXP, n_pars));
  double* grad_data = REAL(grad);
  
  if (!Rf_isNull(par_names)) {
    Rf_setAttrib(grad, R_NamesSymbol, par_names);
  }
  
  // Gaussian part: 2 * wr^T * dwrdp
  for (int j = 0; j < n_pars; ++j) {
    double sum_wr = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_wr += wr_vec[i] * dwrdp[i + j * n];
    }
    grad_data[j] = 2.0 * sum_wr;
  }
  
  // Sigma/log part: 2 * sum(dlogsdp)
  for (int j = 0; j < n_pars; ++j) {
    double sum_log = 0.0;
    for (int i = 0; i < n; ++i) {
      sum_log += dlogsdp[i + j * n];
    }
    grad_data[j] += 2.0 * sum_log;
  }
  
  // M4BEAL gradient
  if (use_M4BEAL) {
    for (int j = 0; j < n_pars; ++j) {
      double sum_g = 0.0;
      for (int i = 0; i < n; ++i) {
        sum_g += G_w0[i] * dw0dp[i + j * n];
      }
      grad_data[j] += 2.0 * sum_g;
    }
  }
  
  // -------------------------
  // 9. Hessian
  // -------------------------
  SEXP hessian = PROTECT(Rf_allocMatrix(REALSXP, n_pars, n_pars));
  double* H = REAL(hessian);
  std::fill(H, H + n_pars * n_pars, 0.0);
  
  if (!Rf_isNull(par_names)) {
    SEXP dimnames_H = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames_H, 0, par_names);
    SET_VECTOR_ELT(dimnames_H, 1, par_names);
    Rf_setAttrib(hessian, R_DimNamesSymbol, dimnames_H);
    UNPROTECT(1);
  }
  
  // ===========================================================
  // CASE A: NO second derivatives
  // ===========================================================
  if (Rf_isNull(deriv2)) {
    bool ALOQ_part1 = get_named_flag(opt_hessian, "ALOQ_part1", true);
    bool ALOQ_part2 = get_named_flag(opt_hessian, "ALOQ_part2", true);
    bool ALOQ_part3 = get_named_flag(opt_hessian, "ALOQ_part3", true);
    
    // Base Gaussian term: 2 * t(dwrdp) %*% dwrdp
    for (int j = 0; j < n_pars; ++j) {
      for (int k = 0; k < n_pars; ++k) {
        double sum_dw = 0.0;
        for (int i = 0; i < n; ++i) {
          sum_dw += dwrdp[i + j * n] * dwrdp[i + k * n];
        }
        H[j + k * n_pars] = 2.0 * sum_dw;
      }
    }
    
    if (has_dsdp) {
      std::vector<double> inv_s2(n), a(n), bvec(n);
      for (int i = 0; i < n; ++i) {
        double si = s[i];
        double si2 = si * si;
        inv_s2[i] = 1.0 / si2;
        a[i] = -wr_vec[i] * inv_s2[i];
        bvec[i] = 2.0 * wr_vec[i] * wr_vec[i] * inv_s2[i];
      }
      
      // ALOQ_part1
      if (ALOQ_part1) {
        for (int j = 0; j < n_pars; ++j) {
          for (int k = 0; k < n_pars; ++k) {
            double sum1 = 0.0, sum2 = 0.0;
            for (int i = 0; i < n; ++i) {
              double ai = a[i];
              double dxj = dxdp_data[i + j * n];
              double dxk = dxdp_data[i + k * n];
              double dsj = dsdp_data[i + j * n];
              double dsk = dsdp_data[i + k * n];
              sum1 += ai * dxj * dsk;
              sum2 += ai * dsj * dxk;
            }
            H[j + k * n_pars] += 2.0 * (sum1 + sum2);
          }
        }
      }
      
      // ALOQ_part2
      if (ALOQ_part2) {
        for (int j = 0; j < n_pars; ++j) {
          for (int k = 0; k < n_pars; ++k) {
            double sum_b = 0.0;
            for (int i = 0; i < n; ++i) {
              sum_b += bvec[i] * dsdp_data[i + j * n] * dsdp_data[i + k * n];
            }
            H[j + k * n_pars] += 2.0 * sum_b;
          }
        }
      }
      
      // ALOQ_part3
      if (ALOQ_part3) {
        for (int j = 0; j < n_pars; ++j) {
          for (int k = 0; k < n_pars; ++k) {
            double sum_log = 0.0;
            for (int i = 0; i < n; ++i) {
              sum_log += dlogsdp[i + j * n] * dlogsdp[i + k * n];
            }
            H[j + k * n_pars] -= 2.0 * sum_log;
          }
        }
      }
      
      // M4BEAL Hessian
      if (use_M4BEAL) {
        // Term 1
        for (int j = 0; j < n_pars; ++j) {
          for (int k = 0; k < n_pars; ++k) {
            double sum1 = 0.0;
            for (int i = 0; i < n; ++i) {
              double Gi = G_w0[i];
              double w0i = w0_vec[i];
              double coef = (-w0i * Gi - Gi * Gi);
              sum1 += coef * dw0dp[i + j * n] * dw0dp[i + k * n];
            }
            H[j + k * n_pars] += 2.0 * sum1;
          }
        }
        
        // Term 2
        for (int j = 0; j < n_pars; ++j) {
          for (int k = 0; k < n_pars; ++k) {
            double sum_dxds = 0.0;
            double sum_dsdx = 0.0;
            for (int i = 0; i < n; ++i) {
              double Gi = G_w0[i];
              double coef = -Gi * inv_s2[i];
              sum_dxds += coef * dxdp_data[i + j * n] * dsdp_data[i + k * n];
              sum_dsdx += coef * dsdp_data[i + j * n] * dxdp_data[i + k * n];
            }
            H[j + k * n_pars] += 2.0 * (sum_dxds + sum_dsdx);
          }
        }
        
        // Term 3
        if (ALOQ_part1) {
          for (int j = 0; j < n_pars; ++j) {
            for (int k = 0; k < n_pars; ++k) {
              double sum3 = 0.0;
              for (int i = 0; i < n; ++i) {
                double Gi = G_w0[i];
                double coef = 2.0 * Gi * w0_vec[i] * inv_s2[i];
                sum3 += coef * dsdp_data[i + j * n] * dsdp_data[i + k * n];
              }
              H[j + k * n_pars] += 2.0 * sum3;
            }
          }
        }
      }
    }
    
    // ===========================================================
    // CASE B: 2nd derivatives provided
    // ===========================================================
  } else {
    SEXP d2x_dim = Rf_getAttrib(deriv2, R_DimSymbol);
    if (Rf_length(d2x_dim) != 3)
      Rf_error("'deriv2' must be a 3D array");
    
    int* dims = INTEGER(d2x_dim);
    if (dims[0] != n || dims[1] != n_pars || dims[2] != n_pars)
      Rf_error("'deriv2' dimensions must be (n_data, n_pars, n_pars)");
    
    double* d2x_data = REAL(deriv2);
    
    bool has_d2s = false;
    double* d2s_data = nullptr;
    if (!Rf_isNull(deriv2_err)) {
      SEXP d2s_dim = Rf_getAttrib(deriv2_err, R_DimSymbol);
      int* dims_s = INTEGER(d2s_dim);
      if (Rf_length(d2s_dim) != 3 || dims_s[0] != n || dims_s[1] != n_pars || dims_s[2] != n_pars)
        Rf_error("'deriv2_err' must be a 3D array with dim = c(n_data, n_pars, n_pars)");
      
      d2s_data = REAL(deriv2_err);
      has_d2s = true;
    }
    
    // Lambda for 3D array access (column-major)
    auto d2x = [&](int i, int j, int k) -> double {
      return d2x_data[i + n * (j + n_pars * k)];
    };
    auto d2s = [&](int i, int j, int k) -> double {
      if (!has_d2s) return 0.0;
      return d2s_data[i + n * (j + n_pars * k)];
    };
    
    // Exact Hessian computation
    for (int j = 0; j < n_pars; ++j) {
      for (int k = 0; k < n_pars; ++k) {
        double Hjk = 0.0;
        
        for (int i = 0; i < n; ++i) {
          double si = s[i];
          double si2 = si * si;
          double inv_s = 1.0 / si;
          double inv_s2 = 1.0 / si2;
          double wr_i = wr_vec[i];
          double w0_i = w0_vec[i];
          
          int idx_j = i + j * n;
          int idx_k = i + k * n;
          
          double dx_j = dxdp_data[idx_j];
          double dx_k = dxdp_data[idx_k];
          double ds_j = dsdp_data[idx_j];
          double ds_k = dsdp_data[idx_k];
          
          double dwr_j = dwrdp[idx_j];
          double dwr_k = dwrdp[idx_k];
          
          double d2x_jk = d2x(i, j, k);
          double d2s_jk = d2s(i, j, k);
          
          // Second derivative of wr
          double d2wr_jk =
            (b * inv_s) * d2x_jk
          - (b * inv_s2) * (dx_j * ds_k + dx_k * ds_j)
            + (2.0 * wr_i * inv_s2) * ds_j * ds_k
            - (wr_i * inv_s) * d2s_jk;
            
            // Second derivative of log s
            double d2logs_jk = inv_s * d2s_jk - inv_s2 * ds_j * ds_k;
            
            // Gaussian ALOQ contribution
            double H_i = 2.0 * (dwr_j * dwr_k + wr_i * d2wr_jk) + 2.0 * d2logs_jk;
            
            // M4BEAL contribution
            if (use_M4BEAL) {
              double G_i = G_w0[i];
              double dw0_j = dw0dp[idx_j];
              double dw0_k = dw0dp[idx_k];
              
              double d2w0_jk =
                (b * inv_s) * d2x_jk
              - (b * inv_s2) * (dx_j * ds_k + dx_k * ds_j)
                + (2.0 * w0_i * inv_s2) * ds_j * ds_k
                - (w0_i * inv_s) * d2s_jk;
                
                double term_M4 =
                2.0 * ((-w0_i * G_i - G_i * G_i) * dw0_j * dw0_k
                         + G_i * d2w0_jk);
                
                H_i += term_M4;
            }
            
            Hjk += H_i;
        }
        
        H[j + k * n_pars] = Hjk;
      }
    }
  }
  
  // -------------------------
  // 10. Build output list
  // -------------------------
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  
  SET_STRING_ELT(names, 0, Rf_mkChar("value"));
  SET_STRING_ELT(names, 1, Rf_mkChar("gradient"));
  SET_STRING_ELT(names, 2, Rf_mkChar("hessian"));
  
  SEXP value = PROTECT(Rf_ScalarReal(obj));
  SET_VECTOR_ELT(out, 0, value);
  SET_VECTOR_ELT(out, 1, grad);
  SET_VECTOR_ELT(out, 2, hessian);
  
  Rf_setAttrib(out, R_NamesSymbol, names);
  
  // Attributes
  SEXP attr_chisq = PROTECT(Rf_ScalarReal(chisquare_ml));
  Rf_setAttrib(out, Rf_install("chisquare"), attr_chisq);
  
  SEXP attr_neg2ll = PROTECT(Rf_ScalarReal(neg2ll_ml));
  Rf_setAttrib(out, Rf_install("neg2ll"), attr_neg2ll);
  
  SEXP attr_bessel = PROTECT(Rf_ScalarLogical(use_bessel));
  Rf_setAttrib(out, Rf_install("besselcorrected"), attr_bessel);
  
  UNPROTECT(8);
  return out;
}

// =============================================================================
// C_nll_BLOQ
// =============================================================================

extern "C" SEXP C_nll_BLOQ(SEXP nout_bloq, SEXP derivs_bloq, SEXP derivs_err_bloq,
                          SEXP opt_BLOQ, SEXP opt_hessian) {
  
  // Extract required vectors
  int n;
  double* value = get_numeric_column(nout_bloq, "value", &n);
  
  int n_check;
  double* wr = get_numeric_column(nout_bloq, "weighted.residual", &n_check);
  if (n_check != n) Rf_error("Column lengths mismatch");
  
  double* w0 = get_numeric_column(nout_bloq, "weighted.0", &n_check);
  if (n_check != n) Rf_error("Column lengths mismatch");
  
  double* s = get_numeric_column(nout_bloq, "sigma", &n_check);
  if (n_check != n) Rf_error("Column lengths mismatch");
  
  // Objective value
  std::vector<double> objvals(n);
  const char* opt_BLOQ_str = CHAR(STRING_ELT(opt_BLOQ, 0));
  
  if (strcmp(opt_BLOQ_str, "M3") == 0) {
    for (int i = 0; i < n; ++i)
      objvals[i] = -2.0 * R::pnorm5(-wr[i], 0, 1, 1, 1);
  }
  else if (strcmp(opt_BLOQ_str, "M4NM") == 0 || strcmp(opt_BLOQ_str, "M4BEAL") == 0) {
    for (int i = 0; i < n; ++i) {
      double Phi_wr = R::pnorm5(wr[i], 0, 1, 1, 0);
      double Phi_w0 = R::pnorm5(w0[i], 0, 1, 1, 0);
      double r = Phi_wr / Phi_w0;
      objvals[i] = (r >= 1) ? R_PosInf : -2 * std::log(1 - r);
    }
    
    // Numeric fix
    for (int i = 0; i < n; ++i) {
      if (!R_finite(objvals[i])) {
        double d = w0[i] - wr[i];
        double ld = std::log(d);
        double intercept = (ld > 0 ? 1.8 : -1.9 * ld + 0.9);
        double lin = (ld > 0 ? 0.9 : 0.5);
        objvals[i] = intercept + lin * w0[i] + 0.95 * w0[i] * w0[i];
      }
    }
  }
  
  double obj = 0;
  for (int i = 0; i < n; ++i) obj += objvals[i];
  
  // No derivatives -> return early
  if (Rf_isNull(derivs_bloq)) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    
    SET_STRING_ELT(names, 0, Rf_mkChar("value"));
    SET_STRING_ELT(names, 1, Rf_mkChar("gradient"));
    SET_STRING_ELT(names, 2, Rf_mkChar("hessian"));
    
    SEXP value_out = PROTECT(Rf_ScalarReal(obj));
    SET_VECTOR_ELT(out, 0, value_out);
    SET_VECTOR_ELT(out, 1, R_NilValue);
    SET_VECTOR_ELT(out, 2, R_NilValue);
    
    Rf_setAttrib(out, R_NamesSymbol, names);
    
    SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(class_attr, 0, Rf_mkChar("objlist"));
    SET_STRING_ELT(class_attr, 1, Rf_mkChar("list"));
    Rf_setAttrib(out, R_ClassSymbol, class_attr);
    
    UNPROTECT(4);
    return out;
  }
  
  // Extract dxdp, dsdp
  int n_pars = Rf_ncols(derivs_bloq);
  double* dxdp_data = REAL(derivs_bloq);
  
  std::vector<double> dsdp_data(n * n_pars, 0.0);
  if (!Rf_isNull(derivs_err_bloq)) {
    if (Rf_nrows(derivs_err_bloq) != n || Rf_ncols(derivs_err_bloq) != n_pars)
      Rf_error("derivs_err_bloq wrong dims");
    
    double* derivs_err_data = REAL(derivs_err_bloq);
    std::copy(derivs_err_data, derivs_err_data + n * n_pars, dsdp_data.begin());
  }
  
  SEXP dimnames = Rf_getAttrib(derivs_bloq, R_DimNamesSymbol);
  SEXP par_names = R_NilValue;
  if (!Rf_isNull(dimnames) && Rf_length(dimnames) >= 2) {
    par_names = VECTOR_ELT(dimnames, 1);
  }
  
  // First-order derivatives
  std::vector<double> dwrdp(n * n_pars);
  std::vector<double> dw0dp(n * n_pars);
  std::vector<double> dlogsdp(n * n_pars);
  
  for (int i = 0; i < n; ++i) {
    double invs = 1.0 / s[i];
    for (int j = 0; j < n_pars; ++j) {
      int idx = i + j * n;
      double dx = dxdp_data[idx];
      double ds = dsdp_data[idx];
      dwrdp[idx] = invs * dx - (wr[i] * invs) * ds;
      dw0dp[idx] = invs * dx - (w0[i] * invs) * ds;
      dlogsdp[idx] = invs * ds;
    }
  }
  
  // Gradient
  SEXP grad = PROTECT(Rf_allocVector(REALSXP, n_pars));
  double* grad_data = REAL(grad);
  
  if (!Rf_isNull(par_names)) {
    Rf_setAttrib(grad, R_NamesSymbol, par_names);
  }
  
  if (strcmp(opt_BLOQ_str, "M3") == 0) {
    std::vector<double> Gm(n);
    for (int i = 0; i < n; ++i)
      Gm[i] = G_single(-wr[i], -wr[i]);
    
    for (int j = 0; j < n_pars; ++j) {
      double s_ = 0;
      for (int i = 0; i < n; ++i)
        s_ += Gm[i] * dwrdp[i + j * n];
      grad_data[j] = 2 * s_;
    }
  }
  else if (strcmp(opt_BLOQ_str, "M4NM") == 0 || strcmp(opt_BLOQ_str, "M4BEAL") == 0) {
    // Use vector versions
    std::vector<double> wr_vec(wr, wr + n);
    std::vector<double> w0_vec(w0, w0 + n);
    
    std::vector<double> G_wr_w0 = G_vec(wr_vec, w0_vec);
    std::vector<double> G_wr_wr = G_vec(wr_vec, wr_vec);
    std::vector<double> G_w0_w0 = G_vec(w0_vec, w0_vec);
    std::vector<double> G_w0_wr = G_vec(w0_vec, wr_vec);
    std::vector<double> G0 = G_vec(w0_vec);
    
    std::vector<double> c1(n), c2(n), c3(n);
    for (int i = 0; i < n; ++i) {
      c1[i] = 2.0 / (1.0 / G_wr_w0[i] - 1.0 / G_wr_wr[i]);
      c2[i] = 2.0 / (1.0 / G_w0_w0[i] - 1.0 / G_w0_wr[i]);
      c3[i] = 2.0 * G0[i];
    }
    
    for (int j = 0; j < n_pars; ++j) {
      double S1 = 0, S2 = 0, S3 = 0;
      for (int i = 0; i < n; ++i) {
        S1 += c1[i] * dwrdp[i + j * n];
        S2 += c2[i] * dw0dp[i + j * n];
        S3 += c3[i] * dw0dp[i + j * n];
      }
      grad_data[j] = S1 - S2 + S3;
    }
  }
  
  // Hessian
  SEXP hessian = PROTECT(Rf_allocMatrix(REALSXP, n_pars, n_pars));
  double* H = REAL(hessian);
  std::fill(H, H + n_pars * n_pars, 0.0);
  
  if (!Rf_isNull(par_names)) {
    SEXP dimnames_H = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames_H, 0, par_names);
    SET_VECTOR_ELT(dimnames_H, 1, par_names);
    Rf_setAttrib(hessian, R_DimNamesSymbol, dimnames_H);
    UNPROTECT(1);
  }
  
  bool P1 = get_named_flag(opt_hessian, "BLOQ_part1", true);
  bool P2 = get_named_flag(opt_hessian, "BLOQ_part2", true);
  bool P3 = get_named_flag(opt_hessian, "BLOQ_part3", true);
  
  if (strcmp(opt_BLOQ_str, "M3") == 0) {
    std::vector<double> Gm(n);
    for (int i = 0; i < n; ++i)
      Gm[i] = G_single(-wr[i], -wr[i]);
    
    // part1
    if (P1) {
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          double sum_ = 0;
          for (int i = 0; i < n; ++i) {
            double coef = (-wr[i] * Gm[i] + Gm[i] * Gm[i]);
            sum_ += coef * dwrdp[i + j * n] * dwrdp[i + k * n];
          }
          H[j + k * n_pars] += 2 * sum_;
        }
      }
    }
    
    // part2
    if (P2) {
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          double s1 = 0, s2 = 0;
          for (int i = 0; i < n; ++i) {
            double invs2 = 1.0 / (s[i] * s[i]);
            double coef = Gm[i] * invs2;
            s1 += coef * dxdp_data[i + j * n] * dsdp_data[i + k * n];
            s2 += coef * dsdp_data[i + j * n] * dxdp_data[i + k * n];
          }
          H[j + k * n_pars] -= 2 * (s1 + s2);
        }
      }
    }
    
    // part3
    if (P3) {
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          double sum3 = 0;
          for (int i = 0; i < n; ++i) {
            double invs2 = 1.0 / (s[i] * s[i]);
            double coef = Gm[i] * (2 * (-wr[i])) * invs2;
            sum3 += coef * dsdp_data[i + j * n] * dsdp_data[i + k * n];
          }
          H[j + k * n_pars] -= 2 * sum3;
        }
      }
    }
  }
  else if (strcmp(opt_BLOQ_str, "M4NM") == 0 || strcmp(opt_BLOQ_str, "M4BEAL") == 0) {
    // Convert to vectors for helper functions
    std::vector<double> wr_vec(wr, wr + n);
    std::vector<double> w0_vec(w0, w0 + n);
    std::vector<double> s_vec(s, s + n);
    
    std::vector<double> stable_wr = stable_fun(wr_vec, w0_vec, wr_vec);
    std::vector<double> stable_w0 = stable_fun(w0_vec, w0_vec, wr_vec);
    std::vector<double> G0 = G_vec(w0_vec);
    
    std::vector<double> A1(n), A2(n), A3(n), A4(n), A5(n), A6(n);
    for (int i = 0; i < n; ++i) {
      A1[i] = -wr[i] * stable_wr[i];
      A2[i] = stable_wr[i];
      A3[i] = -w0[i] * stable_w0[i];
      A4[i] = stable_w0[i];
      A5[i] = -w0[i] * G0[i] - G0[i] * G0[i];
      A6[i] = G0[i];
    }
    
    // Temporary matrices for Hessian parts
    std::vector<double> tmp11(n_pars * n_pars);
    std::vector<double> tmp12(n_pars * n_pars);
    std::vector<double> tmp21(n_pars * n_pars);
    std::vector<double> tmp22(n_pars * n_pars);
    std::vector<double> tmp31(n_pars * n_pars);
    std::vector<double> tmp32(n_pars * n_pars);
    
    // part1 contributions
    d_dp_sq_cpp(A1.data(), wr, s, dxdp_data, dsdp_data.data(), 
                tmp11.data(), n, n_pars);
    d2_dp2_cpp(A2.data(), wr, s, dxdp_data, dsdp_data.data(), 
               tmp12.data(), n, n_pars);
    d_dp_sq_cpp(A3.data(), w0, s, dxdp_data, dsdp_data.data(), 
                tmp21.data(), n, n_pars);
    d2_dp2_cpp(A4.data(), w0, s, dxdp_data, dsdp_data.data(), 
               tmp22.data(), n, n_pars);
    
    if (P1) {
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          int idx = j + k * n_pars;
          H[idx] += 2.0 * (tmp11[idx] + tmp12[idx] - tmp21[idx] - tmp22[idx]);
        }
      }
    }
    
    // part2
    if (P2) {
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          double sum = 0.0;
          for (int i = 0; i < n; ++i) {
            double aij = stable_wr[i] * dwrdp[i + j * n] 
            - stable_w0[i] * dw0dp[i + j * n];
            double aik = stable_wr[i] * dwrdp[i + k * n] 
            - stable_w0[i] * dw0dp[i + k * n];
            sum += aij * aik;
          }
          H[j + k * n_pars] += -2.0 * sum;
        }
      }
    }
    
    // part3
    if (P3) {
      d_dp_sq_cpp(A5.data(), w0, s, dxdp_data, dsdp_data.data(), 
                  tmp31.data(), n, n_pars);
      d2_dp2_cpp(A6.data(), w0, s, dxdp_data, dsdp_data.data(), 
                 tmp32.data(), n, n_pars);
      
      for (int j = 0; j < n_pars; ++j) {
        for (int k = 0; k < n_pars; ++k) {
          int idx = j + k * n_pars;
          H[idx] += 2.0 * (tmp31[idx] + tmp32[idx]);
        }
      }
    }
  }
  
  // Build output
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  
  SET_STRING_ELT(names, 0, Rf_mkChar("value"));
  SET_STRING_ELT(names, 1, Rf_mkChar("gradient"));
  SET_STRING_ELT(names, 2, Rf_mkChar("hessian"));
  
  SEXP value_out = PROTECT(Rf_ScalarReal(obj));
  SET_VECTOR_ELT(out, 0, value_out);
  SET_VECTOR_ELT(out, 1, grad);
  SET_VECTOR_ELT(out, 2, hessian);
  
  Rf_setAttrib(out, R_NamesSymbol, names);
  
  SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(class_attr, 0, Rf_mkChar("objlist"));
  SET_STRING_ELT(class_attr, 1, Rf_mkChar("list"));
  Rf_setAttrib(out, R_ClassSymbol, class_attr);
  
  UNPROTECT(6);
  return out;
}
