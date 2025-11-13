#ifndef NLL_HELPERS_HPP
#define NLL_HELPERS_HPP

#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// ============================================================================
// 1) get_flag  (named logical flag from R)
// ============================================================================
inline bool get_flag(const LogicalVector& lv,
                     const std::string& name,
                     bool default_value = true)
{
  if (lv.size() == 0) return default_value;
  
  CharacterVector nms = lv.names();
  if (nms.size() != lv.size()) return default_value;
  
  for (int i = 0; i < lv.size(); ++i) {
    if (TYPEOF(nms[i]) == STRSXP &&
        as<std::string>(nms[i]) == name)
    {
      if (LogicalVector::is_na(lv[i])) return false;
      return lv[i];
    }
  }
  return default_value;
}

// ============================================================================
// 2) LOG-SAFE NORMAL HELPERS
// ============================================================================
inline double log_phi(double x) {
  return R::dnorm4(x, 0.0, 1.0, 1);
}

inline double log_Phi(double x) {
  return R::pnorm5(x, 0.0, 1.0, 1, 1);
}

// ============================================================================
// 3) G = φ/Φ  (scalar + vectorized + 2-argument)
// ============================================================================

// Single argument φ(w)/Φ(w)
inline double G_single(double w)
{
  if (w < -35.0) return 0.0;   // avoid overflow
  return std::exp(log_phi(w) - log_Phi(w));
}

// Two-argument φ(w1)/Φ(w2)
inline double G_single(double w1, double w2)
{
  if (w2 < -35.0) return 0.0;
  return std::exp(log_phi(w1) - log_Phi(w2));
}

// Vectorized
inline NumericVector G_vec(const NumericVector& w1,
                           const NumericVector& w2)
{
  int n = w1.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i)
    out[i] = G_single(w1[i], w2[i]);
  return out;
}

inline NumericVector G_vec(const NumericVector& w)
{
  return G_vec(w, w);
}

// ============================================================================
// 4) stable_fun  (used in M4-BEAL / M4-NM Hessian)
// ============================================================================
inline NumericVector stable_fun(const NumericVector& wn,
                                const NumericVector& w0,
                                const NumericVector& wr)
{
  int n = wn.size();
  NumericVector out(n);
  
  for (int i = 0; i < n; ++i) {
    
    double Phi0 = R::pnorm5(w0[i], 0, 1, 1, 0);
    double Phir = R::pnorm5(wr[i], 0, 1, 1, 0);
    double denom = Phi0 - Phir;
    
    // Case wn == w0
    if (wn[i] == w0[i]) {
      if (denom <= 0)
        out[i] = 0.0;
      else
        out[i] = R::dnorm4(w0[i], 0, 1, 0) / denom;
    }
    // Case wn == wr
    else {
      if (denom <= 0)
        out[i] = 1.0/(w0[i] - wr[i]) + wr[i];
      else
        out[i] = R::dnorm4(wr[i], 0, 1, 0) / denom;
    }
  }
  
  return out;
}

// ============================================================================
// 5) d_dp_sq_cpp : A * (dw/dp)(dw/dp)^T
// ============================================================================
inline NumericMatrix d_dp_sq_cpp(const NumericVector& A,
                                 const NumericVector& w,
                                 const NumericVector& s,
                                 const NumericMatrix& dxdp,
                                 const NumericMatrix& dsdp)
{
  int n = w.size();
  int p = dxdp.ncol();
  
  NumericMatrix out(p, p);
  
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      
      double sum = 0.0;
      
      for (int i = 0; i < n; ++i) {
        
        double invs = 1.0 / s[i];
        
        double dw_j = invs*dxdp(i,j) - (w[i]*invs)*dsdp(i,j);
        double dw_k = invs*dxdp(i,k) - (w[i]*invs)*dsdp(i,k);
        
        sum += A[i] * dw_j * dw_k;
      }
      
      out(j,k) = sum;
    }
  }
  
  return out;
}

// ============================================================================
// 6) d2_dp2_cpp: A * second derivatives part
// ============================================================================
inline NumericMatrix d2_dp2_cpp(const NumericVector& A,
                                const NumericVector& w,
                                const NumericVector& s,
                                const NumericMatrix& dxdp,
                                const NumericMatrix& dsdp)
{
  int n = w.size();
  int p = dxdp.ncol();
  
  NumericMatrix out(p, p);
  
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      
      double sum = 0.0;
      
      for (int i = 0; i < n; ++i) {
        
        double invs2 = 1.0 / (s[i] * s[i]);
        
        double term =
          A[i] * (
              - invs2 * dxdp(i,j) * dsdp(i,k)
          - invs2 * dsdp(i,j) * dxdp(i,k)
          + (2.0 * w[i]) * invs2 * dsdp(i,j) * dsdp(i,k)
          );
        
        sum += term;
      }
      
      out(j,k) = sum;
    }
  }
  
  return out;
}

#endif  // NLL_HELPERS_HPP
