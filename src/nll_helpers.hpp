#ifndef NLL_HELPERS_HPP
#define NLL_HELPERS_HPP

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// ============================================================================
// 1) get_flag  (named logical flag from R)
// ============================================================================
inline bool get_flag(const std::vector<int>& lv,
                     const std::vector<std::string>& names,
                     const std::string& name,
                     bool default_value = true)
{
  if (lv.empty()) return default_value;
  if (names.size() != lv.size()) return default_value;
  
  for (size_t i = 0; i < lv.size(); ++i) {
    if (names[i] == name) {
      if (lv[i] == NA_LOGICAL) return false;
      return lv[i] != 0;
    }
  }
  return default_value;
}

// ============================================================================
// 2) LOG-SAFE NORMAL HELPERS
// ============================================================================
inline double log_phi(double x) {
  return dnorm(x, 0.0, 1.0, 1);  // log=TRUE
}

inline double log_Phi(double x) {
  return pnorm(x, 0.0, 1.0, 1, 1);  // lower.tail=TRUE, log.p=TRUE
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

// Vectorized version using std::vector
inline std::vector<double> G_vec(const std::vector<double>& w1,
                                 const std::vector<double>& w2)
{
  size_t n = w1.size();
  std::vector<double> out(n);
  for (size_t i = 0; i < n; ++i)
    out[i] = G_single(w1[i], w2[i]);
  return out;
}

inline std::vector<double> G_vec(const std::vector<double>& w)
{
  return G_vec(w, w);
}

// Pointer-based versions (for direct use with R data)
inline void G_vec_ptr(const double* w1, const double* w2, double* out, int n)
{
  for (int i = 0; i < n; ++i)
    out[i] = G_single(w1[i], w2[i]);
}

inline void G_vec_ptr(const double* w, double* out, int n)
{
  G_vec_ptr(w, w, out, n);
}

// ============================================================================
// 4) stable_fun  (used in M4-BEAL / M4-NM Hessian)
// ============================================================================
inline std::vector<double> stable_fun(const std::vector<double>& wn,
                                      const std::vector<double>& w0,
                                      const std::vector<double>& wr)
{
  size_t n = wn.size();
  std::vector<double> out(n);
  
  for (size_t i = 0; i < n; ++i) {
    
    double Phi0 = pnorm(w0[i], 0, 1, 1, 0);  // lower.tail=TRUE, log.p=FALSE
    double Phir = pnorm(wr[i], 0, 1, 1, 0);
    double denom = Phi0 - Phir;
    
    // Case wn == w0
    if (wn[i] == w0[i]) {
      if (denom <= 0)
        out[i] = 0.0;
      else
        out[i] = dnorm(w0[i], 0, 1, 0) / denom;  // log=FALSE
    }
    // Case wn == wr
    else {
      if (denom <= 0)
        out[i] = 1.0/(w0[i] - wr[i]) + wr[i];
      else
        out[i] = dnorm(wr[i], 0, 1, 0) / denom;
    }
  }
  
  return out;
}

// Pointer-based version
inline void stable_fun_ptr(const double* wn, const double* w0, const double* wr,
                           double* out, int n)
{
  for (int i = 0; i < n; ++i) {
    
    double Phi0 = pnorm(w0[i], 0, 1, 1, 0);
    double Phir = pnorm(wr[i], 0, 1, 1, 0);
    double denom = Phi0 - Phir;
    
    // Case wn == w0
    if (wn[i] == w0[i]) {
      if (denom <= 0)
        out[i] = 0.0;
      else
        out[i] = dnorm(w0[i], 0, 1, 0) / denom;
    }
    // Case wn == wr
    else {
      if (denom <= 0)
        out[i] = 1.0/(w0[i] - wr[i]) + wr[i];
      else
        out[i] = dnorm(wr[i], 0, 1, 0) / denom;
    }
  }
}

// ============================================================================
// 5) d_dp_sq_cpp : A * (dw/dp)(dw/dp)^T
//    Matrix stored in column-major order
// ============================================================================
inline void d_dp_sq_cpp(const double* A,      // [n]
                        const double* w,      // [n]
                        const double* s,      // [n]
                        const double* dxdp,   // [n x p] column-major
                        const double* dsdp,   // [n x p] column-major
                        double* out,          // [p x p] column-major
                        int n,
                        int p)
{
  std::fill(out, out + p * p, 0.0);
  
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      
      double sum = 0.0;
      
      for (int i = 0; i < n; ++i) {
        
        double invs = 1.0 / s[i];
        
        double dw_j = invs * dxdp[i + j * n] - (w[i] * invs) * dsdp[i + j * n];
        double dw_k = invs * dxdp[i + k * n] - (w[i] * invs) * dsdp[i + k * n];
        
        sum += A[i] * dw_j * dw_k;
      }
      
      out[j + k * p] = sum;
    }
  }
}

// ============================================================================
// 6) d2_dp2_cpp: A * second derivatives part
//    Matrix stored in column-major order
// ============================================================================
inline void d2_dp2_cpp(const double* A,      // [n]
                       const double* w,      // [n]
                       const double* s,      // [n]
                       const double* dxdp,   // [n x p] column-major
                       const double* dsdp,   // [n x p] column-major
                       double* out,          // [p x p] column-major
                       int n,
                       int p)
{
  std::fill(out, out + p * p, 0.0);
  
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      
      double sum = 0.0;
      
      for (int i = 0; i < n; ++i) {
        
        double invs2 = 1.0 / (s[i] * s[i]);
        
        double term =
          A[i] * (
              - invs2 * dxdp[i + j * n] * dsdp[i + k * n]
        - invs2 * dsdp[i + j * n] * dxdp[i + k * n]
        + (2.0 * w[i]) * invs2 * dsdp[i + j * n] * dsdp[i + k * n]
          );
        
        sum += term;
      }
      
      out[j + k * p] = sum;
    }
  }
}

// ============================================================================
// 7) Matrix addition utility (for combining Hessian parts)
// ============================================================================
inline void matrix_add(double* dest, const double* src, int n_elem)
{
  for (int i = 0; i < n_elem; ++i)
    dest[i] += src[i];
}

inline void matrix_add_scaled(double* dest, const double* src, double scale, int n_elem)
{
  for (int i = 0; i < n_elem; ++i)
    dest[i] += scale * src[i];
}

#endif  // NLL_HELPERS_HPP