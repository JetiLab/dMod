#ifndef NLL_H
#define NLL_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif
  
  SEXP C_nll_ALOQ(SEXP nout, SEXP derivs, SEXP derivs_err, 
                  SEXP opt_BLOQ, SEXP opt_hessian, SEXP bessel_correction,
                  SEXP deriv2, SEXP deriv2_err);
  
  SEXP C_nll_BLOQ(SEXP nout_bloq, SEXP derivs_bloq, SEXP derivs_err_bloq,
                  SEXP opt_BLOQ, SEXP opt_hessian);
  
#ifdef __cplusplus
}
#endif

#endif // NLL_H