// trust-region optimiser with L1 penalty
//
// Minimises  objfun(theta) + sum_i lambda_i * |theta_i - mu_i|
// (or one-sided  lambda_i * max(0, mu_i - theta_i)) over a subset of
// parameters named in mu, on top of the Moré-Sorensen subproblem in
// trust_subproblem.h. Two L1-specific interventions sit in the inner loop:
//
//  1. L1 active set: a penalised coordinate sitting exactly on its kink
//     theta_i = mu_i is dropped from the reduced subproblem whenever
//     |grad_obj_i| <= lambda_i (two-sided) resp. -grad_obj_i <= lambda_i
//     (one-sided). Below that threshold the L1 force dominates and the
//     coordinate stays pinned.
//
//  2. Kink clamping: after the trust step is added to theta, any penalised
//     coordinate that crossed its kink is snapped back to mu_i, so the next
//     iteration can re-examine the active set. The assignment is verbatim --
//     downstream sparsity tests read an exact equality.
//
// Box bounds are handled by `boundary`, exactly as in trust_kernel.cpp:
// "reflective" applies the Coleman-Li scaling to the coordinates that survive
// the L1 active set, "clip" is the frozen historical scheme. The kink active
// set is orthogonal to that choice and is used by both -- L1 sparsity needs
// coordinates to land exactly on mu, which an interior method cannot deliver.

#include <Rcpp.h>
#include "trust_subproblem.h"
#include "trust_driver.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace Rcpp;
using dmod::trust_internal::affine_scaling;
using dmod::trust_internal::eigen_sym_local;
using dmod::trust_internal::model_value;
using dmod::trust_internal::stepback;
using dmod::trust_internal::trust_sub;
using dmod::trust_driver::Blather;
using dmod::trust_driver::Reporter;
using dmod::trust_driver::eval_objfun;
using dmod::trust_driver::fill_bound;
using dmod::trust_driver::fill_parscale;
using dmod::trust_driver::kInf;
using dmod::trust_driver::kStallLimit;
using dmod::trust_driver::push_interior;
using dmod::trust_driver::subproblem_label;

namespace {

// Per-parameter L1 metadata, resolved from the named (mu, lambda) pair.
struct L1Spec {
  std::vector<unsigned char> has;
  std::vector<double>        mu;
  std::vector<double>        lambda;
  bool                       one_sided = false;

  void build(const NumericVector& mu_in, const NumericVector& lambda_in,
             const CharacterVector& parnames, int K, bool one_sided_) {
    has.assign(K, 0);
    mu.assign(K, 0.0);
    lambda.assign(K, 0.0);
    one_sided = one_sided_;
    if (mu_in.size() == 0) return;
    if (!mu_in.hasAttribute("names"))
      stop("trustL1: mu must be a named numeric vector");
    if (lambda_in.size() != mu_in.size())
      stop("trustL1: lambda must have the same length as mu");
    CharacterVector mnames = mu_in.names();
    for (int j = 0; j < mu_in.size(); ++j) {
      std::string nm = as<std::string>(mnames[j]);
      for (int i = 0; i < K; ++i)
        if (as<std::string>(parnames[i]) == nm) {
          has[i] = 1; mu[i] = mu_in[j]; lambda[i] = lambda_in[j];
          break;
        }
    }
  }

  double value(const std::vector<double>& th) const {
    double s = 0.0;
    for (std::size_t i = 0; i < th.size(); ++i) {
      if (!has[i]) continue;
      const double d = th[i] - mu[i];
      if (one_sided) { if (d < 0.0) s += lambda[i] * (-d); }
      else           { s += lambda[i] * std::fabs(d); }
    }
    return s;
  }
  // Subgradient of the penalty. At the kink the coordinate is pinned by the
  // active set and dropped, so the zero returned there is never used.
  double grad(int i, double th_i) const {
    if (!has[i]) return 0.0;
    const double d = th_i - mu[i];
    if (one_sided) return (d < 0.0) ? -lambda[i] : 0.0;
    if (d > 0.0) return  lambda[i];
    if (d < 0.0) return -lambda[i];
    return 0.0;
  }
  bool pinned(int i, double th_i, double grad_obj_i) const {
    if (!has[i] || th_i != mu[i]) return false;
    if (one_sided) return (-grad_obj_i) <= lambda[i];
    return std::fabs(grad_obj_i) <= lambda[i];
  }
  // Snap any coordinate that crossed its kink back onto mu, verbatim.
  void clamp(const std::vector<double>& th, std::vector<double>& th_try) const {
    for (std::size_t i = 0; i < th.size(); ++i) {
      if (!has[i]) continue;
      if ((th[i] - mu[i]) * (th_try[i] - mu[i]) < 0.0) th_try[i] = mu[i];
      if (one_sided && th_try[i] < mu[i])              th_try[i] = mu[i];
    }
  }
};

// -------------------------------------------------------------------------
// Coleman-Li interior trust-region-reflective on the L1-active coordinates
// -------------------------------------------------------------------------
List trustL1_reflective(Function objfun, NumericVector parinit,
                        const L1Spec& l1,
                        double rinit, double rmax,
                        Nullable<NumericVector> parscale,
                        int iterlim,
                        double ftol, double mtol,
                        double gtol, double xtol,
                        double rmin, double thetamax,
                        bool minimize, bool blather_on,
                        Nullable<NumericVector> parupper,
                        Nullable<NumericVector> parlower,
                        bool printIter, Nullable<CharacterVector> traceFile) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  for (int i = 0; i < K; ++i) {
    if (!(pl[i] < pu[i]))
      stop("trustL1: parlower must lie strictly below parupper");
    // The kink is a target the iterate must be able to sit on exactly, so it
    // has to be reachable from inside the box.
    if (l1.has[i] && !(pl[i] < l1.mu[i] && l1.mu[i] < pu[i]))
      stop("trustL1: mu must lie strictly inside [parlower, parupper]");
  }
  fill_parscale(parscale, ps, K, "trustL1");

  std::vector<double> lbz(K), ubz(K), z(K);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    lbz[i] = ps[i] * pl[i];
    ubz[i] = ps[i] * pu[i];
    double zi = ps[i] * parinit[i];
    if (zi > ubz[i]) { any_above = true; zi = ubz[i]; }
    if (zi < lbz[i]) { any_below = true; zi = lbz[i]; }
    z[i] = zi;
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");
  push_interior(K, z, lbz, ubz);

  std::vector<double> theta(K);
  for (int i = 0; i < K; ++i) theta[i] = z[i] / ps[i];

  Reporter report;
  report.init(printIter, iterlim, traceFile, K, parnames);

  NumericVector x_named(theta.begin(), theta.end());
  x_named.names() = parnames;
  List out_init;
  if (!eval_objfun(objfun, x_named, out_init))
    stop("parinit not feasible: objfun failed");
  double val_obj = as<double>(out_init["value"]);
  if (!std::isfinite(val_obj)) stop("parinit not feasible: value is not finite");
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);

  std::vector<double> grad_obj(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  double val = val_obj + l1.value(theta);
  int neval = 1;
  report(neval, val, x_named, /*head=*/true);

  Blather trace;
  std::vector<int>    active;
  std::vector<double> zr, lbr, ubr, gr, Hr, absv, jv, sqrtv, ghat, Bhat;
  std::vector<double> eigvals, eigvecs, shat, shat_step, s_step, shat_real;
  std::vector<double> z_try(K), theta_try(K);
  std::vector<unsigned char> at_bound(K, 0);

  double r = rinit;
  double f_used = minimize ? val : -val;
  double opt_measure = kInf;
  bool   accept = true, converged = false, bail = false;
  int    n_iter = 0, n_fail = 0, n_stall = 0;
  std::string stop_reason = "iterlim";

  for (int iter = 1; iter <= iterlim; ++iter) {
    R_CheckUserInterrupt();

    const double sgn = minimize ? 1.0 : -1.0;
    int Kred = 0;

    if (accept) {
      f_used = minimize ? val : -val;

      // Reduced space: drop only the coordinates pinned at their L1 kink. Box
      // bounds are handled by the scaling, not by dropping.
      active.clear();
      for (int i = 0; i < K; ++i)
        if (!l1.pinned(i, theta[i], grad_obj[i])) active.push_back(i);
      Kred = static_cast<int>(active.size());

      zr.assign(Kred, 0.0); lbr.assign(Kred, 0.0); ubr.assign(Kred, 0.0);
      gr.assign(Kred, 0.0); Hr.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        const int i = active[ii];
        zr[ii]  = z[i];
        lbr[ii] = lbz[i];
        ubr[ii] = ubz[i];
        gr[ii]  = sgn * (grad_obj[i] + l1.grad(i, theta[i])) / ps[i];
        for (int jj = 0; jj < Kred; ++jj) {
          const int j = active[jj];
          Hr[ii + (std::size_t) jj * Kred] =
              sgn * H_full[i + (std::size_t) j * K] / (ps[i] * ps[j]);
        }
      }

      absv.assign(Kred, 0.0); jv.assign(Kred, 0.0); sqrtv.assign(Kred, 0.0);
      affine_scaling(Kred, zr.data(), gr.data(), lbr.data(), ubr.data(),
                     absv.data(), jv.data());
      for (int ii = 0; ii < Kred; ++ii) sqrtv[ii] = std::sqrt(absv[ii]);

      // Pinned coordinates are stationary by the active-set test, so the
      // optimality measure only has to cover the reduced space.
      opt_measure = 0.0;
      for (int ii = 0; ii < Kred; ++ii)
        opt_measure = std::max(opt_measure, std::fabs(absv[ii] * gr[ii]));
      // Set before the convergence break, which is exactly when it matters.
      const double btol = std::max(gtol, 1e-10);
      std::fill(at_bound.begin(), at_bound.end(), 0);
      for (int ii = 0; ii < Kred; ++ii)
        at_bound[active[ii]] = (jv[ii] > 0.0 &&
                                std::fabs(absv[ii] * gr[ii]) <= btol &&
                                std::fabs(gr[ii]) > btol) ? 1 : 0;

      if (opt_measure <= gtol) { converged = true; stop_reason = "gradient"; break; }

      ghat.assign(Kred, 0.0);
      Bhat.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) ghat[ii] = sqrtv[ii] * gr[ii];
      for (int jj = 0; jj < Kred; ++jj)
        for (int ii = 0; ii < Kred; ++ii)
          Bhat[ii + (std::size_t) jj * Kred] =
              sqrtv[ii] * sqrtv[jj] * Hr[ii + (std::size_t) jj * Kred];
      for (int ii = 0; ii < Kred; ++ii)
        Bhat[ii + (std::size_t) ii * Kred] += std::fabs(gr[ii]) * jv[ii];

      eigvals.assign(Kred, 0.0);
      eigvecs.assign((std::size_t) Kred * Kred, 0.0);
      if (Kred > 0)
        eigen_sym_local(Bhat.data(), Kred, eigvals.data(), eigvecs.data());
    }

    Kred = static_cast<int>(active.size());

    if (blather_on) {
      trace.argpath.insert(trace.argpath.end(), theta.begin(), theta.end());
      trace.r.push_back(r);
      trace.valpath.push_back(val);
    }

    bool is_newton = false, is_hard = false, is_easy = false;
    const char* sb_label = "full";
    double m_value = 0.0;
    shat.assign(Kred, 0.0);
    shat_step.assign(Kred, 0.0);
    s_step.assign(Kred, 0.0);
    if (Kred > 0) {
      double pred_unused = 0.0;
      trust_sub(Kred, ghat.data(), eigvals.data(), eigvecs.data(), r,
                shat.data(), &pred_unused, &is_newton, &is_hard, &is_easy);
      const double theta_frac =
          std::min(std::max(thetamax, 1.0 - opt_measure), 1.0 - 1e-12);
      stepback(Kred, zr.data(), lbr.data(), ubr.data(), sqrtv.data(),
               ghat.data(), Bhat.data(), shat.data(), r, theta_frac,
               shat_step.data(), s_step.data(), &m_value, &sb_label);
    }

    z_try = z;
    for (int ii = 0; ii < Kred; ++ii) z_try[active[ii]] += s_step[ii];
    // See trust_kernel.cpp. mu is validated strictly inside the box, so this
    // never disturbs a pinned kink.
    push_interior(K, z_try, lbz, ubz);
    for (int i = 0; i < K; ++i) theta_try[i] = z_try[i] / ps[i];
    l1.clamp(theta, theta_try);
    for (int i = 0; i < K; ++i) z_try[i] = ps[i] * theta_try[i];

    // Rescore the model at the step actually taken -- the kink clamp shortens
    // individual coordinates after the stepback has chosen a candidate.
    shat_real.assign(Kred, 0.0);
    bool rescore = true;
    for (int ii = 0; ii < Kred; ++ii) {
      if (!(sqrtv[ii] > 0.0)) { rescore = false; break; }
      shat_real[ii] = (z_try[active[ii]] - zr[ii]) / sqrtv[ii];
    }
    if (rescore && Kred > 0)
      m_value = model_value(Kred, ghat.data(), Bhat.data(), shat_real.data());
    const std::vector<double>& shat_used = rescore ? shat_real : shat_step;

    double stepnorm = 0.0;
    for (int ii = 0; ii < Kred; ++ii) stepnorm += shat_used[ii] * shat_used[ii];
    stepnorm = std::sqrt(stepnorm);

    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try);
    double val_obj_try = kInf, val_try = kInf;
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_obj_try = as<double>(out_try["value"]);
      grad_try    = as<NumericVector>(out_try["gradient"]);
      Htry_mat    = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_obj_try)) eval_ok = false;
      else val_try = val_obj_try + l1.value(theta_try);
    }
    neval++;
    report(neval, val_try, x_try, /*head=*/false);

    const double pred_pos  = -m_value;
    const double ftry_used = minimize ? val_try : -val_try;
    const double dval      = std::fabs(ftry_used - f_used);
    const double rho = (eval_ok && pred_pos > 0.0)
                         ? (f_used - ftry_used) / pred_pos : -kInf;

    if (!eval_ok) {
      n_fail++;
      accept = false;
      r *= 0.25;
      if (n_fail >= 3) { bail = true; stop_reason = "objfun"; }
    } else {
      n_fail = 0;
      if (rho < 0.25) {
        accept = false;
        r = std::min(0.25 * r, 0.25 * stepnorm);
      } else {
        accept = true;
        if (rho > 0.75 && stepnorm >= 0.9 * r) r = std::min(2.0 * r, rmax);
      }
    }

    // Count rejected steps that also left the objective flat.
    if (accept || !eval_ok || dval >= ftol) n_stall = 0;
    else                                                              n_stall++;

    if (accept && eval_ok) {
      theta = theta_try;
      z     = z_try;
      val   = val_try;
      grad_obj.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
    }

    if (blather_on) {
      trace.argtry.insert(trace.argtry.end(), theta_try.begin(), theta_try.end());
      trace.valtry.push_back(val_try);
      trace.accept.push_back(accept ? 1 : 0);
      trace.preddiff.push_back(m_value);
      trace.stepnorm.push_back(stepnorm);
      trace.rho.push_back(rho);
      trace.steptype.push_back(subproblem_label(is_newton, is_hard, is_easy));
      trace.stepback.push_back(sb_label);
    }
    n_iter = iter;

    if (bail) break;
    if (accept && eval_ok) {
      if (dval < ftol) { converged = true; stop_reason = "fvalue"; break; }
      if (std::fabs(m_value) < mtol)             { converged = true; stop_reason = "preddiff";  break; }
      if (xtol > 0.0 && stepnorm < xtol)         { converged = true; stop_reason = "step";   break; }
    }
    if (r < rmin) { stop_reason = "radius"; break; }
    if (n_stall >= kStallLimit) { converged = true; stop_reason = "stagnation"; break; }
  }

  if (stop_reason == "objfun")
    Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
  else if (stop_reason == "radius")
    Rf_warning("Trust radius fell below rmin. Fit is not converged.");
  else if (!converged && n_iter >= iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  NumericVector arg_out(theta.begin(), theta.end());
  arg_out.names() = parnames;
  // The combined gradient, so the result is self-consistent with `value`.
  NumericVector grad_out(K);
  for (int i = 0; i < K; ++i) grad_out[i] = grad_obj[i] + l1.grad(i, theta[i]);
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) Hess_out(i, j) = H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (at_bound[i] != 0);
  at_bound_out.names() = parnames;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = n_iter,
      Named("converged")  = converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = stop_reason);

  if (blather_on) trace.attach(result, n_iter, K, parnames, minimize);
  return result;
}

// -------------------------------------------------------------------------
// Legacy: active-set reduction plus componentwise clipping
// -------------------------------------------------------------------------
List trustL1_clip(Function objfun, NumericVector parinit,
                  const L1Spec& l1,
                  double rinit, double rmax,
                  Nullable<NumericVector> parscale,
                  int iterlim, double fterm, double mterm,
                  bool minimize, bool blather_on,
                  Nullable<NumericVector> parupper,
                  Nullable<NumericVector> parlower,
                  bool printIter, Nullable<CharacterVector> traceFile) {

  const int K = parinit.size();
  CharacterVector parnames = parinit.names();

  std::vector<double> pl(K, -kInf), pu(K, kInf), ps(K, 1.0);
  fill_bound(parupper, pu, parnames, K);
  fill_bound(parlower, pl, parnames, K);
  const bool rescale = parscale.isNotNull();
  fill_parscale(parscale, ps, K, "trustL1");

  std::vector<double> theta(K);
  bool any_above = false, any_below = false;
  for (int i = 0; i < K; ++i) {
    if (parinit[i] > pu[i])      { any_above = true; theta[i] = pu[i]; }
    else if (parinit[i] < pl[i]) { any_below = true; theta[i] = pl[i]; }
    else                          theta[i] = parinit[i];
  }
  if (any_above) Rf_warning("init above range");
  if (any_below) Rf_warning("init below range");

  std::vector<unsigned char> at_upper(K, 0), at_lower(K, 0);
  for (int i = 0; i < K; ++i) {
    at_upper[i] = (theta[i] >= pu[i]) ? 1 : 0;
    at_lower[i] = (theta[i] <= pl[i]) ? 1 : 0;
  }

  Reporter report;
  report.init(printIter, iterlim, traceFile, K, parnames);

  NumericVector x_named(theta.begin(), theta.end());
  x_named.names() = parnames;
  List out_init = as<List>(objfun(x_named));
  double val_obj = as<double>(out_init["value"]);
  NumericVector grad0 = as<NumericVector>(out_init["gradient"]);
  NumericMatrix Hmat0 = as<NumericMatrix>(out_init["hessian"]);
  if (!std::isfinite(val_obj)) stop("parinit not feasible: value is not finite");

  std::vector<double> grad_obj(grad0.begin(), grad0.end());
  std::vector<double> H_full((std::size_t) K * K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i)
      H_full[i + (std::size_t) j * K] = Hmat0(i, j);

  double val = val_obj + l1.value(theta);
  int neval = 1;
  report(neval, val, x_named, /*head=*/true);

  Blather trace;
  std::vector<int>    active;
  std::vector<double> g_red, H_red, eigvals_red, eigvecs_red;
  double f_used = val;

  bool accept = true, converged = false, is_terminate = false;
  double r = rinit;
  int iter = 0, n_fail = 0;
  std::vector<double> theta_try(K), p_full(K), p_red;

  for (iter = 1; iter <= iterlim; ++iter) {

    if (blather_on) {
      trace.argpath.insert(trace.argpath.end(), theta.begin(), theta.end());
      trace.r.push_back(r);
      trace.valpath.push_back(val);
    }

    if (accept) {
      active.clear();
      for (int i = 0; i < K; ++i) {
        double gi = grad_obj[i];
        bool drop_bound;
        if (minimize) drop_bound = (at_upper[i] && gi < 0.0) || (at_lower[i] && gi > 0.0);
        else          drop_bound = (at_upper[i] && gi > 0.0) || (at_lower[i] && gi < 0.0);
        if (!drop_bound && !l1.pinned(i, theta[i], gi)) active.push_back(i);
      }
      int Kred = static_cast<int>(active.size());
      g_red.assign(Kred, 0.0);
      H_red.assign((std::size_t) Kred * Kred, 0.0);
      for (int ii = 0; ii < Kred; ++ii) {
        int i = active[ii];
        g_red[ii] = grad_obj[i] + l1.grad(i, theta[i]);
        for (int jj = 0; jj < Kred; ++jj) {
          int j = active[jj];
          H_red[ii + (std::size_t) jj * Kred] = H_full[i + (std::size_t) j * K];
        }
      }
      if (rescale) {
        for (int ii = 0; ii < Kred; ++ii) g_red[ii] /= ps[active[ii]];
        for (int jj = 0; jj < Kred; ++jj)
          for (int ii = 0; ii < Kred; ++ii)
            H_red[ii + (std::size_t) jj * Kred] /= ps[active[ii]] * ps[active[jj]];
      }
      f_used = val;
      if (!minimize) {
        for (auto& x : g_red) x = -x;
        for (auto& x : H_red) x = -x;
        f_used = -val;
      }
      if (Kred > 0) {
        eigvals_red.assign(Kred, 0.0);
        eigvecs_red.assign((std::size_t) Kred * Kred, 0.0);
        eigen_sym_local(H_red.data(), Kred, eigvals_red.data(), eigvecs_red.data());
      }
    }

    int Kred = static_cast<int>(active.size());

    p_red.assign(Kred, 0.0);
    bool is_newton = false, is_hard = false, is_easy = false;
    double m_value = 0.0, pred_pos = 0.0;
    if (Kred > 0) {
      trust_sub(Kred, g_red.data(), eigvals_red.data(), eigvecs_red.data(), r,
                p_red.data(), &pred_pos, &is_newton, &is_hard, &is_easy);
      m_value  = model_value(Kred, g_red.data(), H_red.data(), p_red.data());
      pred_pos = -m_value;
    }

    std::fill(p_full.begin(), p_full.end(), 0.0);
    for (int ii = 0; ii < Kred; ++ii) {
      int i = active[ii];
      double s = p_red[ii];
      if (rescale) s /= ps[i];
      p_full[i] = s;
    }
    double stepnorm = 0.0;
    for (int i = 0; i < K; ++i) stepnorm += p_full[i] * p_full[i];
    stepnorm = std::sqrt(stepnorm);

    for (int i = 0; i < K; ++i) theta_try[i] = theta[i] + p_full[i];
    l1.clamp(theta, theta_try);

    for (int i = 0; i < K; ++i) {
      at_upper[i] = !(theta_try[i] < pu[i]) ? 1 : 0;
      at_lower[i] = !(theta_try[i] > pl[i]) ? 1 : 0;
      if (at_upper[i]) theta_try[i] = pu[i];
      if (at_lower[i]) theta_try[i] = pl[i];
    }

    NumericVector x_try(theta_try.begin(), theta_try.end());
    x_try.names() = parnames;
    List out_try;
    bool eval_ok = eval_objfun(objfun, x_try, out_try);
    double val_obj_try = kInf, val_try = kInf;
    NumericVector grad_try;
    NumericMatrix Htry_mat;
    if (eval_ok) {
      val_obj_try = as<double>(out_try["value"]);
      grad_try    = as<NumericVector>(out_try["gradient"]);
      Htry_mat    = as<NumericMatrix>(out_try["hessian"]);
      if (!std::isfinite(val_obj_try)) eval_ok = false;
      else val_try = val_obj_try + l1.value(theta_try);
    }
    neval++;
    report(neval, val_try, x_try, /*head=*/false);

    double ftry_used = minimize ? val_try : -val_try;
    double rho;
    if (eval_ok && pred_pos > 0.0) rho = (f_used - ftry_used) / pred_pos;
    else                           rho = -kInf;

    is_terminate = eval_ok && (std::fabs(ftry_used - f_used) < fterm ||
                               std::fabs(m_value) < mterm);

    if (!eval_ok) {
      n_fail++;
      accept = false;
      r *= 0.25;
      if (n_fail >= 3) {
        Rf_warning("trustL1: objfun evaluation failed 3 times in a row");
        break;
      }
    } else {
      n_fail = 0;
      if (is_terminate) {
        accept = (ftry_used < f_used);
      } else if (rho < 0.25) {
        accept = false;
        r *= 0.25;
      } else {
        accept = true;
        if (rho > 0.75 && !is_newton) r = std::min(2.0 * r, rmax);
      }
    }

    if (accept && eval_ok) {
      theta = theta_try;
      val   = val_try;
      grad_obj.assign(grad_try.begin(), grad_try.end());
      for (int j = 0; j < K; ++j)
        for (int i = 0; i < K; ++i)
          H_full[i + (std::size_t) j * K] = Htry_mat(i, j);
    }

    if (blather_on) {
      trace.argtry.insert(trace.argtry.end(), theta_try.begin(), theta_try.end());
      trace.valtry.push_back(val_try);
      trace.accept.push_back(accept ? 1 : 0);
      trace.preddiff.push_back(m_value);
      trace.stepnorm.push_back(stepnorm);
      trace.rho.push_back(rho);
      trace.steptype.push_back(subproblem_label(is_newton, is_hard, is_easy));
      trace.stepback.push_back("full");
    }

    if (is_terminate) { converged = true; break; }
  }

  int final_iter = (iter <= iterlim) ? iter : iterlim;
  if (!converged && final_iter == iterlim)
    Rf_warning("Maximum number of iterations exceeded. Fit is not converged.");

  NumericVector arg_out(theta.begin(), theta.end());
  arg_out.names() = parnames;
  NumericVector grad_out(K);
  for (int i = 0; i < K; ++i) grad_out[i] = grad_obj[i] + l1.grad(i, theta[i]);
  grad_out.names() = parnames;
  NumericMatrix Hess_out(K, K);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) Hess_out(i, j) = H_full[i + (std::size_t) j * K];
  Hess_out.attr("dimnames") = List::create(parnames, parnames);
  LogicalVector at_bound_out(K);
  for (int i = 0; i < K; ++i) at_bound_out[i] = (at_upper[i] || at_lower[i]);
  at_bound_out.names() = parnames;

  List result = List::create(
      Named("argument")   = arg_out,
      Named("value")      = val,
      Named("gradient")   = grad_out,
      Named("hessian")    = Hess_out,
      Named("iterations") = final_iter,
      Named("converged")  = converged,
      Named("atBound")    = at_bound_out,
      Named("stopReason") = std::string(converged ? "fvalue" : "iterlim"));

  if (blather_on) trace.attach(result, final_iter, K, parnames, minimize);
  return result;
}

}  // namespace


// [[Rcpp::export]]
List trustL1_impl(Function objfun,
                  NumericVector parinit,
                  NumericVector mu,
                  NumericVector lambda,
                  bool   one_sided,
                  double rinit,
                  double rmax,
                  Nullable<NumericVector> parscale = R_NilValue,
                  int    iterlim   = 100,
                  double ftol      = 1e-6,
                  double mtol      = 1e-6,
                  double gtol      = 1e-6,
                  double xtol      = 0.0,
                  double rmin      = 0.0,
                  double thetamax  = 0.99995,
                  std::string boundary = "reflective",
                  bool   minimize  = true,
                  bool   blather   = false,
                  Nullable<NumericVector>  parupper  = R_NilValue,
                  Nullable<NumericVector>  parlower  = R_NilValue,
                  bool   printIter = false,
                  Nullable<CharacterVector> traceFile = R_NilValue) {

  const int K = parinit.size();
  if (K == 0) stop("trustL1: parinit must be non-empty");
  if (!parinit.hasAttribute("names"))
    stop("trustL1: parinit must be a named numeric vector");
  for (int i = 0; i < K; ++i)
    if (!std::isfinite(parinit[i])) stop("trustL1: parinit not all finite");

  L1Spec l1;
  l1.build(mu, lambda, parinit.names(), K, one_sided);

  if (boundary == "clip")
    return trustL1_clip(objfun, parinit, l1, rinit, rmax, parscale, iterlim,
                        ftol, mtol, minimize, blather,
                        parupper, parlower, printIter, traceFile);
  if (boundary != "reflective")
    stop("trustL1: boundary must be one of \"reflective\", \"clip\"");

  return trustL1_reflective(objfun, parinit, l1, rinit, rmax, parscale, iterlim,
                            ftol, mtol, gtol, xtol, rmin, thetamax,
                            minimize, blather, parupper, parlower,
                            printIter, traceFile);
}
