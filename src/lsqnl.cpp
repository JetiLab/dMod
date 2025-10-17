// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <lsqcpp.hpp>

// Simplified Objective Function - Always expects full Jacobian
class RObjectiveFunctionAnalytical
{
public:
  static constexpr bool ComputesJacobian = true;
  
  // Default constructor (required for lsqcpp initialization)
  RObjectiveFunctionAnalytical() 
    : _obj_fn(R_NilValue), _param_names()
  {}
  
  RObjectiveFunctionAnalytical(SEXP obj_fn, Rcpp::CharacterVector param_names)
    : _obj_fn(obj_fn), _param_names(param_names)
  {}
  
  // 3-parameter version: ONLY this is called by lsqcpp
  template<typename Input, typename Output, typename Jacobian>
  void operator()(const Eigen::MatrixBase<Input> &xval,
                Eigen::MatrixBase<Output> &fval,
                Eigen::MatrixBase<Jacobian> &jacobian) const
  {
    // Check if initialized
    if(_obj_fn == R_NilValue) {
      Rcpp::stop("Objective function not initialized!");
    }
    
    // Convert to R
    Rcpp::NumericVector x_r = Rcpp::wrap(xval);
    x_r.names() = _param_names;
    
    // Call R function
    Rcpp::Function fn(_obj_fn);
    Rcpp::List result = fn(x_r);
    
    // Extract residuals
    Rcpp::NumericVector residuals = result["residuals"];
    int n_res = residuals.size();
    int n_par = xval.size();
    
    fval.derived().resize(n_res);
    Eigen::Map<const Eigen::VectorXd> res_map(residuals.begin(), n_res);
    fval.derived() = res_map;
    
    // Extract Jacobian - MUST be provided!
    if(!result.containsElementNamed("jacobian") || Rf_isNull(result["jacobian"])) {
      Rcpp::stop("Jacobian must be provided by objective function!");
    }
    
    Rcpp::NumericMatrix jac_r = result["jacobian"];
    
    // Validate dimensions
    if(jac_r.nrow() != n_res || jac_r.ncol() != n_par) {
      Rcpp::stop("Jacobian dimensions mismatch. Expected " + 
        std::to_string(n_res) + "x" + std::to_string(n_par) + 
        ", got " + std::to_string(jac_r.nrow()) + "x" + std::to_string(jac_r.ncol()));
    }
    
    // Direct memory copy - maximum efficiency!
    jacobian.derived().resize(n_res, n_par);
    std::memcpy(jacobian.derived().data(), jac_r.begin(), n_res * n_par * sizeof(double));
  }
  
private:
  SEXP _obj_fn;
  Rcpp::CharacterVector _param_names;
};

//' Nonlinear Least Squares with Analytical Jacobian
//'
//' High-performance C++ implementation of nonlinear least squares optimization
//' using lsqcpp library. Requires analytical Jacobian matrix.
//'
//' @param obj R function that returns list(residuals = ..., jacobian = ...)
//' @param init Named numeric vector of initial parameter values
//' @param method Optimization method: "lm", "gn", "gd"
//' @param solver Linear solver: "cholesky" (fast, default) or "svd" (robust)
//' @param max_iter Maximum iterations (0 = unlimited)
//' @param min_step_len Minimum step length for convergence
//' @param min_grad_len Minimum gradient length for convergence  
//' @param min_error Minimum error for convergence
//' @param min_error_delta Minimum error change for stagnation detection
//' @param max_stagnation_iter Maximum stagnating iterations (0 = disabled)
//' @param verbosity Output level (0 = silent)
//' @param refine_method Step refinement: "constant", "armijo", "wolfe", "dogleg", "barzilai_borwein"
//' @param step_factor Step factor for constant refinement
//' @param lm_lambda Initial lambda for Levenberg-Marquardt
//' @param lm_lambda_inc Lambda increase factor
//' @param lm_lambda_dec Lambda decrease factor
//' @param lm_max_iter Maximum LM lambda search iterations
//' @param armijo_decrease Backtracking decrease factor
//' @param armijo_c1 Armijo condition constant
//' @param armijo_min_step Minimum step bound
//' @param armijo_max_step Maximum step bound
//' @param armijo_max_iter Maximum line search iterations
//' @param wolfe_decrease Backtracking decrease factor
//' @param wolfe_c1 Wolfe Armijo constant
//' @param wolfe_c2 Wolfe curvature constant
//' @param wolfe_min_step Minimum step bound
//' @param wolfe_max_step Maximum step bound
//' @param wolfe_max_iter Maximum line search iterations
//' @param dogleg_radius Initial trust region radius
//' @param dogleg_max_radius Maximum trust region radius
//' @param dogleg_radius_eps Radius epsilon for increase trigger
//' @param dogleg_accept_fitness Minimum fitness for acceptance
//' @param dogleg_max_iter Maximum trust region iterations
//' @param bb_mode Barzilai-Borwein mode: "direct" or "inverse"
//' @param bb_const_step Initial constant step
//'
//' @return List with: par, residuals, error, iterations, converged, succeeded
//' 
//' @details
//' This function provides a high-performance interface to the lsqcpp C++ library.
//' It requires full analytical Jacobian matrices for maximum efficiency.
//' 
//' Solver selection:
//' - "cholesky": Fast Cholesky decomposition (3-5x faster than SVD). Use for 
//'   well-conditioned problems (default and recommended).
//' - "svd": Robust SVD-based solver. Use only for ill-conditioned or rank-deficient problems.
//' 
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List lsqnl_cpp(
   SEXP obj,
   Rcpp::NumericVector init,
   std::string method = "lm",
   std::string solver = "cholesky",
   int max_iter = 100,
   double min_step_len = 1e-9,
   double min_grad_len = 1e-9,
   double min_error = 0.0,
   double min_error_delta = 1e-12,
   int max_stagnation_iter = 10,
   int verbosity = 0,
   std::string refine_method = "constant",
   double step_factor = 1.0,
   double lm_lambda = 1.0,
   double lm_lambda_inc = 2.0,
   double lm_lambda_dec = 0.5,
   int lm_max_iter = 0,
   double armijo_decrease = 0.8,
   double armijo_c1 = 1e-4,
   double armijo_min_step = 1e-10,
   double armijo_max_step = 1.0,
   int armijo_max_iter = 0,
   double wolfe_decrease = 0.8,
   double wolfe_c1 = 1e-4,
   double wolfe_c2 = 0.9,
   double wolfe_min_step = 1e-10,
   double wolfe_max_step = 1.0,
   int wolfe_max_iter = 0,
   double dogleg_radius = 1.0,
   double dogleg_max_radius = 2.0,
   double dogleg_radius_eps = 1e-6,
   double dogleg_accept_fitness = 0.25,
   int dogleg_max_iter = 0,
   std::string bb_mode = "direct",
   double bb_const_step = 1e-2)
{
 using Scalar = double;
 
 // Handle parameter names
 Rcpp::CharacterVector param_names = init.names();
 if(param_names.size() == 0) {
   param_names = Rcpp::CharacterVector(init.size());
   for(int i = 0; i < init.size(); ++i) {
     param_names[i] = "p" + std::to_string(i + 1);
   }
 }
 
 // Initialize
 Eigen::VectorXd init_eigen = Rcpp::as<Eigen::VectorXd>(init);
 auto objective = RObjectiveFunctionAnalytical(obj, param_names);
 
 // Result variables
 Eigen::VectorXd result_vec;
 Eigen::VectorXd residuals;
 double error;
 int iterations;
 bool converged, succeeded;
 
 // Common algorithm setup lambda
 auto setup_common = [&](auto& algo) {
   algo.setObjective(objective);
   algo.setMaximumIterations(max_iter);
   algo.setMinimumStepLength(min_step_len);
   algo.setMinimumGradientLength(min_grad_len);
   algo.setMinimumError(min_error);
   algo.setMinimumErrorDelta(min_error_delta);
   algo.setMaximumStagnationIterations(max_stagnation_iter);
   algo.setVerbosity(verbosity);
 };
 
 auto extract_result = [&](auto& result) {
   result_vec = result.xval;
   residuals = result.fval;
   error = result.error;
   iterations = result.iterations;
   converged = result.converged;
   succeeded = result.succeeded;
 };
 
 // Solver selection
 bool use_cholesky = (solver == "cholesky");
 
 // ============ LEVENBERG-MARQUARDT ============
 if(method == "lm") {
   if(use_cholesky) {
     using Algo = lsqcpp::LevenbergMarquardtX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::DenseCholeskySolver>;
     Algo algo;
     lsqcpp::LevenbergMarquardtParameter<Scalar> lm_param(lm_lambda, lm_lambda_inc, lm_lambda_dec, lm_max_iter);
     algo.setMethodParameters(lm_param);
     setup_common(algo);
     auto result = algo.minimize(init_eigen);
     extract_result(result);
   } else {
     using Algo = lsqcpp::LevenbergMarquardtX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::DenseSVDSolver>;
     Algo algo;
     lsqcpp::LevenbergMarquardtParameter<Scalar> lm_param(lm_lambda, lm_lambda_inc, lm_lambda_dec, lm_max_iter);
     algo.setMethodParameters(lm_param);
     setup_common(algo);
     auto result = algo.minimize(init_eigen);
     extract_result(result);
   }
 }
 // ============ GAUSS-NEWTON ============
 else if(method == "gn") {
   
   if(refine_method == "armijo") {
     if(use_cholesky) {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ArmijoBacktracking, lsqcpp::DenseCholeskySolver>;
       Algo algo;
       lsqcpp::ArmijoBacktrackingParameter<Scalar> param(armijo_decrease, armijo_c1, armijo_min_step, armijo_max_step, armijo_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     } else {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ArmijoBacktracking, lsqcpp::DenseSVDSolver>;
       Algo algo;
       lsqcpp::ArmijoBacktrackingParameter<Scalar> param(armijo_decrease, armijo_c1, armijo_min_step, armijo_max_step, armijo_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     }
   }
   else if(refine_method == "wolfe") {
     if(use_cholesky) {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::WolfeBacktracking, lsqcpp::DenseCholeskySolver>;
       Algo algo;
       lsqcpp::WolfeBacktrackingParameter<Scalar> param(wolfe_decrease, wolfe_c1, wolfe_c2, wolfe_min_step, wolfe_max_step, wolfe_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     } else {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::WolfeBacktracking, lsqcpp::DenseSVDSolver>;
       Algo algo;
       lsqcpp::WolfeBacktrackingParameter<Scalar> param(wolfe_decrease, wolfe_c1, wolfe_c2, wolfe_min_step, wolfe_max_step, wolfe_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     }
   }
   else if(refine_method == "dogleg") {
     if(use_cholesky) {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::DoglegMethod, lsqcpp::DenseCholeskySolver>;
       Algo algo;
       lsqcpp::DoglegMethodParameter<Scalar> param(dogleg_radius, dogleg_max_radius, dogleg_radius_eps, dogleg_accept_fitness, dogleg_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     } else {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::DoglegMethod, lsqcpp::DenseSVDSolver>;
       Algo algo;
       lsqcpp::DoglegMethodParameter<Scalar> param(dogleg_radius, dogleg_max_radius, dogleg_radius_eps, dogleg_accept_fitness, dogleg_max_iter);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     }
   }
   else if(refine_method == "barzilai_borwein") {
     if(use_cholesky) {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::BarzilaiBorwein, lsqcpp::DenseCholeskySolver>;
       Algo algo;
       auto mode = bb_mode == "inverse" ? lsqcpp::BarzilaiBorwein::Mode::Inverse : lsqcpp::BarzilaiBorwein::Mode::Direct;
       lsqcpp::BarzilaiBorweinParameter<Scalar> param(mode, bb_const_step);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     } else {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::BarzilaiBorwein, lsqcpp::DenseSVDSolver>;
       Algo algo;
       auto mode = bb_mode == "inverse" ? lsqcpp::BarzilaiBorwein::Mode::Inverse : lsqcpp::BarzilaiBorwein::Mode::Direct;
       lsqcpp::BarzilaiBorweinParameter<Scalar> param(mode, bb_const_step);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     }
   }
   else {
     if(use_cholesky) {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ConstantStepFactor, lsqcpp::DenseCholeskySolver>;
       Algo algo;
       lsqcpp::ConstantStepFactorParameter<Scalar> param(step_factor);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     } else {
       using Algo = lsqcpp::GaussNewtonX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ConstantStepFactor, lsqcpp::DenseSVDSolver>;
       Algo algo;
       lsqcpp::ConstantStepFactorParameter<Scalar> param(step_factor);
       algo.setRefinementParameters(param);
       setup_common(algo);
       auto result = algo.minimize(init_eigen);
       extract_result(result);
     }
   }
 }
 // ============ GRADIENT DESCENT ============
 else if(method == "gd") {
   
   if(refine_method == "armijo") {
     using Algo = lsqcpp::GradientDescentX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ArmijoBacktracking>;
     Algo algo;
     lsqcpp::ArmijoBacktrackingParameter<Scalar> param(armijo_decrease, armijo_c1, armijo_min_step, armijo_max_step, armijo_max_iter);
     algo.setRefinementParameters(param);
     setup_common(algo);
     auto result = algo.minimize(init_eigen);
     extract_result(result);
   }
   else if(refine_method == "barzilai_borwein") {
     using Algo = lsqcpp::GradientDescentX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::BarzilaiBorwein>;
     Algo algo;
     auto mode = bb_mode == "inverse" ? lsqcpp::BarzilaiBorwein::Mode::Inverse : lsqcpp::BarzilaiBorwein::Mode::Direct;
     lsqcpp::BarzilaiBorweinParameter<Scalar> param(mode, bb_const_step);
     algo.setRefinementParameters(param);
     setup_common(algo);
     auto result = algo.minimize(init_eigen);
     extract_result(result);
   }
   else {
     using Algo = lsqcpp::GradientDescentX<Scalar, RObjectiveFunctionAnalytical, lsqcpp::ConstantStepFactor>;
     Algo algo;
     lsqcpp::ConstantStepFactorParameter<Scalar> param(step_factor);
     algo.setRefinementParameters(param);
     setup_common(algo);
     auto result = algo.minimize(init_eigen);
     extract_result(result);
   }
 }
 else {
   Rcpp::stop("Unknown method. Use 'lm', 'gn', or 'gd'");
 }
 
 // Return results
 Rcpp::NumericVector par_out = Rcpp::wrap(result_vec);
 par_out.names() = param_names;
 
 return Rcpp::List::create(
   Rcpp::Named("par") = par_out,
   Rcpp::Named("residuals") = Rcpp::wrap(residuals),
   Rcpp::Named("error") = error,
   Rcpp::Named("iterations") = iterations,
   Rcpp::Named("converged") = converged,
   Rcpp::Named("succeeded") = succeeded
 );
}