#' Nonlinear Least Squares Optimization with Analytical Jacobian
#'
#' Solves nonlinear least squares problems using various optimization methods.
#' Requires full Jacobian matrix.
#'
#' @param obj Objective function that takes a named numeric vector of parameters
#'   and returns a list with components \code{residuals} and \code{jacobian}.
#'   The function signature should be:
#'   \code{function(par) { ... return(list(residuals = ..., jacobian = ...)) }}
#'   
#'   The Jacobian MUST be provided as a matrix with dimensions
#'   [n_residuals x n_parameters]. All derivatives must be computed analytically.
#'
#' @param init Named numeric vector of initial parameter values. Parameter names
#'   are used for output. If unnamed, default names (p1, p2, ...) are assigned
#'   with a warning.
#'
#' @param method Character string specifying the optimization method:
#'   \itemize{
#'     \item \code{"lm"} (default): Levenberg-Marquardt algorithm. Most robust,
#'           handles ill-conditioned problems well, automatically adjusts between
#'           gradient descent and Gauss-Newton.
#'     \item \code{"gn"}: Gauss-Newton algorithm. Faster convergence than LM but
#'           requires a good initial guess. May fail for ill-conditioned problems.
#'     \item \code{"gd"}: Gradient Descent. Simplest method, guaranteed to decrease
#'           error at each step but slowest convergence.
#'   }
#' @param solver Character string specifying the linear solver:
#'   \itemize{
#'     \item \code{"cholesky"} (default): Fast Cholesky decomposition for 
#'           positive-definite systems. Recommended for most problems (3-5x faster than SVD).
#'     \item \code{"svd"}: Robust SVD-based solver. Use for ill-conditioned or
#'           rank-deficient problems. Slower but handles singular matrices.
#'   }
#' @param control Named list of control parameters:
#'   \describe{
#'     \item{\code{max_iter}}{Maximum number of iterations (default: 100).
#'           Set to 0 for unlimited iterations.}
#'     \item{\code{min_step_len}}{Minimum step length for convergence (default: sqrt(.Machine$double.eps)).
#'           Algorithm stops when ||step|| < min_step_len.}
#'     \item{\code{min_grad_len}}{Minimum gradient length for convergence (default: sqrt(.Machine$double.eps)).
#'           Algorithm stops when ||gradient|| < min_grad_len.}
#'     \item{\code{min_error}}{Minimum error for convergence (default: 0).
#'           Algorithm stops when error < min_error.}
#'     \item{\code{min_error_delta}}{Minimum error change for stagnation detection
#'           (default: sqrt(.Machine$double.eps)). If |error_new - error_old| < min_error_delta for
#'           max_stagnation_iter consecutive iterations, optimization stops.}
#'     \item{\code{max_stagnation_iter}}{Maximum consecutive stagnating iterations
#'           (default: 2). Set to 0 or negative to disable stagnation detection.}
#'     \item{\code{verbosity}}{Output level (default: 0):
#'           \itemize{
#'             \item 0: Silent
#'             \item 1: Print iteration, step length, gradient length, error
#'             \item 2+: More detailed output
#'           }}
#'     \item{\code{refine_method}}{Step refinement method (default: "constant"):
#'           \itemize{
#'             \item \code{"constant"}: Fixed step size
#'             \item \code{"armijo"}: Armijo backtracking line search
#'             \item \code{"wolfe"}: Wolfe line search (Armijo + curvature condition)
#'             \item \code{"dogleg"}: Powell's dogleg method (trust region)
#'             \item \code{"barzilai_borwein"}: Barzilai-Borwein step size
#'           }
#'           Note: Only for Gauss-Newton and Gradient Descent. LM has built-in control.}
#'     \item{\code{step_factor}}{Step scaling factor for constant refinement (default: 1.0).
#'           Values < 1 make the algorithm more conservative.}
#'     \item{\code{lm_lambda}}{Initial damping parameter for LM (default: 1.0).}
#'     \item{\code{lm_lambda_inc}}{Lambda increase factor for LM (default: 2.0).}
#'     \item{\code{lm_lambda_dec}}{Lambda decrease factor for LM (default: 0.5).}
#'     \item{\code{lm_max_iter}}{Maximum inner iterations for LM lambda search (default: 0 = unlimited).}
#'     \item{\code{armijo_decrease}}{Backtracking decrease factor (default: 0.8).}
#'     \item{\code{armijo_c1}}{Armijo condition constant (default: 1e-4).}
#'     \item{\code{armijo_min_step}}{Minimum step bound for Armijo (default: 1e-10).}
#'     \item{\code{armijo_max_step}}{Maximum step bound for Armijo (default: 1.0).}
#'     \item{\code{armijo_max_iter}}{Maximum line search iterations (default: 0 = unlimited).}
#'     \item{\code{wolfe_decrease}}{Backtracking decrease factor (default: 0.8).}
#'     \item{\code{wolfe_c1}}{Wolfe Armijo constant (default: 1e-4).}
#'     \item{\code{wolfe_c2}}{Wolfe curvature constant (default: 0.9).}
#'     \item{\code{wolfe_min_step}}{Minimum step bound for Wolfe (default: 1e-10).}
#'     \item{\code{wolfe_max_step}}{Maximum step bound for Wolfe (default: 1.0).}
#'     \item{\code{wolfe_max_iter}}{Maximum line search iterations (default: 0 = unlimited).}
#'     \item{\code{dogleg_radius}}{Initial trust region radius (default: 1.0).}
#'     \item{\code{dogleg_max_radius}}{Maximum trust region radius (default: 2.0).}
#'     \item{\code{dogleg_radius_eps}}{Radius epsilon for increase trigger (default: 1e-6).}
#'     \item{\code{dogleg_accept_fitness}}{Minimum fitness for acceptance (default: 0.25).}
#'     \item{\code{dogleg_max_iter}}{Maximum trust region iterations (default: 0 = unlimited).}
#'     \item{\code{bb_mode}}{Barzilai-Borwein mode: "direct" or "inverse" (default: "direct").}
#'     \item{\code{bb_const_step}}{Initial constant step for BB (default: 1e-2).}
#'   }
#'
#' @return Object of class \code{lsqnl} (a list) with components:
#'   \describe{
#'     \item{\code{par}}{Named numeric vector of optimized parameter values}
#'     \item{\code{residuals}}{Final residual vector at the optimum}
#'     \item{\code{error}}{Final error value: 0.5 * sum(residuals^2)}
#'     \item{\code{iterations}}{Number of iterations performed}
#'     \item{\code{converged}}{Logical: TRUE if a convergence criterion was met
#'           (step length, gradient length, error threshold, or stagnation)}
#'     \item{\code{succeeded}}{Logical: TRUE if optimization succeeded without
#'           numerical issues (e.g., singular matrices)}
#'   }
#'
#' @details
#' \strong{Analytical Jacobian Required:}
#'
#' This function requires a complete analytical Jacobian matrix. The objective
#' function must return a list with:
#' \itemize{
#'   \item \code{residuals}: Numeric vector of length n_residuals
#'   \item \code{jacobian}: Numeric matrix of dimension [n_residuals x n_parameters]
#' }
#'
#' The Jacobian is checked at the first call and must have correct dimensions.
#'
#' \strong{Algorithm Selection:}
#'
#' \itemize{
#'   \item \strong{Levenberg-Marquardt} (default): Best general-purpose choice.
#'         Handles ill-conditioning, automatically adapts step size, most robust.
#'   \item \strong{Gauss-Newton}: Fastest for well-conditioned problems with good
#'         initial guess. Use with Armijo line search or Dogleg for robustness.
#'   \item \strong{Gradient Descent}: Use only for very ill-conditioned problems
#'         or as a fallback. Very slow but guaranteed descent direction.
#' }
#'
#' \strong{Convergence Criteria:}
#'
#' Optimization stops when ANY of these conditions is met:
#' \itemize{
#'   \item Maximum iterations reached
#'   \item Step length < min_step_len (local minimum found)
#'   \item Gradient length < min_grad_len (stationary point reached)
#'   \item Error < min_error (target error achieved)
#'   \item Stagnation detected (error unchanged for max_stagnation_iter iterations)
#' }
#'
#' \strong{Stagnation Detection:}
#'
#' The optimizer tracks error changes between iterations. If the absolute change
#' |error_new - error_old| < min_error_delta for max_stagnation_iter consecutive
#' iterations, the optimization terminates. This is more robust than checking
#' only once, as it avoids premature termination due to numerical noise.
#'
#' @examples
#' \dontrun{
#' # Example 1: Exponential decay with full Jacobian
#' set.seed(42)
#' x <- seq(0, 5, length.out = 50)
#' y <- 2.5 * exp(-0.5 * x) + rnorm(50, sd = 0.1)
#'
#' obj <- function(par) {
#'   a <- par["a"]
#'   b <- par["b"]
#'   exp_term <- exp(-b * x)
#'   residuals <- y - a * exp_term
#'   
#'   # Full analytical Jacobian
#'   jacobian <- cbind(
#'     -exp_term,           # d(residuals)/da
#'     a * x * exp_term     # d(residuals)/db
#'   )
#'   
#'   list(residuals = residuals, jacobian = jacobian)
#' }
#'
#' result <- lsqnl(obj, init = c(a = 1, b = 1))
#' print(result$par)
#'
#' # Example 2: Gauss-Newton with Dogleg trust region
#' result_gn <- lsqnl(obj, c(a = 1, b = 1), method = "gn",
#'                    control = list(refine_method = "dogleg",
#'                                   dogleg_radius = 0.5,
#'                                   dogleg_max_radius = 5.0))
#'
#' # Example 3: Custom LM with tight stagnation control
#' result_lm <- lsqnl(obj, c(a = 1, b = 1), method = "lm",
#'                    control = list(lm_lambda = 0.1, 
#'                                   lm_lambda_inc = 3,
#'                                   min_error_delta = 1e-15,
#'                                   max_stagnation_iter = 20))
#'
#' # Example 4: Verbose output
#' result_verbose <- lsqnl(obj, c(a = 1, b = 1), 
#'                         control = list(verbosity = 2))
#' }
#'
#' @export
lsqnl <- function(obj, init, method = c("lm", "gn", "gd"), 
                  solver = c("cholesky", "svd"), control = list()) {
  
  if(!is.function(obj)) stop("'obj' must be a function")
  if(!is.numeric(init)) stop("'init' must be numeric")
  
  if(is.null(names(init))) {
    names(init) <- paste0("p", seq_along(init))
    warning("'init' has no names. Using default names: p1, p2, ...")
  }
  
  method <- match.arg(method)
  solver <- match.arg(solver)
  
  # Default control parameters
  control <- modifyList(
    list(
      # Convergence criteria
      max_iter = 100,
      min_step_len = 1e-4,
      min_grad_len = 1e-4,
      min_error = 0,
      min_error_delta = 1e-4,
      max_stagnation_iter = 2,
      verbosity = 0,
      
      # Step refinement
      refine_method = "constant",
      step_factor = 1.0,
      
      # Levenberg-Marquardt
      lm_lambda = 1.0,
      lm_lambda_inc = 2.0,
      lm_lambda_dec = 0.5,
      lm_max_iter = 0,
      
      # Armijo backtracking
      armijo_decrease = 0.8,
      armijo_c1 = 1e-4,
      armijo_min_step = 1e-10,
      armijo_max_step = 1.0,
      armijo_max_iter = 0,
      
      # Wolfe line search
      wolfe_decrease = 0.8,
      wolfe_c1 = 1e-4,
      wolfe_c2 = 0.9,
      wolfe_min_step = 1e-10,
      wolfe_max_step = 1.0,
      wolfe_max_iter = 0,
      
      # Dogleg trust region
      dogleg_radius = 1.0,
      dogleg_max_radius = 2.0,
      dogleg_radius_eps = 1e-6,
      dogleg_accept_fitness = 0.25,
      dogleg_max_iter = 0,
      
      # Barzilai-Borwein
      bb_mode = "direct",
      bb_const_step = 1e-2
    ),
    control
  )
  
  # Test objective function
  test <- tryCatch(obj(init), error = function(e) 
    stop("Error evaluating objective function: ", e$message))
  
  if(!is.list(test)) 
    stop("Objective function must return a list")
  if(!("residuals" %in% names(test))) 
    stop("Objective function must return 'residuals' component")
  if(!("jacobian" %in% names(test))) 
    stop("Objective function must return 'jacobian' component")
  if(!is.numeric(test$residuals)) 
    stop("Residuals must be numeric")
  if(!is.matrix(test$jacobian) && !is.null(test$jacobian))
    stop("Jacobian must be a matrix or NULL")
  if(is.null(test$jacobian))
    stop("Jacobian cannot be NULL. Analytical Jacobian is required.")
  if(ncol(test$jacobian) != length(init))
    stop("Jacobian must have ", length(init), " columns (one per parameter)")
  if(nrow(test$jacobian) != length(test$residuals))
    stop("Jacobian must have ", length(test$residuals), " rows (one per residual)")
  
  # Call C++ implementation
  result <- lsqnl_cpp(
    obj = obj,
    init = init,
    method = method,
    solver = solver,
    max_iter = control$max_iter,
    min_step_len = control$min_step_len,
    min_grad_len = control$min_grad_len,
    min_error = control$min_error,
    min_error_delta = control$min_error_delta,
    max_stagnation_iter = control$max_stagnation_iter,
    verbosity = control$verbosity,
    refine_method = control$refine_method,
    step_factor = control$step_factor,
    lm_lambda = control$lm_lambda,
    lm_lambda_inc = control$lm_lambda_inc,
    lm_lambda_dec = control$lm_lambda_dec,
    lm_max_iter = control$lm_max_iter,
    armijo_decrease = control$armijo_decrease,
    armijo_c1 = control$armijo_c1,
    armijo_min_step = control$armijo_min_step,
    armijo_max_step = control$armijo_max_step,
    armijo_max_iter = control$armijo_max_iter,
    wolfe_decrease = control$wolfe_decrease,
    wolfe_c1 = control$wolfe_c1,
    wolfe_c2 = control$wolfe_c2,
    wolfe_min_step = control$wolfe_min_step,
    wolfe_max_step = control$wolfe_max_step,
    wolfe_max_iter = control$wolfe_max_iter,
    dogleg_radius = control$dogleg_radius,
    dogleg_max_radius = control$dogleg_max_radius,
    dogleg_radius_eps = control$dogleg_radius_eps,
    dogleg_accept_fitness = control$dogleg_accept_fitness,
    dogleg_max_iter = control$dogleg_max_iter,
    bb_mode = control$bb_mode,
    bb_const_step = control$bb_const_step
  )
  
  class(result) <- c("lsqnl", "list")
  result
}

#' Print method for lsqnl objects
#'
#' @param x An object of class \code{lsqnl}
#' @param digits Number of significant digits to print (default: 6)
#' @param ... Additional arguments (ignored)
#'
#' @return The input object (invisibly)
#'
#' @export
print.lsqnl <- function(x, digits = 6, ...) {
  cat("Nonlinear Least Squares Result\n")
  cat("==============================\n\n")
  cat("Status:", if(x$converged) "Converged" else "Not converged", "\n")
  cat("       ", if(x$succeeded) "Succeeded" else "Failed", "\n\n")
  
  cat("Parameters:\n")
  par_df <- data.frame(
    Value = signif(x$par, digits),
    row.names = names(x$par)
  )
  print(par_df)
  
  cat("\n")
  cat("Error:      ", signif(x$error, digits), 
      " (0.5 * sum(residuals^2))\n")
  cat("RSS:        ", signif(sum(x$residuals^2), digits), 
      " (residual sum of squares)\n")
  cat("RMSE:       ", signif(sqrt(mean(x$residuals^2)), digits), 
      " (root mean squared error)\n")
  cat("Iterations: ", x$iterations, "\n")
  
  invisible(x)
}

#' Summary method for lsqnl objects
#'
#' Provides detailed summary statistics for the optimization result.
#'
#' @param object An object of class \code{lsqnl}
#' @param ... Additional arguments (ignored)
#'
#' @return Object of class \code{summary.lsqnl} with additional summary statistics
#'
#' @export
summary.lsqnl <- function(object, ...) {
  structure(
    list(
      par = object$par,
      error = object$error,
      rss = sum(object$residuals^2),
      rmse = sqrt(mean(object$residuals^2)),
      n_residuals = length(object$residuals),
      n_parameters = length(object$par),
      iterations = object$iterations,
      converged = object$converged,
      succeeded = object$succeeded,
      residuals = object$residuals
    ),
    class = "summary.lsqnl"
  )
}

#' Print method for summary.lsqnl objects
#'
#' @param x An object of class \code{summary.lsqnl}
#' @param digits Number of significant digits (default: 6)
#' @param ... Additional arguments (ignored)
#'
#' @return The input object (invisibly)
#'
#' @export
print.summary.lsqnl <- function(x, digits = 6, ...) {
  cat("Nonlinear Least Squares Summary\n")
  cat("===============================\n\n")
  
  cat("Optimization Status:\n")
  cat("  Converged: ", x$converged, "\n")
  cat("  Succeeded: ", x$succeeded, "\n")
  cat("  Iterations:", x$iterations, "\n\n")
  
  cat("Problem Dimensions:\n")
  cat("  Parameters: ", x$n_parameters, "\n")
  cat("  Residuals:  ", x$n_residuals, "\n\n")
  
  cat("Parameters:\n")
  par_df <- data.frame(
    Estimate = signif(x$par, digits),
    row.names = names(x$par)
  )
  print(par_df)
  
  cat("\n")
  cat("Fit Statistics:\n")
  cat("  Error (0.5*RSS): ", signif(x$error, digits), "\n")
  cat("  RSS:             ", signif(x$rss, digits), "\n")
  cat("  RMSE:            ", signif(x$rmse, digits), "\n")
  cat("  Degrees of freedom:", x$n_residuals - x$n_parameters, "\n\n")
  
  cat("Residual Statistics:\n")
  res_summary <- signif(summary(x$residuals), digits)
  cat("  Min:    ", res_summary[1], "\n")
  cat("  1st Qu: ", res_summary[2], "\n")
  cat("  Median: ", res_summary[3], "\n")
  cat("  Mean:   ", res_summary[4], "\n")
  cat("  3rd Qu: ", res_summary[5], "\n")
  cat("  Max:    ", res_summary[6], "\n")
  
  invisible(x)
}