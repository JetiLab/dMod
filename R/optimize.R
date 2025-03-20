#' Optimize an Objective Function using Various Optimization Methods using optimx
#'
#' This function extends the existing optimization capabilities by integrating the `optimx` package,
#' allowing the use of various optimization algorithms such as "L-BFGS-B", "BFGS", "Nelder-Mead", etc.
#'
#' @param objfn Object of class obsfn, i.e. a function obj(..., fixed, deriv, conditions, env) that returns an objective list, objlist.
#' @param parinit A numeric vector of initial parameter values.
#' @param method A character string specifying the optimization method. Defaults to "L-BFGS-B".
#'               Available methods include those supported by the `optimx` package, such as "BFGS",
#'               "Nelder-Mead", "L-BFGS-B", etc.
#' @param lower A numeric vector of lower bounds for the parameters (used only by methods that support
#'              box constraints, e.g., "L-BFGS-B"). Defaults to `-Inf`.
#' @param upper A numeric vector of upper bounds for the parameters (used only by methods that support
#'              box constraints, e.g., "L-BFGS-B"). Defaults to `Inf`.
#' @param control A list of control parameters to pass to the optimization algorithm.
#' @param ... Additional arguments to pass to optimx::optimr
#' 
#' @return A list containing:
#'         - `value`: The value of the objective function at the optimum.
#'         - `gradient`: The gradient at the optimum.
#'         - `hessian`: The Hessian at the optimum.
#'         - `argument`: The optimized parameters.
#'         - `converged`: Logical indicating if the optimizer converged.
#'         - `iterations`: The number of function and gradient evaluations.
#' @import optimx
#' @export
optimize <- function(objfn, parinit, method = "L-BFGS-B", lower = -Inf, upper = Inf, control = list(), ...) {
  
  # Sanitize the initial parameters.
  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit <- sanePars$pars
  
  # Ensure lower/upper bounds are vectors with proper names.
  if (length(lower) == 1) {
    lower <- rep(lower, length(parinit))
    names(lower) <- names(parinit)
  } else if (is.null(names(lower))) {
    names(lower) <- names(parinit)
  }
  
  if (length(upper) == 1) {
    upper <- rep(upper, length(parinit))
    names(upper) <- names(parinit)
  } else if (is.null(names(upper))) {
    names(upper) <- names(parinit)
  }
  
  # Initialize cache to store only the most recently evaluated parameter result.
  cache <- new.env(hash = TRUE, parent = emptyenv())
  last_key <- NULL  # Track the key of the last cached parameter set
  
  # Helper function to generate a unique key for each parameter set.
  generate_key <- function(par) paste0(format(par, digits = 8), collapse = ",")
  
  # Define the function wrappers required by optimx with error handling and caching.
  fn <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      
      # Clear previous cache entry if it exists
      if (!is.null(last_key) && exists(last_key, envir = cache)) {
        rm(list = last_key, envir = cache)
      }
      
      # Store only the current result
      assign(key, result, envir = cache)
      last_key <<- key
    }
    
    if (inherits(result, "try-error") || !is.finite(result$value)) {
      stop("Optimization stopped: Objective function evaluation failed.")  # Stop optimization with message
    }
    result$value
  }
  
  gr <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      
      # Clear previous cache entry if it exists
      if (!is.null(last_key) && exists(last_key, envir = cache)) {
        rm(list = last_key, envir = cache)
      }
      
      # Store only the current result
      assign(key, result, envir = cache)
      last_key <<- key
    }
    
    if (inherits(result, "try-error") || !all(is.finite(result$gradient))) {
      stop("Optimization stopped: Gradient evaluation failed.")  # Stop optimization with message
    }
    unname(result$gradient)
  }
  
  hess <- function(par, ...) {
    names(par) <- names(parinit)
    key <- generate_key(par)
    
    if (exists(key, envir = cache)) {
      result <- get(key, envir = cache)
    } else {
      result <- try(objfn(par, ...), silent = TRUE)
      
      # Clear previous cache entry if it exists
      if (!is.null(last_key) && exists(last_key, envir = cache)) {
        rm(list = last_key, envir = cache)
      }
      
      # Store only the current result
      assign(key, result, envir = cache)
      last_key <<- key
    }
    
    if (inherits(result, "try-error") || !all(is.finite(result$hessian))) {
      stop("Optimization stopped: Hessian evaluation failed.")  # Stop optimization with message
    }
    unname(result$hessian)
  }
  
  # Perform the optimization.
  optim_result <- optimx::optimr(
    par = as.numeric(parinit),
    fn = fn,
    gr = gr,
    hess = hess,
    method = method,
    lower = unname(lower),
    upper = unname(upper),
    control = control,
    ...
  )
  
  # Extract the optimized parameters.
  final_par <- structure(optim_result$par, names = names(parinit))
  
  # Remove any unwanted attributes from the optimized parameters
  attributes(final_par) <- list(names = names(parinit))
  
  # Evaluate the objective function at the optimum.
  final_result <- objfn(final_par, ...)
  
  # Attach optimization metadata.
  final_result$argument <- final_par
  final_result$converged <- !as.logical(optim_result$convergence)
  final_result$iterations <- optim_result$counts["function"]
  
  attr(final_result, 'status') <- attr(optim_result$par, 'status') 
  return(final_result)
}

#' Perform Multi-Start Optimization using Various Sampling Methods
#'
#' This function conducts multiple optimization runs starting from different points
#' generated using customizable sampling functions. It's useful for finding global optima 
#' in complex objective functions that may have multiple local minima.
#'
#' @param obj Object of class obsfn, i.e. a function obj(..., fixed, deriv, conditions, env) that returns an objective list, objlist.
#' @param center Parameter values around which the initial values for each fit are randomly sampled.
#'        The initial values are the sum of center and the output of samplefun.
#' @param fits The number of optimization runs to perform. Defaults to 20.
#' @param cores The number of cores to use for parallel processing. Defaults to 1.
#' @param method A character string specifying the optimization method. Defaults to "L-BFGS-B".
#' @param lower A numeric vector of lower bounds for the parameters. Defaults to `-Inf`.
#' @param upper A numeric vector of upper bounds for the parameters. Defaults to `Inf`.
#' @param control A list of control parameters to pass to the optimization algorithm.
#' @param samplefun Function to sample random initial values. It is assumed that samplefun has 
#'        a named parameter "n" which defines how many random numbers are to be returned, 
#'        such as for rnorm or runif. By default rnorm is used.
#' @param seed A seed for the random number generator. If NULL, no seed is set.
#' @param start1stfromCenter Logical, if TRUE the first optimization starts from the center without adding random values.
#' @param ... Additional arguments to pass to the objective function or to samplefun.
#'
#' @return A list where each element contains the result of an individual optimization run.
#'         Each element has the same structure as the return value of the `optimize` function.
#'
#' @import parallel
#' @export
msoptimize <- function(obj, center, fits = 20, cores = 1, 
                      method = "L-BFGS-B", lower = -Inf, upper = Inf, control = list(), 
                      samplefun = stats::rnorm, seed = NULL, start1stfromCenter = TRUE, ...) {
  
  if (!inherits(obj, "objfn")) {
    stop("The object 'obj' must inherit from class 'objfn'.")
  }

  parNames <- attr(obj, 'parameters')                      
  npar <- length(parNames)
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Gather all function arguments
  varargslist <- list(...)
  
  # Determine which arguments should go to samplefun
  namessample <- intersect(names(formals(samplefun)), names(varargslist))
  
  # Create argument list for samplefun
  argssample <- structure(vector("list", length = length(namessample)), names = namessample)
  for (name in namessample) {
    argssample[[name]] <- varargslist[[name]]
  }
  
  # Add n parameter for samplefun
  argssample$n <- npar
  
  # Generate starting points
  start_points <- matrix(0, nrow = fits, ncol = npar)
  colnames(start_points) <- parNames
  
  for (i in 1:fits) {
    if (i == 1 && start1stfromCenter) {
      # First fit starts from center
      start_points[i, ] <- center
    } else {
      # Other fits start from random positions
      start_points[i, ] <- center + do.call(samplefun, argssample)
    }
  }
  
  # Function to perform optimization from a single starting point
  run_optimization <- function(start_point) {
    optimize(
      objfn = obj, 
      parinit = start_point, 
      method = method, 
      lower = lower, 
      upper = upper, 
      control = control, 
      ...
    )
  }
  
  # Run optimizations in parallel if cores > 1
  if (cores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply for parallel processing on Unix-like systems
    results <- parallel::mclapply(
      1:nrow(start_points),
      function(i) {
        result <- run_optimization(start_points[i,])
        result$parinit <- start_points[i,]  # Add initial parameters
        result$index <- i                   # Add index
        return(result)
      },
      mc.cores = min(cores, fits)
    )
  } else {
    # Fall back to sequential processing on Windows or if cores <= 1
    results <- lapply(1:nrow(start_points), function(i) {
      result <- run_optimization(start_points[i,])
      result$parinit <- start_points[i,]  # Add initial parameters
      result$index <- i                   # Add index
      return(result)
    })
  }
  class(results) <- "parlist"
  # Convert results to a list
  return(results)
}


