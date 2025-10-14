#' @export
print.resObjlist <- function(x, n = 10, ...) {
  n <- min(n, length(x$residuals))
  
  cat("Residual object\n")
  cat("===============\n\n")
  cat("Number of residuals:", length(x$residuals), "\n")
  cat("RSS (Residual Sum of Squares):", sum(x$residuals^2), "\n")
  
  # Show individual contributions if available
  all_attrs <- names(attributes(x))
  contrib_attrs <- setdiff(all_attrs, c("names", "class", "residuals", "jacobian", "env", "conditions"))
  
  if (length(contrib_attrs) > 0) {
    cat("\nSpecific Contributions:\n")
    for (a in contrib_attrs) {
      val <- attr(x, a)
      if (is.numeric(val) && length(val) == 1) {
        cat("  ", a, ": ", val, "\n", sep = "")
      }
    }
  }
  
  cat("\nResiduals [1:", n, "]:\n", sep = "")
  print(x$residuals[1:n])
  
  if (!is.null(x$jacobian)) {
    cat("\nJacobian:", nrow(x$jacobian), "×", ncol(x$jacobian), " matrix\n")
    cat("\nParameters:", paste0(head(colnames(x$jacobian), 10), collapse = ", "))
    if (ncol(x$jacobian) > 10) cat(", ...")
    cat("\n")
  } else {
    cat("\nJacobian: NULL\n")
  }
  
  invisible(x)
}


#' @export
print.resObjfn <- function(x, ...) {
  parameters <- attr(x, "parameters")
  conditions <- attr(x, "conditions")
  modelname <- attr(x, "modelname")
  
  cat("Residual-based objective function\n")
  cat("==================================\n\n")
  
  cat("Function signature:\n")
  print(args(x))
  
  if (!is.null(parameters)) {
    cat("\nParameters (", length(parameters), "):\n  ", sep = "")
    cat(paste0(head(parameters, 10), collapse = ", "))
    if (length(parameters) > 10) cat(", ...")
    cat("\n")
  }
  
  if (!is.null(conditions)) {
    cat("\nConditions (", length(conditions), "):\n  ", sep = "")
    cat(paste0(conditions, collapse = ", "), "\n")
  }
  
  if (!is.null(modelname)) {
    cat("\nModel: ", paste0(modelname, collapse = ", "), "\n", sep = "")
  }
  
  invisible(x)
}


#' Residual-based Objective Function for lsqcpp optimization
#'
#' @description
#' Creates an objective function that returns residuals and Jacobians compatible
#' with lsqcpp C++ optimizers (Gauss-Newton, Levenberg-Marquardt, etc.).
#' 
#' Unlike `normL2()` which returns an objlist with value/gradient/hessian,
#' `resL2()` returns residuals suitable for least-squares optimization.
#'
#' @param data object of class \link{datalist}
#' @param x object of class \link{prdfn}
#' @param errmodel object of class \link{obsfn}. If provided, error model parameters
#'   are included in the optimization and artificial residuals for log(sigma) terms
#'   are added.
#' @param times numeric vector, additional time points
#' @param optBLOQ Character: "none", "M1", "M3", or "M4NM" for BLOQ handling
#' @param add_c Numeric: additive constant for error residuals (default: 50).
#'   This constant ensures that the error residuals are positive and properly
#'   scaled for the least-squares optimization.
#' @param attr.name Character: name of the attribute to store -2*log(L) value (default: "data")
#' @param useFitErrorCorrection Logical: apply bias correction for error parameter estimation
#'   (default: TRUE). Applies correction factor sqrt(n/(n-p)) to residuals
#'   to obtain unbiased variance estimates.
#'
#' @return Object of class `resObjfn` - a function that returns a list with:
#'   \item{residuals}{Numeric vector of weighted residuals (plus artificial sigma residuals if errmodel given)}
#'   \item{jacobian}{Matrix with derivatives of residuals w.r.t. all parameters}
#'
#' @details
#' The returned function can be used with lsqcpp optimizers or combined with
#' constraints using the `+` operator.
#' 
#' **Residuals for data points:**
#' For each data point, the weighted residual is computed as:
#' \deqn{r_i = \frac{y_i - f_i}{\sigma_i} \sqrt{c_{corr}}}
#' where \eqn{c_{corr}} is the bias correction factor (if enabled).
#' 
#' **Artificial sigma residuals:**
#' If an error model is provided, artificial residuals are added for each ALOQ 
#' (above limit of quantification) data point, following
#' \deqn{r_{\sigma,i} = \sqrt{2 \log(\sigma_i) + c}}
#' where \eqn{c} is the additive constant `add_c` (default: 50).
#' 
#' **Bias correction:**
#' When `useFitErrorCorrection = TRUE` and error parameters are estimated, a correction
#' factor is applied to ensure unbiased variance estimates:
#' \deqn{c_{corr} = \frac{n_{data}}{n_{data} - n_{params}}}
#' where \eqn{n_{params}} is the number of fitted non-error parameters.
#' This is the Bessel correction factor that accounts for degrees of freedom lost
#' due to parameter estimation.
#' 
#' **Jacobian matrix:**
#' The Jacobian matrix contains the derivatives of residuals with respect to parameters,
#' including the correction factor when applicable. The positive sign ensures correct
#' gradient direction for minimizing the negative log-likelihood via \eqn{\nabla f = J^T r}.
#' 
#' **Likelihood calculation:**
#' The attribute specified by `attr.name` contains:
#' \deqn{-2 \log(L) = \chi^2_{data} + 2 \sum \log(\sigma_i) + n \log(2\pi)}
#' where \eqn{\chi^2_{data}} is the sum of squared weighted residuals (without correction factor).
#'
#' @export
#' @family dMod interface
#' @examples
#' \dontrun{
#' obj <- resL2(data, model, errmodel, optBLOQ = "M3", add_c = 50)
#' result <- obj(pars = initial_pars)
#' # Access -2*log(L)
#' neg2logL <- attr(result, "data")
#' # Use with lsqcpp optimizer
#' }
resL2 <- function(data, x, errmodel = NULL, times = NULL, optBLOQ = "none", add_c = 50, 
                  attr.name = "data", useFitErrorCorrection = TRUE) {
  
  timesD <- sort(unique(c(0, unlist(lapply(data, function(d) d$time)))))
  if (!is.null(times)) timesD <- sort(union(times, timesD))
  
  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  if (!all(data.conditions %in% x.conditions)) 
    stop("Prediction function does not provide predictions for all conditions in data.")
  
  e.conditions <- names(attr(errmodel, "mappings"))
  force(errmodel)
  has_errmodel <- !is.null(errmodel)
  
  err_params <- character(0)
  if(has_errmodel) {
    err_params <- intersect(getParameters(errmodel), getSymbols(getEquations(errmodel)))
  }
  
  controls <- list(times = timesD, conditions = intersect(x.conditions, data.conditions), 
                   optBLOQ = optBLOQ, add_c = add_c, attr.name = attr.name,
                   useFitErrorCorrection = useFitErrorCorrection,
                   err_params = err_params)
  
  # Thread-safe warning flag stored in closure
  correction_warning_issued <- FALSE
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions) {
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    prediction <- x(times = controls$times, pars = pouter, fixed = fixed, 
                    deriv = deriv, conditions = conditions)
    
    # Compute residuals for each condition
    res.list <- lapply(conditions, function(cn) {
      err <- if ((has_errmodel & is.null(e.conditions)) | 
                 (!is.null(e.conditions) && cn %in% e.conditions)) {
        errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
      } else NULL
      resCpp(data[[cn]], prediction[[cn]], err[[cn]], optBLOQ = controls$optBLOQ)
    })
    
    # Extract weighted residuals
    residuals <- unlist(lapply(res.list, `[[`, "weighted.residual"), use.names = FALSE)
    n_data <- length(residuals)
    
    # Calculate correction factor
    fiterrors_correction_factor <- 1
    if (has_errmodel && controls$useFitErrorCorrection) {
      # Count number of fitted non-error parameters (like D2D: qFit==1 & qError!=1)
      all_pars <- names(pouter)
      n_params_fitted <- length(setdiff(all_pars, controls$err_params))
      
      # Bessel correction: n / (n - p)
      if (n_data > n_params_fitted) {
        fiterrors_correction_factor <- n_data / (n_data - n_params_fitted)
      } else {
        fiterrors_correction_factor <- 1
        # Thread-safe warning using closure variable
        if (!correction_warning_issued) {
          correction_warning_issued <<- TRUE
          warning("Not enough data for bias correction (n_data=", n_data, 
                  " <= n_params=", n_params_fitted, 
                  "). Setting correction factor to 1.", 
                  call. = FALSE)
        }
      }
    }
    
    # Apply correction factor to data residuals
    residuals <- residuals * sqrt(fiterrors_correction_factor)
    
    # Add artificial sigma residuals if error model is estimated
    if (has_errmodel) {
      sigma_list <- vector("list", length(res.list))
      for (i in seq_along(res.list)) {
        r <- res.list[[i]]
        reserr <- 2 * log(r$sigma) + controls$add_c
        sigma_list[[i]] <- ifelse(r$bloq | reserr <= 0, 0, sqrt(reserr))
      }
      sigma_res <- unlist(sigma_list, use.names = FALSE)
      residuals <- c(residuals, sigma_res)
    }
    
    jacobian <- NULL
    if (deriv) {
      # Extract and process Jacobians
      jac_data <- vector("list", length(res.list))
      
      for (i in seq_along(res.list)) {
        r <- res.list[[i]]
        dxdp <- attr(r, "deriv")
        dsdp <- attr(r, "deriv.err")
        
        if (!is.null(dxdp)) dxdp <- as.matrix(dxdp[, -(1:2), drop = FALSE])
        if (!is.null(dsdp)) dsdp <- as.matrix(dsdp[, -(1:2), drop = FALSE])
        
        s <- r$sigma
        wr <- r$weighted.residual
        
        if (is.null(dsdp) || ncol(dsdp) == 0) {
          # Apply correction factor to jacobian
          jac_data[[i]] <- sweep(dxdp, 1, sqrt(fiterrors_correction_factor)/s, "*")
        } else {
          all_cols <- union(colnames(dxdp), colnames(dsdp))
          n_cols <- length(all_cols)
          
          if (ncol(dxdp) < n_cols) {
            tmp <- matrix(0, nrow(dxdp), n_cols, dimnames = list(NULL, all_cols))
            tmp[, colnames(dxdp)] <- dxdp
            dxdp <- tmp
          }
          if (ncol(dsdp) < n_cols) {
            tmp <- matrix(0, nrow(dsdp), n_cols, dimnames = list(NULL, all_cols))
            tmp[, colnames(dsdp)] <- dsdp
            dsdp <- tmp
          }
          
          # Apply correction factor to jacobian
          jac_data[[i]] <- sweep(dxdp, 1, sqrt(fiterrors_correction_factor)/s, "*") + 
            sweep(dsdp, 1, wr * sqrt(fiterrors_correction_factor)/s, "*")
        }
      }
      
      # Sigma residuals jacobian
      if (has_errmodel) {
        jac_sigma <- vector("list", length(res.list))
        
        for (i in seq_along(res.list)) {
          r <- res.list[[i]]
          dsdp <- attr(r, "deriv.err")
          
          if (is.null(dsdp)) {
            jac_sigma[[i]] <- matrix(0, nrow(r), 0)
          } else {
            dsdp <- as.matrix(dsdp[, -(1:2), drop = FALSE])
            if (ncol(dsdp) == 0) {
              jac_sigma[[i]] <- dsdp
            } else {
              reserr <- 2 * log(r$sigma) + controls$add_c
              scale <- ifelse(r$bloq | reserr <= 0, 0, 1 / (r$sigma * sqrt(reserr)))
              jac_sigma[[i]] <- sweep(dsdp, 1, scale, "*")
            }
          }
        }
        
        all_pars <- unique(c(unlist(lapply(jac_data, colnames)), 
                             unlist(lapply(jac_sigma, colnames))))
        n_pars <- length(all_pars)
        
        expand <- function(jac) {
          if (nrow(jac) == 0) return(matrix(0, 0, n_pars, dimnames = list(NULL, all_pars)))
          if (ncol(jac) == n_pars && identical(colnames(jac), all_pars)) return(jac)
          out <- matrix(0, nrow(jac), n_pars, dimnames = list(NULL, all_pars))
          out[, colnames(jac)] <- jac
          out
        }
        
        jacobian <- rbind(do.call(rbind, lapply(jac_data, expand)),
                          do.call(rbind, lapply(jac_sigma, expand)))
      } else {
        all_pars <- unique(unlist(lapply(jac_data, colnames)))
        n_pars <- length(all_pars)
        
        expand <- function(jac) {
          if (nrow(jac) == 0) return(matrix(0, 0, n_pars, dimnames = list(NULL, all_pars)))
          if (ncol(jac) == n_pars && identical(colnames(jac), all_pars)) return(jac)
          out <- matrix(0, nrow(jac), n_pars, dimnames = list(NULL, all_pars))
          out[, colnames(jac)] <- jac
          out
        }
        jacobian <- do.call(rbind, lapply(jac_data, expand))
      }
    }
    
    # Calculate chi-square WITHOUT correction factor (for likelihood calculation)
    chisquare <- sum(residuals[seq_len(n_data)]^2) / fiterrors_correction_factor
    
    # Calculate -2*log(L) = chi2_data + 2*Σlog(σ) + n*log(2π)
    neg2logL <- chisquare + 2 * sum(vapply(res.list, function(r) sum(log(r$sigma)), numeric(1))) + n_data * log(2 * pi)
    
    out <- resObjlist(residuals = residuals, jacobian = jacobian)
    attr(out, "conditions") <- conditions
    attr(out, "chisquare") <- chisquare
    attr(out, controls$attr.name) <- neg2logL
    out
  }
  
  class(myfn) <- c("resObjfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(x, "parameters")
  attr(myfn, "modelname") <- modelname(x, errmodel)
  myfn
}


#' Soft L2 constraint on parameters (residual form)
#'
#' @description
#' Creates a constraint function that returns residuals for use with lsqcpp.
#' Equivalent to `constraintL2()` but returns residuals instead of objlist.
#'
#' @param mu Named numeric, the prior values
#' @param sigma Named numeric or character, the prior standard deviations.
#'   If character, sigma parameters are estimated (on log scale) and artificial
#'   residuals are added for the log-likelihood term.
#' @param add_c Numeric: additive constant for sigma residuals (default: 50).
#'   Only used when sigma is estimated (character).
#' @param attr.name character. The constraint value is additionally returned in an 
#'   attribute with this name
#' @param condition Character, condition for which constraint applies
#'
#' @return Object of class `resObjfn`
#'
#' @details
#' Returns residuals `(p - mu) / sigma` for each parameter.
#' 
#' **Fixed sigma (numeric):**
#' Only parameter residuals are computed: `r = (p - mu) / σ`
#' 
#' **Estimated sigma (character):**
#' Additional artificial residuals are added for each sigma parameter, following
#' \deqn{r_{\sigma} = \sqrt{2 \log(\sigma) + c}}
#' where \eqn{c} is the additive constant `add_c` (default: 50).
#' 
#' **Jacobian matrix:**
#' For parameter residuals with respect to parameters:
#' \deqn{\frac{\partial r}{\partial p} = \frac{1}{\sigma}}
#' 
#' For parameter residuals with respect to sigma (on log scale):
#' \deqn{\frac{\partial r}{\partial \log \sigma} = -\frac{p - \mu}{\sigma}}
#' 
#' For sigma residuals with respect to sigma (on log scale):
#' \deqn{\frac{\partial r_{\sigma}}{\partial \log \sigma} = \frac{1}{\sqrt{2 \log(\sigma) + c}}}
#' 
#' **Constraint value:**
#' The constraint value (stored in the attribute) is computed as:
#' \deqn{constraint = \sum r_{param}^2 + 2 \sum \log(\sigma) + c \cdot n_{sigma}}
#' 
#' This ensures proper likelihood evaluation where:
#' - `Σr_{param}²` contributes the data fit term
#' - `2*Σlog(σ)` contributes the log-likelihood term for variances
#' - The additive constant ensures numerical stability
#'
#' @export
#' @family dMod interface
#' @examples
#' \dontrun{
#' # Fixed sigma
#' constraint1 <- constraintResL2(c(k1 = 1, k2 = 2), sigma = 0.5)
#' 
#' # Estimated sigma
#' constraint2 <- constraintResL2(c(k1 = 1, k2 = 2), sigma = "sigma_prior", add_c = 50)
#' obj <- resL2(data, model) + constraint2
#' }
constraintResL2 <- function(mu, sigma = 1, add_c = 50, attr.name = "prior", condition = NULL) {
  
  estimateSigma <- is.character(sigma)
  
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  if (is.null(names(sigma)))
    names(sigma) <- names(mu)
  if (!all(names(mu) %in% names(sigma)))
    stop("Names of sigma and names of mu do not match.")
  
  sigma <- sigma[names(mu)]
  controls <- list(mu = mu, sigma = sigma, estimateSigma = estimateSigma, 
                   add_c = add_c, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition, env = NULL) {
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    mu <- controls$mu
    sigma <- controls$sigma
    estimateSigma <- controls$estimateSigma
    add_c <- controls$add_c
    attr.name <- controls$attr.name
    
    # Handle list of parameters (multiple conditions)
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      if (length(available) == 0 || (!is.null(condition) && !condition %in% conditions)) 
        return(NULL)
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    # Pre-allocate for efficiency
    n_cond <- length(pouter)
    res_list <- vector("list", n_cond)
    jac_list <- if (deriv) vector("list", n_cond) else NULL
    chi2_param <- 0
    sum_log_sigma_total <- 0
    n_sigma_total <- 0L
    
    for (idx in seq_len(n_cond)) {
      p <- pouter[[idx]]
      
      pars <- c(p, fixed)[names(mu)]
      p1 <- setdiff(intersect(names(mu), names(p)), names(fixed))
      
      # Get sigma values
      if (estimateSigma) {
        sigmapars <- sigma
        sigma_vals <- exp(c(p, fixed)[sigma])
        names(sigma_vals) <- names(mu)
        p2 <- setdiff(intersect(unique(sigmapars), names(p)), names(fixed))
      } else {
        sigma_vals <- sigma
      }
      
      # Compute parameter residuals: (p - mu) / sigma
      residuals <- ((pars - mu) / sigma_vals)[p1]
      n_param_res <- length(residuals)
      
      # Calculate chi2 for parameters
      chi2_param <- chi2_param + sum(residuals^2)
      
      # Add artificial residuals for log(sigma) if estimated
      if (estimateSigma && length(p2) > 0) {
        log_sigma <- log(sigma_vals[p1])
        sum_log_sigma_total <- sum_log_sigma_total + sum(log_sigma)
        
        reserr <- 2 * log_sigma + add_c
        sigma_residuals <- ifelse(reserr <= 0, 0, sqrt(reserr))
        n_sigma_total <- n_sigma_total + sum(sigma_residuals > 0)
        
        residuals <- c(residuals, structure(sigma_residuals, names = paste0("log_", names(residuals))))
      }
      
      res_list[[idx]] <- residuals
      
      if (deriv) {
        n_sigma_res <- length(residuals) - n_param_res
        n_total <- n_param_res + n_sigma_res
        
        jacobian <- matrix(0, nrow = n_total, ncol = length(p),
                           dimnames = list(names(residuals), names(p)))
        
        # Jacobian for parameter residuals w.r.t. parameters
        diag(jacobian[seq_len(n_param_res), p1]) <- 1 / sigma_vals[p1]
        
        if (estimateSigma && n_sigma_res > 0) {
          # Jacobian for parameter residuals w.r.t. sigma parameters (on log scale)
          for (i in seq_len(n_param_res)) {
            par_name <- p1[i]
            sigma_name <- sigmapars[par_name]
            if (sigma_name %in% p2) {
              jacobian[i, sigma_name] <- -(pars[par_name] - mu[par_name]) / sigma_vals[par_name]
            }
          }
          
          # Jacobian for sigma residuals w.r.t. sigma parameters (on log scale)
          log_sigma <- log(sigma_vals[p1])
          reserr <- 2 * log_sigma + add_c
          for (i in seq_len(n_param_res)) {
            par_name <- p1[i]
            sigma_name <- sigmapars[par_name]
            if (sigma_name %in% p2 && reserr[i] > 0) {
              jacobian[n_param_res + i, sigma_name] <- 1 / sqrt(reserr[i])
            }
          }
        }
        
        # Apply chain rule if parameter transformation exists
        dP <- attr(p, "deriv")
        if (!is.null(dP)) jacobian <- jacobian %*% dP
        
        jac_list[[idx]] <- jacobian
      }
    }
    
    # Combine results
    out <- resObjlist(
      residuals = unlist(res_list, use.names = FALSE),
      jacobian = if (deriv) do.call(rbind, jac_list) else NULL
    )
    
    # Calculate constraint value: chi2_param + 2*sum(log(sigma)) + add_c * n_sigma
    constraint_value <- chi2_param + 2 * sum_log_sigma_total + add_c * n_sigma_total
    
    attr(out, controls$attr.name) <- constraint_value
    attr(out, "env") <- env
    
    return(out)
  }
  
  class(myfn) <- c("resObjfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  return(myfn)
}


#' L2 objective function for validation data point (residual form)
#'
#' @description
#' Creates residual for a single validation point.
#' Equivalent to `datapointL2()` but returns residual instead of objlist.
#'
#' @param name Character, prediction name (e.g., state name)
#' @param time Numeric, time point
#' @param value Character, parameter name containing the validation value
#' @param sigma Numeric, uncertainty
#' @param condition Character, condition name
#'
#' @return Object of class `resObjfn`
#' 
#' @details
#' Returns a single residual `(prediction - value) / sigma` for the validation point.
#' The prediction is extracted from the environment passed via `env` argument.
#'
#' @export
#' @family dMod interface
#' @examples
#' \dontrun{
#' vali <- datapointResL2("Drug", time = 24, value = "Drug_24h", 
#'                        sigma = 0.1, condition = "treatment")
#' obj <- resL2(data, model) + vali
#' }
datapointResL2 <- function(name, time, value, sigma = 1, condition) {
  
  controls <- list(mu = structure(name, names = value)[1], time = time[1], sigma = sigma[1])
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {
    
    mu <- controls$mu
    time <- controls$time
    sigma <- controls$sigma
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    
    if (!is.null(env)) {
      prediction <- as.list(env)$prediction 
    } else {
      stop("No prediction available. Use env argument to pass environment with prediction.")
    }
    
    # Check condition overlap
    if (!is.null(conditions) && !condition %in% conditions) return(NULL)
    if (is.null(conditions) && !condition %in% names(prediction))
      stop("datapointResL2 requests unavailable condition.")
    
    # Extract prediction at time point
    time.index <- which(prediction[[condition]][,"time"] == time)
    if (length(time.index) == 0) 
      stop("datapointResL2 requests time point without prediction. Add time in resL2().")
    
    withDeriv <- !is.null(attr(prediction[[condition]], "deriv"))
    pred <- prediction[[condition]][time.index, mu]
    
    # Compute residual: (pred - value) / sigma
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    
    residual <- (pred - c(fixed, pouter)[names(mu)]) / sigma
    names(residual) <- paste0("datapoint_", mu)
    
    jacobian <- NULL
    if (deriv && withDeriv) {
      deriv_pred <- attr(prediction[[condition]], "deriv")[time.index, ]
      mu.para <- intersect(paste(mu, parapar, sep = "."), names(deriv_pred))
      
      jacobian <- matrix(0, nrow = 1, ncol = length(pouter),
                         dimnames = list(names(residual), names(pouter)))
      
      if (length(parapar) > 0) {
        jacobian[1, parapar] <- as.numeric(deriv_pred[mu.para]) / sigma
      }
      if (length(datapar) > 0) {
        jacobian[1, datapar] <- -1 / sigma
      }
    }
    
    out <- resObjlist(residuals = residual, jacobian = jacobian)
    attr(out, "prediction") <- pred
    return(out)
  }
  
  class(myfn) <- c("resObjfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  return(myfn)
}