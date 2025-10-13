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
#'
#' @return Object of class `resObjfn` - a function that returns a list with:
#'   \item{residuals}{Numeric vector of weighted residuals (plus artificial sigma residuals if errmodel given)}
#'   \item{jacobian}{Matrix with derivatives of residuals w.r.t. all parameters}
#'
#' @details
#' The returned function can be used with lsqcpp optimizers or combined with
#' constraints using the `+` operator.
#' 
#' **Artificial sigma residuals:** If an error model is provided, artificial residuals
#' `sign(log σ) * sqrt(2*|log σ|)` are added for each ALOQ data point to represent
#' the log-likelihood term. When squared, these give `2*log(σ)`, ensuring that
#' `½·Σ[r²] = NLL`.
#'
#' @export
#' @family dMod interface
#' @examples
#' \dontrun{
#' obj <- resL2(data, model, errmodel, optBLOQ = "M3")
#' result <- obj(pars = initial_pars)
#' # Use with lsqcpp optimizer
#' }
resL2 <- function(data, x, errmodel = NULL, times = NULL, optBLOQ = "none") {
  
  timesD <- sort(unique(c(0, unlist(lapply(data, function(d) d$time)))))
  if (!is.null(times)) timesD <- sort(union(times, timesD))
  
  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  if (!all(data.conditions %in% x.conditions)) 
    stop("Prediction function does not provide predictions for all conditions in data.")
  
  e.conditions <- names(attr(errmodel, "mappings"))
  controls <- list(times = timesD, conditions = intersect(x.conditions, data.conditions), optBLOQ = optBLOQ)
  force(errmodel)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = controls$conditions) {
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    prediction <- x(times = controls$times, pars = pouter, fixed = fixed, 
                    deriv = deriv, conditions = conditions)
    
    # Compute residuals for each condition
    res.list <- lapply(conditions, function(cn) {
      err <- if ((!is.null(errmodel) & is.null(e.conditions)) | 
                 (!is.null(e.conditions) && cn %in% e.conditions)) {
        errmodel(out = prediction[[cn]], pars = getParameters(prediction[[cn]]), conditions = cn)
      } else NULL
      resCpp(data[[cn]], prediction[[cn]], err[[cn]], optBLOQ = controls$optBLOQ)
    })
    
    # Combine weighted residuals
    residuals <- unlist(lapply(res.list, function(r) r$weighted.residual))
    
    # Add artificial residuals for log(sigma) if error model is estimated
    # r = sign(log σ) * sqrt(2*|log σ|), then r² = 2*log σ
    if (!is.null(errmodel)) {
      sigma_residuals <- unlist(lapply(res.list, function(r) {
        log_sigma <- log(r$sigma)
        ifelse(r$bloq, 0, sign(log_sigma) * sqrt(2 * abs(log_sigma)))
      }))
      residuals <- c(residuals, sigma_residuals)
      has_sigma_residuals <- TRUE
    } else {
      has_sigma_residuals <- FALSE
    }
    
    jacobian <- NULL
    if (deriv) {
      # Helper to extract jacobians
      get_jac <- function(r, attr_name) {
        d <- attr(r, attr_name)
        if (is.null(d)) return(matrix(0, nrow = nrow(r), ncol = 0))
        as.matrix(d[, -(1:2), drop = FALSE])
      }
      
      # Data residuals: compute dwrdp from dxdp and dsdp
      jac_data <- lapply(res.list, function(r) {
        dxdp <- get_jac(r, "deriv")      # ∂f/∂p
        dsdp <- get_jac(r, "deriv.err")  # ∂σ/∂p
        wr <- r$weighted.residual
        s <- r$sigma
        
        # Compute dwrdp = ∂((y-f)/σ)/∂p = -1/σ · ∂f/∂p - wr/σ · ∂σ/∂p
        # Handle case where dsdp might be empty
        if (ncol(dsdp) == 0) {
          # No error parameters: dwrdp = -1/σ · ∂f/∂p
          dwrdp <- sweep(dxdp, 1, -1/s, "*")
        } else {
          # Expand to common parameter set if needed
          all_cols <- union(colnames(dxdp), colnames(dsdp))
          if (ncol(dxdp) != length(all_cols)) {
            dxdp_new <- matrix(0, nrow = nrow(dxdp), ncol = length(all_cols),
                               dimnames = list(NULL, all_cols))
            dxdp_new[, colnames(dxdp)] <- dxdp
            dxdp <- dxdp_new
          }
          if (ncol(dsdp) != length(all_cols)) {
            dsdp_new <- matrix(0, nrow = nrow(dsdp), ncol = length(all_cols),
                               dimnames = list(NULL, all_cols))
            dsdp_new[, colnames(dsdp)] <- dsdp
            dsdp <- dsdp_new
          }
          # Now compute dwrdp with both terms
          dwrdp <- sweep(dxdp, 1, -1/s, "*") - sweep(dsdp, 1, wr/s, "*")
        }
        
        return(dwrdp)
      })
      
      # Sigma residuals jacobian if error model present
      if (has_sigma_residuals) {
        jac_sigma <- lapply(res.list, function(r) {
          dsdp <- get_jac(r, "deriv.err")
          if (ncol(dsdp) == 0) return(matrix(0, nrow = nrow(r), ncol = 0))
          
          log_sigma <- log(r$sigma)
          s <- r$sigma
          
          # dlogsdp = 1/σ · ∂σ/∂p
          dlogsdp <- sweep(dsdp, 1, 1/s, "*")
          
          # ∂r_sigma/∂p where r_sigma = sign(log σ) * sqrt(2*|log σ|)
          # Scale factor: 1/sqrt(|log σ|/2)
          scale_factor <- ifelse(r$bloq | abs(log_sigma) < 1e-10, 0,
                                 1 / sqrt(abs(log_sigma) / 2))
          
          jac_err <- sweep(dlogsdp, 1, scale_factor, "*")
          
          # Expand to match columns of jac_data
          return(jac_err)
        })
      } else {
        jac_sigma <- list()
      }
      
      # Expand all jacobians to same parameter set
      all_pars <- unique(unlist(lapply(c(jac_data, jac_sigma), colnames)))
      expand <- function(jac) {
        if (nrow(jac) == 0) return(matrix(0, 0, length(all_pars), dimnames = list(NULL, all_pars)))
        out <- matrix(0, nrow(jac), length(all_pars), dimnames = list(NULL, all_pars))
        out[, colnames(jac)] <- jac
        out
      }
      
      jacobian <- if (has_sigma_residuals) {
        rbind(do.call(rbind, lapply(jac_data, expand)),
              do.call(rbind, lapply(jac_sigma, expand)))
      } else {
        do.call(rbind, lapply(jac_data, expand))
      }
    }
    
    out <- resObjlist(residuals = residuals, jacobian = jacobian)
    attr(out, "conditions") <- conditions
    return(out)
  }
  
  class(myfn) <- c("resObjfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(x, "parameters")
  attr(myfn, "modelname") <- modelname(x, errmodel)
  return(myfn)
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
#' @param attr.name character. The constraint value is additionally returned in an 
#'   attributed with this name
#' @param condition Character, condition for which constraint applies
#'
#' @return Object of class `resObjfn`
#'
#' @details
#' Returns residuals `(p - mu) / sigma` for each parameter.
#' 
#' If `sigma` is character (estimated parameters), additional artificial residuals
#' `sign(log σ) * sqrt(2*|log σ|)` are added to represent the log(sigma^2) term.
#' When squared, these give `2*log(σ)`, ensuring that `Σ[r²]` equals the 
#' constraint value from `constraintL2()`.
#'
#' @export
#' @family dMod interface
#' @examples
#' \dontrun{
#' # Fixed sigma
#' constraint1 <- constraintResL2(c(k1 = 1, k2 = 2), sigma = 0.5)
#' 
#' # Estimated sigma
#' constraint2 <- constraintResL2(c(k1 = 1, k2 = 2), sigma = "sigma_prior")
#' obj <- resL2(data, model) + constraint2
#' }
constraintResL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {
  
  estimateSigma <- is.character(sigma)
  
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  if (is.null(names(sigma)))
    names(sigma) <- names(mu)
  if (!all(names(mu) %in% names(sigma)))
    stop("Names of sigma and names of mu do not match.")
  
  sigma <- sigma[names(mu)]
  controls <- list(mu = mu, sigma = sigma, estimateSigma = estimateSigma, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition, env = NULL) {
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    mu <- controls$mu
    sigma <- controls$sigma
    estimateSigma <- controls$estimateSigma
    attr.name <- controls$attr.name
    
    # Handle list of parameters (multiple conditions)
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      if (length(available) == 0 || (!is.null(condition) && !condition %in% conditions)) 
        return(NULL)
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      pars <- c(p, fixed)[names(mu)]
      p1 <- setdiff(intersect(names(mu), names(p)), names(fixed))
      
      # Get sigma values
      if (estimateSigma) {
        sigmapars <- sigma
        sigma <- exp(c(p, fixed)[sigma])
        names(sigma) <- names(mu)
        p2 <- setdiff(intersect(unique(sigmapars), names(p)), names(fixed))
      }
      
      # Compute parameter residuals: (p - mu) / sigma
      residuals <- ((pars - mu) / sigma)[p1]
      
      # Add artificial residuals for log(sigma^2) = 2*log(sigma) if estimated
      # r = sign(log σ) * sqrt(2*|log σ|), then r² = 2*log σ
      if (estimateSigma && length(p2) > 0) {
        log_sigma <- log(sigma[p1])
        sigma_residuals <- sign(log_sigma) * sqrt(2 * abs(log_sigma))
        residuals <- c(residuals, structure(sigma_residuals, names = paste0("log_", names(residuals))))
      }
      
      jacobian <- NULL
      if (deriv) {
        n_param_res <- length(p1)
        n_sigma_res <- if (estimateSigma && length(p2) > 0) length(p1) else 0
        n_total <- n_param_res + n_sigma_res
        
        jacobian <- matrix(0, nrow = n_total, ncol = length(p),
                           dimnames = list(names(residuals), names(p)))
        
        # Jacobian for parameter residuals w.r.t. parameters
        diag(jacobian[1:n_param_res, p1]) <- 1 / sigma[p1]
        
        if (estimateSigma && n_sigma_res > 0) {
          # Jacobian for parameter residuals w.r.t. sigma parameters (on log scale)
          # ∂((p-mu)/σ)/∂(log σ) = -(p-mu)/σ
          for (i in seq_along(p1)) {
            par_name <- p1[i]
            sigma_name <- sigmapars[par_name]
            if (sigma_name %in% p2) {
              jacobian[i, sigma_name] <- -(pars[par_name] - mu[par_name]) / sigma[par_name]
            }
          }
          
          # Jacobian for sigma residuals w.r.t. sigma parameters (on log scale)
          # r = sign(log σ) * sqrt(2*|log σ|)
          # ∂r/∂(log σ) = 1 / sqrt(|log σ|/2)
          log_sigma <- log(sigma[p1])
          for (i in seq_along(p1)) {
            par_name <- p1[i]
            sigma_name <- sigmapars[par_name]
            if (sigma_name %in% p2) {
              jacobian[n_param_res + i, sigma_name] <- 1 / sqrt(abs(log_sigma[i]) / 2)
            }
          }
        }
        
        # Apply chain rule if parameter transformation exists
        dP <- attr(p, "deriv")
        if (!is.null(dP)) jacobian <- jacobian %*% dP
      }
      
      list(residuals = residuals, jacobian = jacobian)
    })
    
    # Combine results
    out <- resObjlist(
      residuals = unlist(lapply(outlist, `[[`, "residuals")),
      jacobian = if (deriv) do.call(rbind, lapply(outlist, `[[`, "jacobian")) else NULL
    )
    
    # Add constraint value attribute: sum of squared residuals
    attr(out, controls$attr.name) <- sum(out$residuals^2)
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