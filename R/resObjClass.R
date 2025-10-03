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
#' @param errmodel object of class \link{obsfn}
#' @param times numeric vector, additional time points
#' @param optBLOQ Character: "none", "M1", "M3", or "M4NM" for BLOQ handling
#'
#' @return Object of class `resObjfn` - a function that returns a list with:
#'   \item{residuals}{Numeric vector of weighted residuals plus artificial residuals for log(sigma)}
#'   \item{jacobian}{Matrix with derivatives of residuals w.r.t. all parameters (model + error)}
#'
#' @details
#' The returned function can be used with lsqcpp optimizers or combined with
#' constraints using `+` operator:
#' 
#' ```
#' obj <- resL2(data, model) + constraintResL2(mu, sigma)
#' ```
#' 
#' **Important:** This function includes artificial residuals for the log(sigma) terms
#' from the likelihood. For each ALOQ data point, an additional residual 
#' `sqrt(log(2*pi*sigma^2))` is added to ensure the likelihood is correctly 
#' represented as a sum of squared residuals.
#'
#' @export
#' @family dMod interface
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
    
    # Add artificial residuals for log(sigma) terms (only for ALOQ)
    sigma_residuals <- unlist(lapply(res.list, function(r) {
      ifelse(r$bloq, 0, sqrt(log(2 * pi * r$sigma^2)))
    }))
    residuals <- c(residuals, sigma_residuals)
    
    jacobian <- NULL
    if (deriv) {
      # Helper to extract and combine jacobians
      get_jac <- function(r, attr_name, scale_by_sigma = FALSE) {
        d <- attr(r, attr_name)
        if (is.null(d)) return(matrix(0, nrow = nrow(r), ncol = 0))
        jac <- as.matrix(d[, -(1:2), drop = FALSE])
        if (scale_by_sigma) jac <- sweep(jac, 1, r$sigma, "/")
        jac
      }
      
      # Data residuals: combine model + error derivatives
      jac_data <- lapply(res.list, function(r) {
        cbind(get_jac(r, "deriv"), get_jac(r, "deriv.err"))
      })
      
      # Sigma residuals: only error derivatives, scaled by 1/sigma, set to 0 for BLOQ
      jac_sigma <- lapply(res.list, function(r) {
        jac <- cbind(matrix(0, nrow = nrow(r), ncol = ncol(get_jac(r, "deriv"))),
                     get_jac(r, "deriv.err", scale_by_sigma = TRUE))
        jac[r$bloq, ] <- 0
        jac
      })
      
      # Get all unique parameter names and expand
      all_pars <- unique(unlist(lapply(c(jac_data, jac_sigma), colnames)))
      expand <- function(jac) {
        out <- matrix(0, nrow(jac), length(all_pars), dimnames = list(NULL, all_pars))
        out[, colnames(jac)] <- jac
        out
      }
      
      jacobian <- rbind(
        do.call(rbind, lapply(jac_data, expand)),
        do.call(rbind, lapply(jac_sigma, expand))
      )
    }
    
    out <- resObj(residuals = residuals, jacobian = jacobian)
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
#' @param sigma Named numeric or character, the prior standard deviations
#' @param condition Character, condition for which constraint applies
#'
#' @return Object of class `resObjfn`
#'
#' @details
#' Returns residuals `(p - mu) / sigma` for each parameter.
#' When squared and summed, these give the L2 constraint value.
#'
#' @export
#' @family dMod interface
constraintResL2 <- function(mu, sigma = 1, condition = NULL) {
  
  estimateSigma <- is.character(sigma)
  
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  if (is.null(names(sigma)))
    names(sigma) <- names(mu)
  if (!all(names(mu) %in% names(sigma)))
    stop("Names of sigma and names of mu do not match.")
  
  sigma <- sigma[names(mu)]
  controls <- list(mu = mu, sigma = sigma)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = condition) {
    
    pouter <- list(...)[match.fnargs(list(...), "pars")][[1]]
    mu <- controls$mu
    sigma <- controls$sigma
    
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
      
      # Estimate sigma from parameters if character
      if (estimateSigma) {
        sigmapars <- sigma
        sigma <- exp(c(p, fixed)[sigma])
        names(sigma) <- names(mu)
        p2 <- setdiff(intersect(unique(sigmapars), names(p)), names(fixed))
      }
      
      # Compute residuals: (p - mu) / sigma
      residuals <- ((pars - mu) / sigma)[p1]
      
      jacobian <- NULL
      if (deriv) {
        jacobian <- matrix(0, nrow = length(residuals), ncol = length(p),
                           dimnames = list(names(residuals), names(p)))
        
        # dr/dp for parameter residuals
        diag(jacobian[, p1]) <- 1 / sigma[p1]
        
        # dr/dsigma if estimating sigma
        if (estimateSigma && length(p2) > 0) {
          for (i in seq_along(residuals)) {
            par_name <- names(residuals)[i]
            sigma_name <- sigmapars[par_name]
            if (sigma_name %in% p2)
              jacobian[i, sigma_name] <- -residuals[i]
          }
        }
        
        # Apply chain rule if parameter transformation exists
        dP <- attr(p, "deriv")
        if (!is.null(dP)) jacobian <- jacobian %*% dP
      }
      
      list(residuals = residuals, jacobian = jacobian)
    })
    
    # Combine results
    out <- resObj(
      residuals = unlist(lapply(outlist, `[[`, "residuals")),
      jacobian = if (deriv) do.call(rbind, lapply(outlist, `[[`, "jacobian")) else NULL
    )
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
#' @export
#' @family dMod interface
datapointResL2 <- function(name, time, value, sigma = 1, condition) {
  
  controls <- list(
    mu = structure(name, names = value)[1],
    time = time[1],
    sigma = sigma[1]
  )
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, conditions = NULL, env = NULL) {
    
    mu <- controls$mu
    time <- controls$time
    sigma <- controls$sigma
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    if (!is.null(env)) {
      prediction <- as.list(env)$prediction 
    } else {
      stop("No prediction available. Use env argument to pass environment with prediction.")
    }
    
    # Check condition overlap
    if (!is.null(conditions) && !condition %in% conditions) 
      return(NULL)
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
    
    out <- reslist(residuals = residual, jacobian = jacobian)
    attr(out, "prediction") <- pred
    
    return(out)
  }
  
  class(myfn) <- c("resObjfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  
  return(myfn)
}