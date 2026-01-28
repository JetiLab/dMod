## Methods for class "objfn" -----------------------------------------------



## Class "objlist" and its constructors ------------------------------------



#' Generate objective list from numeric vector
#' 
#' @param p Named numeric vector
#' @return list with entries value (\code{0}), 
#' gradient (\code{rep(0, length(p))}) and 
#' hessian (\code{matrix(0, length(p), length(p))}) of class \code{obj}.
#' @examples
#' p <- c(A = 1, B = 2)
#' as.objlist(p)
#' @export
as.objlist <- function(p) {
  
  objlist(value = 0,
          gradient = structure(rep(0, length(p)), names = names(p)),
          hessian = matrix(0, length(p), length(p), dimnames = list(names(p), names(p))))
  
}


#' Compute a differentiable box prior
#' 
#' @param p Named numeric, the parameter value
#' @param mu Named numeric, the prior values, means of boxes
#' @param sigma Named numeric, half box width
#' @param k Named numeric, shape of box; if 0 a quadratic prior is obtained, the higher k the more box shape, gradient at border of the box (-sigma, sigma) is equal to sigma*k
#' @param fixed Named numeric with fixed parameter values (contribute to the prior value but not to gradient and Hessian)
#' @return list with entries: value (numeric, the weighted residual sum of squares), 
#' gradient (numeric, gradient) and 
#' hessian (matrix of type numeric). Object of class \code{objlist}.
constraintExp2 <- function(p, mu, sigma = 1, k = 0.05, fixed=NULL) {
  
  
  ##
  ## This function need to be extended according to constraintL2()
  ## The parameters sigma and k need to be replaced by more
  ## meaningful parameters.
  ##
  
  kmin <- 1e-5
  
  ## Augment sigma if length = 1
  if(length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu)) 
  ## Augment k if length = 1
  if(length(k) == 1) 
    k <- structure(rep(k, length(mu)), names = names(mu))
  
  k <- sapply(k, function(ki){
    if(ki < kmin){
      kmin
    } else ki
  })
  
  
  ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
  par.fixed <- intersect(names(mu), names(fixed))
  sumOfFixed <- 0
  if(!is.null(par.fixed)) sumOfFixed <- sum(0.5*(exp(k[par.fixed]*((fixed[par.fixed] - mu[par.fixed])/sigma[par.fixed])^2)-1)/(exp(k[par.fixed])-1))
  
  
  par <- intersect(names(mu), names(p))
  t <- p[par]
  mu <- mu[par]
  s <- sigma[par]
  k <- k[par]
  
  # Compute prior value and derivatives 
  
  gr <- rep(0, length(t)); names(gr) <- names(t)
  hs <- matrix(0, length(t), length(t), dimnames = list(names(t), names(t)))
  
  val <- sum(0.5*(exp(k*((t-mu)/s)^2)-1)/(exp(k)-1)) + sumOfFixed
  gr <- (k*(t-mu)/(s^2)*exp(k*((t-mu)/s)^2)/(exp(k)-1))
  diag(hs)[par] <- k/(s*s)*exp(k*((t-mu)/s)^2)/(exp(k)-1)*(1+2*k*(t-mu)/(s^2))
  
  dP <- attr(p, "deriv")
  if(!is.null(dP)) {
    gr <- as.vector(gr%*%dP); names(gr) <- colnames(dP)
    hs <- t(dP)%*%hs%*%dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
  }
  
  objlist(value=val,gradient=gr,hessian=hs)
  
}


#' L2 norm between data and model prediction
#'
#' @description
#' For parameter estimation and optimization, an objective function is needed.
#' `normL2` returns an objective function based on the (negative log-likelihood)
#' L2 norm between data and model prediction. The resulting objective function
#' can be used for optimization with the trust optimizer, see [mstrust].
#'
#' @param data Object of class [datalist].
#' @param x Object of class [prdfn].
#' @param errmodel Optional object of class [obsfn]. The error model does not
#'   need to be defined for all conditions.
#' @param times Numeric vector of additional time points where the prediction
#'   function is evaluated. If NULL, time points are extracted from the data.
#'   If the prediction function uses events, event times should be provided here.
#' @param attr.name Character string. The objective value is additionally returned
#'   as an attribute with this name.
#' @param use.bessel Logical. If TRUE and an error model is present, applies a
#'   global Bessel correction to variance estimates to account for finite-sample
#'   bias. Default is TRUE if an error model is provided, FALSE otherwise.
#'
#' @return
#' An object of class `objfn`, i.e. a function
#' `obj(pars, fixed, deriv, deriv2, conditions, env, cores)` returning an
#' [objlist].
#'
#' @details
#' Objective functions can be combined using the `+` operator, see [sumobjfn].
#'
#' The Bessel correction is applied globally across all conditions and is given by
#' \deqn{\sqrt{n / (n - p)}}
#' where \eqn{n} is the total number of data points and \eqn{p} the number of
#' structural (non-error-model) parameters.
#'
#' Parallel evaluation over conditions is supported via the `cores` argument
#' of the returned objective function. Parallelization is only applied if
#' `cores > 1`. A value of `cores <= 0` is not allowed.
#'
#' @example inst/examples/normL2.R
#' @export
normL2 <- function(data, x, errmodel = NULL, times = NULL,
                   attr.name = "data",
                   use.bessel = ifelse(!is.null(errmodel), TRUE, FALSE)) {
  
  ## --- time grid ------------------------------------------------------------
  timesD <- sort(unique(c(0, do.call(c, lapply(data, function(d) d$time)))))
  if (!is.null(times)) {
    timesD <- sort(union(times, timesD))
  }
  
  ## --- conditions -----------------------------------------------------------
  x.conditions <- names(attr(x, "mappings"))
  data.conditions <- names(data)
  
  if (!all(data.conditions %in% x.conditions)) {
    stop("The prediction function does not provide predictions for all conditions in the data.")
  }
  
  e.conditions <- if (!is.null(errmodel))
    names(attr(errmodel, "mappings"))
  else
    NULL
  
  ## Ensure errmodel is captured in closure (e.g. for runbg)
  force(errmodel)
  
  ## --- Bessel correction ----------------------------------------------------
  if (use.bessel && !is.null(errmodel)) {
    
    n.data <- sum(sapply(data, nrow))
    
    par.structural <- union(
      getParameters(x),
      getParameters(errmodel)
    )
    
    par.err <- setdiff(
      getSymbols(unlist(getEquations(errmodel))),
      names(unlist(getEquations(errmodel)))
    )
    
    p <- length(par.structural) - length(par.err)
    
    bessel.correction <- sqrt(n.data / (n.data - p))
    
  } else {
    bessel.correction <- 1
  }
  
  controls <- list(
    times = timesD,
    attr.name = attr.name,
    conditions = intersect(x.conditions, data.conditions),
    bessel.correction = bessel.correction
  )
  
  ## ========================================================================
  ## Objective function
  ## ========================================================================
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE,
                   conditions = controls$conditions,
                   env = NULL, cores = 1) {
    
    if (!is.numeric(cores) || cores <= 0) {
      stop("`cores` must be a positive integer (cores >= 1).")
    }
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    if (is.null(env)) {
      env <- new.env()
    }
    
    ## Prediction for all conditions at once
    prediction <- x(
      times = controls$times,
      pars = pouter,
      fixed = fixed,
      deriv = deriv,
      deriv2 = deriv2,
      conditions = conditions
    )
    
    pars_pred <- if (length(prediction) > 0)
      getParameters(prediction[[1]])
    else
      NULL
    
    use_errmodel <- !is.null(errmodel) &&
      (is.null(e.conditions) || any(conditions %in% e.conditions))
    
    ## --- condition-wise evaluation -----------------------------------------
    eval_condition <- function(cn) {
      
      if (use_errmodel && (is.null(e.conditions) || cn %in% e.conditions)) {
        
        err <- errmodel(
          out = prediction[[cn]],
          pars = pars_pred,
          conditions = cn
        )
        
        err_cn <- if (!is.null(err) && !is.null(err[[cn]]))
          err[[cn]]
        else
          NULL
        
      } else {
        err_cn <- NULL
      }
      
      nll(
        res(data[[cn]], prediction[[cn]], err_cn),
        pars = pouter,
        deriv = deriv,
        bessel.correction = controls$bessel.correction
      )
    }
    
    ## --- apply over conditions ---------------------------------------------
    if (cores == 1) {
      
      out.data <- lapply(conditions, eval_condition)
      
    } else if (.Platform$OS.type != "windows") {
      
      out.data <- parallel::mclapply(
        conditions, eval_condition, mc.cores = cores
      )
      
    } else {
      
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      out.data <- parallel::parLapply(cl, conditions, eval_condition)
    }
    
    ## --- aggregate ----------------------------------------------------------
    out.data <- Reduce("+", out.data)
    
    out <- out.data
    attr(out, controls$attr.name) <- out.data$value
    assign("prediction", prediction, envir = env)
    attr(out, "env") <- env
    
    out
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- data.conditions
  attr(myfn, "parameters") <- attr(x, "parameters")
  attr(myfn, "modelname") <- modelname(x, errmodel)
  
  myfn
}


#' Soft L2 constraint on parameters
#' 
#' @param mu named numeric, the prior values
#' @param sigma named numeric of length of mu or numeric of length one
#' or character of length of mu or character of length one
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return object of class \code{objfn}
#' @seealso \link{wrss}
#' @details If sigma is numeric, the function computes the constraint value 
#' \deqn{\left(\frac{p-\mu}{\sigma}\right)^2}{(p-mu)^2/sigma^2}
#' and its derivatives with respect to p. If sigma is a character, the 
#' function computes
#' \deqn{\left(\frac{p-\mu}{\sigma}\right)^2 + \log(\sigma^2)}{(p-mu)^2/sigma^2 + log(sigma^2)}
#' and its derivatives with respect to p and sigma. Sigma parameters being
#' passed to the function are ALWAYS assumed to be on a log scale, i.e. internally
#' sigma parameters are converted by \code{exp()}.
#' @export
constraintL2 <- function(mu, sigma = 1, attr.name = "prior", condition = NULL) {
  
  estimateSigma <- is.character(sigma)
  
  if (length(sigma) == 1) 
    sigma <- structure(rep(sigma, length(mu)), names = names(mu))
  if (is.null(names(sigma)))
    names(sigma) <- names(mu)
  sigma <- sigma[names(mu)]
  
  controls <- list(mu = mu, sigma = sigma, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = condition, env = NULL) {
    
    pouter <- list(...)[[match.fnargs(list(...), "pars")]]
    mu <- controls$mu
    sigma <- controls$sigma
    
    # Handle list input (multiple conditions)
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      if (length(available) == 0) return()
      pouter <- pouter[available]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      # Get deriv attributes before any manipulation
      dP <- attr(p, "deriv", exact = TRUE)
      dP2 <- attr(p, "deriv2", exact = TRUE)
      
      # Combine and extract available parameters
      all_pars <- c(as.numeric(p), fixed)
      names(all_pars) <- c(names(p), names(fixed))
      avail <- intersect(names(mu), names(all_pars))
      
      if (length(avail) == 0)
        return(objlist(value = 0, 
                       gradient = if(deriv) setNames(rep(0, length(p)), names(p)) else NULL,
                       hessian = if(deriv) matrix(0, length(p), length(p), dimnames = list(names(p), names(p))) else NULL))
      
      pars <- all_pars[avail]
      mu_a <- mu[avail]
      sig_a <- sigma[avail]
      n_a <- length(avail)
      
      p1 <- intersect(setdiff(names(mu), names(fixed)), names(p))
      p2 <- character(0)
      
      # Handle sigma estimation
      if (estimateSigma) {
        sig_a <- exp(all_pars[sig_a])
        names(sig_a) <- avail
        p2 <- intersect(setdiff(unique(sigma[avail]), names(fixed)), names(p))
      }
      
      # Compute value
      res <- pars - mu_a
      val <- sum(res^2 / sig_a^2) + estimateSigma * sum(2 * log(sig_a))
      
      gr <- hs <- NULL
      if (deriv) {
        gr <- setNames(rep(0, length(p)), names(p))
        hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
        
        p1_a <- intersect(p1, avail)
        if (length(p1_a) > 0) {
          gr[p1_a] <- 2 * res[p1_a] / sig_a[p1_a]^2
          diag(hs)[p1_a] <- 2 / sig_a[p1_a]^2
        }
        
        if (estimateSigma && length(p2) > 0) {
          # Aggregate sigma derivatives by sigma parameter name
          for (sp in p2) {
            idx <- which(sigma[avail] == sp)
            gr[sp] <- sum(-2 * res[idx]^2 / sig_a[idx]^2 + 2)
            hs[sp, sp] <- sum(4 * res[idx]^2 / sig_a[idx]^2)
          }
          # Cross terms p1 x p2
          for (sp in p2) {
            idx <- names(sigma[avail])[sigma[avail] == sp]
            common <- intersect(idx, p1_a)
            if (length(common) > 0) {
              hs[common, sp] <- -4 * res[common] / sig_a[common]^2
              hs[sp, common] <- hs[common, sp]
            }
          }
        }
        
        # Chain rule
        if (!is.null(dP)) {
          gr_inner <- gr
          gr <- as.vector(gr_inner %*% dP)
          names(gr) <- colnames(dP)
          hs <- t(dP) %*% hs %*% dP
          if (!is.null(dP2)) {
            hs <- hs + colSums(gr_inner * matrix(dP2, nrow = length(gr_inner)))
            dim(hs) <- c(ncol(dP), ncol(dP))
          }
          dimnames(hs) <- list(colnames(dP), colnames(dP))
        }
      }
      
      objlist(value = val, gradient = gr, hessian = hs)
    })
    
    out <- Reduce("+", outlist)
    attr(out, controls$attr.name) <- out$value
    attr(out, "env") <- env
    out
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  myfn
}


#' L2 objective function for validation data point
#' 
#' @param name character, the name of the prediction, e.g. a state name.
#' @param time numeric, the time-point associated to the prediction
#' @param value character, the name of the parameter which contains the
#' prediction value.
#' @param sigma numeric, the uncertainty of the introduced test data point
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the prediction is made.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}, \link{constraintL2}
#' @details Computes the constraint value 
#' \deqn{\left(\frac{x(t)-\mu}{\sigma}\right)^2}{(pred-p[names(mu)])^2/sigma^2}
#' and its derivatives with respect to p.
#' @examples
#' prediction <- list(a = matrix(c(0, 1), nrow = 1, dimnames = list(NULL, c("time", "A"))))
#' derivs <- array(c(1, 0.1), dim = c(1, 1, 2), dimnames = list(NULL, "A", c("A", "k1")))
#' attr(prediction$a, "deriv") <- derivs
#' p0 <- c(A = 1, k1 = 2)
#' 
#' vali <- datapointL2(name = "A", time = 0, value = "newpoint", sigma = 1, condition = "a")
#' vali(pars = c(p0, newpoint = 1), env = .GlobalEnv)
#' @export
datapointL2 <- function(name, time, value, sigma = 1, attr.name = "validation", condition) {
  
  
  controls <- list(
    mu = structure(name, names = value)[1], # Only one data point is allowed
    time = time[1],
    sigma = sigma[1],
    attr.name = attr.name
  )
  
  
  myfn <- function(..., fixed = NULL, deriv = TRUE, deriv2 = FALSE, conditions = NULL, env = NULL) {
    
    # Import from controls
    mu <- controls$mu
    time <- controls$time
    sigma <- controls$sigma
    attr.name <- controls$attr.name
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    if (!is.null(env)) {
      prediction <- as.list(env)$prediction 
    } else {
      stop("No prediction available. Use the argument env to pass an environment that contains the prediction.")
    }
    
    # Return result only when the data point condition overlaps with the evaluated conditions
    if (!is.null(conditions) && !condition %in% conditions) 
      return()
    if (is.null(conditions) && !condition %in% names(prediction))
      stop("datapointL2 requests unavailable condition. Call the objective function explicitly stating the conditions argument.")
    
    
    
    # Divide parameter into data point and rest
    datapar <- setdiff(names(mu), names(fixed))
    parapar <- setdiff(names(pouter), c(datapar, names(fixed)))
    n_pouter <- length(pouter)
    
    
    # Get predictions and derivatives at time point
    time.index <- which(prediction[[condition]][, "time"] == time)
    if (length(time.index) == 0) 
      stop("datapointL2() requests time point for which no prediction is available. Please add missing time point by the times argument in normL2()")
    
    # Get prediction value - check if column exists
    pred <- tryCatch(
      prediction[[condition]][time.index, mu],
      error = function(e) {
        available_states <- setdiff(colnames(prediction[[condition]]), "time")
        stop(sprintf(
          "datapointL2() cannot find state '%s' in prediction.\n  Available states: %s\n  If '%s' is an input, make sure to use Y() with attach.input = TRUE!",
          mu, 
          paste(available_states, collapse = ", "),
          mu
        ), call. = FALSE)
      }
    )
    
    # Get first derivatives (deriv is 3D array [time, name, param])
    deriv_array <- attr(prediction[[condition]], "deriv")
    withDeriv <- !is.null(deriv_array)
    dxdp <- NULL
    if (withDeriv) {
      # Extract derivatives: deriv_array[time.index, mu, ]
      dxdp_full <- tryCatch(
        deriv_array[time.index, mu, , drop = TRUE],
        error = function(e) {
          available_states <- if (!is.null(dimnames(deriv_array)[[2]])) 
            dimnames(deriv_array)[[2]] 
          else 
            "unknown"
          stop(sprintf(
            "datapointL2() cannot find derivatives for state '%s'.\n  Available states: %s\n  If '%s' is an input, make sure to use Y() with attach.input = TRUE!",
            mu, 
            paste(available_states, collapse = ", "),
            mu
          ), call. = FALSE)
        }
      )
      
      # Only extract parameters that actually exist in the derivatives
      if (length(parapar) > 0) {
        available_parapar <- intersect(parapar, names(dxdp_full))
        dxdp <- numeric(length(parapar))
        names(dxdp) <- parapar
        if (length(available_parapar) > 0) {
          dxdp[available_parapar] <- dxdp_full[available_parapar]
        }
        # Parameters not in dxdp_full remain 0 (not NA)
      } else {
        dxdp <- numeric(0)
      }
    }
    
    # Get second derivatives (deriv2 is 4D array [time, name, param, param])
    deriv2_array <- attr(prediction[[condition]], "deriv2")
    withDeriv2 <- !is.null(deriv2_array) && deriv2
    d2xdp2 <- NULL
    if (withDeriv2 && length(parapar) > 0) {
      # Extract second derivatives: deriv2_array[time.index, mu, , ]
      d2xdp2_full <- tryCatch(
        deriv2_array[time.index, mu, , , drop = TRUE],
        error = function(e) {
          available_states <- if (!is.null(dimnames(deriv2_array)[[2]])) 
            dimnames(deriv2_array)[[2]] 
          else 
            "unknown"
          stop(sprintf(
            "datapointL2() cannot find second derivatives for state '%s'.\n  Available states: %s\n  If '%s' is an input, make sure to use Y() with attach.input = TRUE!",
            mu, 
            paste(available_states, collapse = ", "),
            mu
          ), call. = FALSE)
        }
      )
      
      # Only extract parameters that actually exist in the second derivatives
      available_parapar <- intersect(parapar, rownames(d2xdp2_full))
      d2xdp2 <- matrix(0, length(parapar), length(parapar), 
                       dimnames = list(parapar, parapar))
      if (length(available_parapar) > 0) {
        d2xdp2[available_parapar, available_parapar] <- 
          d2xdp2_full[available_parapar, available_parapar, drop = FALSE]
      }
      # Parameters not in d2xdp2_full get 0 entries (not NA)
    }
    
    # Compute residual and objective value
    res <- as.numeric(pred - c(fixed, pouter)[names(mu)])
    sigma2_inv <- 1/(sigma * sigma)
    val <- res * res * sigma2_inv
    
    # Initialize gradient and hessian
    gr <- NULL
    hs <- NULL
    
    if (deriv && withDeriv) {
      # Derivative of residual with respect to parameters
      dres.dp <- numeric(n_pouter)
      names(dres.dp) <- names(pouter)
      if (length(parapar) > 0) 
        dres.dp[parapar] <- dxdp
      if (length(datapar) > 0) 
        dres.dp[datapar] <- -1
      
      # Gradient: d(val)/dp = 2 * res * dres.dp / sigma^2
      gr <- 2 * res * sigma2_inv * dres.dp
      
      # Hessian: d2(val)/dp2 = 2 * dres.dp %*% t(dres.dp) / sigma^2 + 2 * res * d2res.dp2 / sigma^2
      # Use tcrossprod for efficiency (C-based BLAS operation)
      hs <- 2 * sigma2_inv * tcrossprod(dres.dp)
      colnames(hs) <- rownames(hs) <- names(pouter)
      
      # Add second derivative contribution if available
      if (withDeriv2 && length(parapar) > 0) {
        hs[parapar, parapar] <- hs[parapar, parapar] + 2 * res * sigma2_inv * d2xdp2
      }
    }
    
    out <- objlist(value = val, gradient = gr, hessian = hs)
    attr(out, controls$attr.name) <- out$value
    attr(out, "prediction") <- pred
    
    attr(out, "env") <- env
    
    return(out)
    
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- value[1]
  
  
  return(myfn)
  
  
}


#' L2 objective function for prior value
#' 
#' @description As a prior function, it returns derivatives with respect to
#' the penalty parameter in addition to parameter derivatives.
#' 
#' @param mu Named numeric, the prior values
#' @param lambda Character of length one. The name of the penalty paramter in \code{p}.
#' @param attr.name character. The constraint value is additionally returned in an 
#' attributed with this name
#' @param condition character, the condition for which the constraint should apply. If
#' \code{NULL}, applies to any condition.
#' @return List of class \code{objlist}, i.e. objective value, gradient and Hessian as list.
#' @seealso \link{wrss}
#' @details Computes the constraint value 
#' \deqn{e^{\lambda} \| p-\mu \|^2}{exp(lambda)*sum((p-mu)^2)}
#' and its derivatives with respect to p and lambda.
#' @examples
#' p <- c(A = 1, B = 2, C = 3, lambda = 0)
#' mu <- c(A = 0, B = 0)
#' obj <- priorL2(mu = mu, lambda = "lambda")
#' obj(pars = p + rnorm(length(p), 0, .1))
#' @export
priorL2 <- function(mu, lambda = "lambda", attr.name = "prior", condition = NULL) {
  
  
  controls <- list(mu = mu, lambda = lambda, attr.name = attr.name)
  
  myfn <- function(..., fixed = NULL, deriv=TRUE, conditions = condition, env = NULL) {
    
    arglist <- list(...)
    arglist <- arglist[match.fnargs(arglist, "pars")]
    pouter <- arglist[[1]]
    
    # Import from controls 
    mu <- controls$mu
    lambda <- controls$lambda
    attr.name <- controls$attr.name
    
    # pouter can be a list (if result from a parameter transformation)
    # In this case match with conditions and evaluate only those
    # If there is no overlap, return NULL
    # If pouter is not a list, evaluate the constraint function 
    # for this pouter.
    
    if (is.list(pouter) && !is.null(conditions)) {
      available <- intersect(names(pouter), conditions)
      defined <- ifelse(is.null(condition), TRUE, condition %in% conditions)
      
      if (length(available) == 0 | !defined) return()
      pouter <- pouter[intersect(available, condition)]
    }
    if (!is.list(pouter)) pouter <- list(pouter)
    
    outlist <- lapply(pouter, function(p) {
      
      
      ## Extract contribution of fixed pars and delete names for calculation of gr and hs  
      par.fixed <- intersect(names(mu), names(fixed))
      sumOfFixed <- 0
      if (!is.null(par.fixed)) sumOfFixed <- sum(exp(c(fixed, p)[lambda])*(fixed[par.fixed] - mu[par.fixed]) ^ 2)
      
      # Compute prior value and derivatives
      par <- intersect(names(mu), names(p))
      par0 <- setdiff(par, lambda)
      
      val <- sum(exp(c(fixed, p)[lambda]) * (p[par] - mu[par]) ^ 2) + sumOfFixed
      
      gr <- hs <- NULL
      if (deriv) {
        gr <- rep(0, length(p)); names(gr) <- names(p)
        gr[par] <- 2*exp(c(fixed, p)[lambda])*(p[par] - mu[par])
        if (lambda %in% names(p)) {
          gr[lambda] <- sum(exp(c(fixed, p)[lambda]) * (p[par0] - mu[par0]) ^ 2) + 
            sum(exp(c(fixed, p)[lambda]) * (fixed[par.fixed] - mu[par.fixed]) ^ 2)
        }
        
        hs <- matrix(0, length(p), length(p), dimnames = list(names(p), names(p)))
        diag(hs)[par] <- 2*exp(c(fixed, p)[lambda])
        if (lambda %in% names(p)) {
          hs[lambda, lambda] <- gr[lambda] 
          hs[lambda, par0] <- hs[par0, lambda] <- gr[par0]
        }
        
        dP <- attr(p, "deriv")
        if (!is.null(dP)) {
          gr <- as.vector(gr %*% dP); names(gr) <- colnames(dP)
          hs <- t(dP) %*% hs %*% dP; colnames(hs) <- colnames(dP); rownames(hs) <- colnames(dP)
        }
      }
      
      objlist(value = val, gradient = gr, hessian = hs)
      
    })
    
    out <- Reduce("+", outlist)
    attr(out, controls$attr.name) <- out$value
    attr(out, "env") <- env
    
    return(out)
    
    
  }
  
  class(myfn) <- c("objfn", "fn")
  attr(myfn, "conditions") <- condition
  attr(myfn, "parameters") <- names(mu)
  return(myfn)
  
}


#' Compute the negative log-likelihood
#' 
#' @description Gaussian Log-likelihood. Supports NONMEM-like BLOQ handling methods M1, M3 and M4 
#' and estimation of error models with optional Bessel correction for variance parameter bias.
#' 
#' @param nout data.frame (result of [res]) or object of class [res].
#' @param pars Named vector of outer parameters to construct the objlist
#' @param deriv Logical. If TRUE, compute gradient and hessian
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data. 
#' One of "M1", "M3", "M4NM", or "M4BEAL".
#' @param opt.hessian Named logical vector to include or exclude various 
#' summands of the hessian matrix. Controls which parts contribute to the 
#' Hessian calculation for both ALOQ and BLOQ data.
#' @param bessel.correction Numeric. Bessel correction factor for variance estimation.
#' Default is 1 (no correction). When use.bessel = TRUE in normL2(), this is set to 
#' sqrt(n/(n-p)) where n is total data points and p is number of structural (non-errormodel) parameters.
#' 
#' @md
#' @return list with entries value (numeric, the negative log-likelihood with Bessel correction applied), 
#' gradient (numeric vector, gradient with respect to parameters) and 
#' hessian (matrix of type numeric, second derivatives). The returned value includes the Bessel correction
#' and is used for optimization. The TRUE maximum likelihood (for BIC/AIC) is stored in attr(*, "neg2ll").
#' 
#' @details The Bessel correction is applied by multiplying weighted residuals with bessel.correction.
#' This changes both the objective value AND the gradient/Hessian, leading to corrected parameter 
#' estimates during optimization. The TRUE maximum likelihood value (without correction) is stored 
#' in the "neg2ll" attribute for use in BIC/AIC calculations.
#' 
#' @export
nll <- function(nout, pars, deriv, opt.BLOQ = "M3", opt.hessian = c(
  ALOQ_part1 = TRUE, ALOQ_part2 = TRUE, ALOQ_part3 = TRUE,
  BLOQ_part1 = TRUE, BLOQ_part2 = TRUE, BLOQ_part3 = TRUE,
  PD = TRUE), bessel.correction = 1) {
  
  # Split residuals into ALOQ and BLOQ
  is.bloq   <- nout$bloq
  nout.bloq <- nout[is.bloq, , drop = FALSE]
  nout.aloq <- nout[!is.bloq, , drop = FALSE]
  
  # Handle derivs (matrices and arrays)
  derivs          <- attr(nout, "deriv")
  derivs.bloq     <- if (!is.null(derivs)) derivs[is.bloq, , drop = FALSE] else NULL
  derivs.aloq     <- if (!is.null(derivs)) derivs[!is.bloq, , drop = FALSE] else NULL
  
  derivs.err      <- attr(nout, "deriv.err")
  derivs.err.bloq <- if (!is.null(derivs.err)) derivs.err[is.bloq, , drop = FALSE] else NULL
  derivs.err.aloq <- if (!is.null(derivs.err)) derivs.err[!is.bloq, , drop = FALSE] else NULL
  
  deriv2          <- attr(nout, "deriv2")
  deriv2.bloq     <- if (!is.null(deriv2)) deriv2[is.bloq, , , drop = FALSE] else NULL
  deriv2.aloq     <- if (!is.null(deriv2)) deriv2[!is.bloq, , , drop = FALSE] else NULL
  
  deriv2.err      <- attr(nout, "deriv2.err")
  deriv2.err.bloq <- if (!is.null(deriv2.err)) deriv2.err[is.bloq, , , drop = FALSE] else NULL
  deriv2.err.aloq <- if (!is.null(deriv2.err)) deriv2.err[!is.bloq, , , drop = FALSE] else NULL
  
  mywrss <- {
    gr <- if (deriv) setNames(numeric(length(pars)), names(pars)) else NULL
    he <- if (deriv) matrix(
      0, length(pars), length(pars),
      dimnames = list(names(pars), names(pars))
    ) else NULL
    
    objlist(value = 0, gradient = gr, hessian = he)
  }
  
  # Apply nll
  nll_ALOQ_result <- NULL
  
  if (!all(is.bloq)) {
    nll_ALOQ_result <- nll_ALOQ(
      nout.aloq, derivs.aloq, derivs.err.aloq, 
      deriv2.aloq, deriv2.err.aloq,
      opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian, 
      bessel.correction = bessel.correction
    )
  }
  
  mywrss <- mywrss + nll_ALOQ_result
  
  if (any(is.bloq) && (!opt.BLOQ == "M1")) {
    mywrss <- mywrss + nll_BLOQ(
      nout.bloq, derivs.bloq, derivs.err.bloq,
      deriv2.bloq, deriv2.err.bloq,
      opt.BLOQ = opt.BLOQ, opt.hessian = opt.hessian
    )
  }
  
  # Extract and store attributes
  chisquare <- attr(nll_ALOQ_result, "chisquare")
  nll <- attr(nll_ALOQ_result, "nll")
  attr(mywrss, "chisquare") <- if (length(chisquare)) chisquare else 0
  attr(mywrss, "nll") <- if (length(nll)) nll else 0 
  
  mywrss
}


#' Non-linear log likelihood for the ALOQ part of the data
#' 
#' @param nout output of [res()]
#' @param derivs,derivs.err matrix of first derivatives
#' @param deriv2,deriv2.err array of second derivatives
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector to include or exclude
#'   various non-convex summands of the hessian matrix
#' @param bessel.correction Numeric. Bessel correction factor applied to the 
#'   weighted residuals. Default is 1 (no correction).
#'   
#' @md
#' @importFrom stats pnorm dnorm
nll_ALOQ <- function(nout,
                     derivs,
                     derivs.err,
                     deriv2 = NULL,
                     deriv2.err = NULL,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(
                       ALOQ_part1 = TRUE,
                       ALOQ_part2 = TRUE,
                       ALOQ_part3 = TRUE
                     ),
                     bessel.correction = 1) {
  
  # .. Residual terms ----#
  wr <- nout$weighted.residual
  w0 <- nout$weighted.0
  s  <- nout$sigma
  
  # .. Compute TRUE negative log-likelihood (for BIC/AIC) ----#
  chisquare_ml <- sum(wr^2)
  neg2ll_ml <- chisquare_ml + sum(log(2*pi*s^2))
  
  # .. Apply Bessel correction to weighted residuals ----#
  use_bessel <- bessel.correction != 1
  if (use_bessel) {
    wr <- wr * bessel.correction
    w0 <- w0 * bessel.correction
  }
  
  # .. Compute corrected objective value (for optimization) ----#
  chisquare <- sum(wr^2)
  obj <- chisquare + sum(log(2*pi*s^2))
  
  if (opt.BLOQ %in% "M4BEAL") {
    bloq_term <- 2 * sum(stats::pnorm(w0, log.p = TRUE))
    obj <- obj + bloq_term
    neg2ll_ml <- neg2ll_ml + bloq_term
  }
  
  grad <- NULL
  hessian <- NULL
  
  if (!is.null(derivs) && nrow(derivs) > 0) {
    # .. Sensitivities terms ----#
    dxdp <- derivs
    dsdp <- if (!is.null(derivs.err)) derivs.err else matrix(0, nrow = nrow(derivs), ncol = ncol(derivs))
    n_pars <- ncol(dxdp)
    
    # Compute base derivatives
    dwrdp <- 1/s * dxdp - wr/s * dsdp
    dw0dp <- 1/s * dxdp - w0/s * dsdp
    
    # Apply Bessel correction to derivatives if needed
    if (use_bessel) {
      dwrdp <- dwrdp * bessel.correction
      dw0dp <- dw0dp * bessel.correction
    }
    
    dlogsdp <- (1/s) * dsdp
    G_by_Phi <- function(w) exp(stats::dnorm(w, log = TRUE) - stats::pnorm(w, log.p = TRUE))
    
    # .. Compute gradient ----#
    grad <- as.vector(2 * colSums(wr * dwrdp) + 2 * colSums(dlogsdp))
    if (opt.BLOQ %in% "M4BEAL") {
      grad <- grad + as.vector(2 * colSums(G_by_Phi(w0) * dw0dp))
    }
    names(grad) <- colnames(dxdp)
    
    # .. Compute hessian ----#
    hessian <- matrix(0, nrow = n_pars, ncol = n_pars, 
                      dimnames = list(colnames(dxdp), colnames(dxdp)))
    
    # Base term: 2 * dwrdp^T * dwrdp
    hessian <- hessian + 2 * crossprod(dwrdp)
    
    if (opt.hessian["ALOQ_part1"]) {
      # Interaction terms: 2 * (-wr/s^2) * (dxdp * dsdp^T + dsdp * dxdp^T)
      hessian <- hessian + 2 * crossprod(-wr/s^2 * dxdp, dsdp)
      hessian <- hessian + 2 * crossprod(-wr/s^2 * dsdp, dxdp)
    }
    
    if (opt.hessian["ALOQ_part2"]) {
      # 2 * (2*wr^2/s^2) * dsdp^T * dsdp
      hessian <- hessian + 2 * crossprod(sqrt(2) * wr/s * dsdp)
    }
    
    if (opt.hessian["ALOQ_part3"]) {
      # Non-convex log(sigma) contribution: -2 * dlogsdp^T * dlogsdp
      hessian <- hessian - 2 * crossprod(dlogsdp)
    }
    
    # Add second derivative terms if available
    if (!is.null(deriv2)) {
      # deriv2 is array [n_data, n_pars, n_pars]
      # einsum("i,ijk->jk", weights, arr) = crossprod of weights with flattened arr
      
      d2xdp2 <- deriv2
      d2sdp2 <- if (!is.null(deriv2.err)) deriv2.err else array(0, dim = dim(deriv2))
      
      n_i <- dim(d2xdp2)[1]
      
      # Weights for the two terms
      w1 <- wr * (1/s)
      w2 <- wr * (wr/s)
      
      # einsum("i,ijk->jk", w, A) = t(w) %*% matrix(A, n_i, n_j*n_k), reshaped
      # Reshape d2xdp2 to [n_i, n_pars * n_pars]
      d2xdp2_flat <- matrix(d2xdp2, nrow = n_i, ncol = n_pars * n_pars)
      d2sdp2_flat <- matrix(d2sdp2, nrow = n_i, ncol = n_pars * n_pars)
      
      contrib_flat <- crossprod(w1, d2xdp2_flat) - crossprod(w2, d2sdp2_flat)
      contrib_mat <- matrix(contrib_flat, nrow = n_pars, ncol = n_pars)
      
      hessian <- hessian + 2 * contrib_mat
    }
    
    # M4BEAL terms
    if (opt.BLOQ %in% "M4BEAL") {
      G_w0 <- G_by_Phi(w0)
      
      # First order terms: 2 * (-w0*G - G^2) * dw0dp^T * dw0dp
      hessian <- hessian + 2 * crossprod(sqrt(pmax(0, -w0 * G_w0 - G_w0^2)) * dw0dp)
      
      # Second order terms: 2 * G * (-1/s^2) * (dxdp * dsdp^T + dsdp * dxdp^T)
      hessian <- hessian + 2 * crossprod(G_w0 * (-1/s^2) * dxdp, dsdp)
      hessian <- hessian + 2 * crossprod(G_w0 * (-1/s^2) * dsdp, dxdp)
      
      if (opt.hessian["ALOQ_part1"]) {
        # 2 * G * (2*w0/s^2) * dsdp^T * dsdp
        hessian <- hessian + 2 * crossprod(sqrt(2 * G_w0 * w0/s^2) * dsdp)
      }
      
      # Add second derivative terms for M4BEAL if available
      if (!is.null(deriv2)) {
        d2xdp2 <- deriv2
        d2sdp2 <- if (!is.null(deriv2.err)) deriv2.err else array(0, dim = dim(deriv2))
        
        n_i <- dim(d2xdp2)[1]
        
        w1 <- G_w0 * (1/s)
        w2 <- G_w0 * (w0/s)
        
        d2xdp2_flat <- matrix(d2xdp2, nrow = n_i, ncol = n_pars * n_pars)
        d2sdp2_flat <- matrix(d2sdp2, nrow = n_i, ncol = n_pars * n_pars)
        
        contrib_flat <- crossprod(w1, d2xdp2_flat) - crossprod(w2, d2sdp2_flat)
        contrib_mat <- matrix(contrib_flat, nrow = n_pars, ncol = n_pars)
        
        hessian <- hessian + 2 * contrib_mat
      }
    }
  }
  
  out <- objlist(value = obj, gradient = grad, hessian = hessian)
  attr(out, "chisquare") <- chisquare_ml
  attr(out, "nll") <- neg2ll_ml
  attr(out, "besselcorrected") <- use_bessel
  out
}



#' Non-linear log likelihood for the BLOQ part of the data
#' @md
#' @param nout.bloq The bloq output of [res()]
#' @param derivs.bloq,derivs.err.bloq matrix of first derivatives
#' @param deriv2.bloq,deriv2.err.bloq array of second derivatives
#' @param opt.BLOQ Character denoting the method to deal with BLOQ data
#' @param opt.hessian Named logical vector to include or exclude
#'   various summands of the hessian matrix
#' @importFrom stats pnorm dnorm
nll_BLOQ <- function(nout.bloq,
                     derivs.bloq,
                     derivs.err.bloq,
                     deriv2.bloq = NULL,
                     deriv2.err.bloq = NULL,
                     opt.BLOQ = c("M3", "M4NM", "M4BEAL", "M1"),
                     opt.hessian = c(
                       BLOQ_part1 = TRUE,
                       BLOQ_part2 = TRUE,
                       BLOQ_part3 = TRUE
                     )) {
  
  # .. Checks -----
  if (opt.BLOQ %in% c("M4NM", "M4BEAL") & any(nout.bloq$value < 0)) {
    stop("M4-Method cannot handle LLOQ < 0. Possible solutions:\n",
         "  * Use M3 which allows negative LLOQ (recommended)\n",
         "  * If you are working with log-transformed DV, exponentiate DV and LLOQ\n")
  }
  
  # .. Residuals ----#
  wr <- nout.bloq$weighted.residual
  w0 <- nout.bloq$weighted.0
  s  <- nout.bloq$sigma
  
  # .. Compute objective value ----#
  if (opt.BLOQ == "M3") {
    objvals.bloq <- -2 * stats::pnorm(-wr, log.p = TRUE)
  }
  
  if (opt.BLOQ %in% c("M4NM", "M4BEAL")) {
    objvals.bloq <- -2 * log(1 - stats::pnorm(wr) / stats::pnorm(w0))
    
    # Catch numerically problematic cases
    intercept <- ifelse(log(w0 - wr) > 0, 1.8, -1.9 * log(w0 - wr) + 0.9)
    lin <- ifelse(log(w0 - wr) > 0, 0.9, 0.5)
    objvals.bloq[!is.finite(objvals.bloq)] <- (intercept + lin * w0 + 0.95 * w0^2)[!is.finite(objvals.bloq)]
  }
  
  obj.bloq <- sum(objvals.bloq)
  grad.bloq <- NULL
  hessian.bloq <- NULL
  
  if (!is.null(derivs.bloq) && nrow(derivs.bloq) > 0) {
    # .. Sensitivities ----#
    dxdp <- derivs.bloq
    dsdp <- if (!is.null(derivs.err.bloq)) derivs.err.bloq else matrix(0, nrow = nrow(derivs.bloq), ncol = ncol(derivs.bloq))
    n_pars <- ncol(dxdp)
    
    dwrdp <- 1/s * dxdp - wr/s * dsdp
    dw0dp <- 1/s * dxdp - w0/s * dsdp
    dlogsdp <- (1/s) * dsdp
    
    G_by_Phi <- function(w1, w2 = w1) {
      exp(stats::dnorm(w1, log = TRUE) - stats::pnorm(w2, log.p = TRUE))
    }
    
    # .. Compute gradient ----#
    if (opt.BLOQ == "M3") {
      grad.bloq <- 2 * as.vector(colSums(G_by_Phi(-wr) * dwrdp))
    }
    
    if (opt.BLOQ %in% c("M4NM", "M4BEAL")) {
      c1 <- 1 / (1/G_by_Phi(wr, w0) - 1/G_by_Phi(wr, wr))
      c2 <- 1 / (1/G_by_Phi(w0, w0) - 1/G_by_Phi(w0, wr))
      c3 <- G_by_Phi(w0)
      
      grad.bloq <- as.vector(2 * colSums(c1 * dwrdp - c2 * dw0dp + c3 * dw0dp))
    }
    names(grad.bloq) <- colnames(dxdp)
    
    # .. Compute hessian ----
    hessian.bloq <- matrix(0, nrow = n_pars, ncol = n_pars, 
                           dimnames = list(colnames(dxdp), colnames(dxdp)))
    
    if (opt.BLOQ %in% "M3") {
      G_neg_wr <- G_by_Phi(-wr)
      
      if (opt.hessian["BLOQ_part1"]) {
        # 2 * (-wr*G + G^2) * dwrdp^T * dwrdp
        coef <- -wr * G_neg_wr + G_neg_wr^2
        # Handle sign: coef can be positive or negative
        hessian.bloq <- hessian.bloq + 2 * crossprod(dwrdp, coef * dwrdp)
      }
      
      if (opt.hessian["BLOQ_part2"]) {
        # -2 * G * (+1/s^2) * (dxdp * dsdp^T + dsdp * dxdp^T)
        hessian.bloq <- hessian.bloq - 2 * crossprod(G_neg_wr * (1/s^2) * dxdp, dsdp)
        hessian.bloq <- hessian.bloq - 2 * crossprod(G_neg_wr * (1/s^2) * dsdp, dxdp)
      }
      
      if (opt.hessian["BLOQ_part3"]) {
        # -2 * G * (2*(-wr)/s^2) * dsdp^T * dsdp
        hessian.bloq <- hessian.bloq - 2 * crossprod(dsdp, G_neg_wr * (2 * (-wr)/s^2) * dsdp)
      }
      
      # Add second derivative terms if available
      if (!is.null(deriv2.bloq)) {
        d2xdp2 <- deriv2.bloq
        d2sdp2 <- if (!is.null(deriv2.err.bloq)) deriv2.err.bloq else array(0, dim = dim(deriv2.bloq))
        
        n_i <- dim(d2xdp2)[1]
        
        # For M3: d(-wr)/dp has sign flip, so d2(-wr)/dp2 = -d2xdp2/s + ...
        w1 <- G_neg_wr * (1/s)
        w2 <- G_neg_wr * (wr/s)
        
        d2xdp2_flat <- matrix(d2xdp2, nrow = n_i, ncol = n_pars * n_pars)
        d2sdp2_flat <- matrix(d2sdp2, nrow = n_i, ncol = n_pars * n_pars)
        
        # -einsum w1 * d2xdp2 + einsum w2 * d2sdp2
        contrib_flat <- -crossprod(w1, d2xdp2_flat) + crossprod(w2, d2sdp2_flat)
        contrib_mat <- matrix(contrib_flat, nrow = n_pars, ncol = n_pars)
        
        hessian.bloq <- hessian.bloq + 2 * contrib_mat
      }
    }
    
    # .. M4 hessian ----
    if (opt.BLOQ %in% c("M4NM", "M4BEAL")) {
      # Helper function for stable division
      stable <- function(wn, w0, wr) {
        out <- stats::dnorm(wn) / (stats::pnorm(w0) - stats::pnorm(wr))
        
        if (identical(wn, w0)) {
          out[is.infinite(out)] <- 0
          return(out)
        }
        if (identical(wn, wr)) {
          out[is.infinite(out)] <- 1/(w0 - wr) + wr
          return(out)
        }
        out
      }
      
      # Coefficient vectors
      A1 <- -wr * stable(wr, w0, wr)
      A2 <- stable(wr, w0, wr)
      A3 <- -w0 * stable(w0, w0, wr)
      A4 <- stable(w0, w0, wr)
      A5 <- -w0 * G_by_Phi(w0) - G_by_Phi(w0)^2
      A6 <- G_by_Phi(w0)
      
      # Part 1: First order derivative terms
      part1 <- 2 * (
        crossprod(dwrdp, A1 * dwrdp) +
          crossprod(dw0dp, A3 * dw0dp) +
          crossprod(dw0dp, A5 * dw0dp)
      )
      
      # Part 2: Cross terms (non-convex)
      part2_vec <- A2 * dwrdp - A4 * dw0dp
      part2 <- -2 * crossprod(part2_vec)
      
      # Part 3: Second order derivative interaction terms
      part3 <- matrix(0, n_pars, n_pars)
      if (opt.hessian["BLOQ_part2"]) {
        part3 <- part3 + 2 * (
          crossprod(A2 * (-1/s^2) * dxdp, dsdp) +
            crossprod(A2 * (-1/s^2) * dsdp, dxdp) +
            crossprod(A4 * (-1/s^2) * dxdp, dsdp) +
            crossprod(A4 * (-1/s^2) * dsdp, dxdp) +
            crossprod(A6 * (-1/s^2) * dxdp, dsdp) +
            crossprod(A6 * (-1/s^2) * dsdp, dxdp)
        )
      }
      
      if (opt.hessian["BLOQ_part3"]) {
        part3 <- part3 + 2 * (
          crossprod(dsdp, A2 * (2 * wr/s^2) * dsdp) +
            crossprod(dsdp, A4 * (2 * w0/s^2) * dsdp) +
            crossprod(dsdp, A6 * (2 * w0/s^2) * dsdp)
        )
      }
      
      if (opt.hessian["BLOQ_part1"]) {
        hessian.bloq <- hessian.bloq + part1
      }
      if (opt.hessian["BLOQ_part2"]) {
        hessian.bloq <- hessian.bloq + part2
      }
      if (opt.hessian["BLOQ_part3"]) {
        hessian.bloq <- hessian.bloq + part3
      }
      
      # Add second derivative terms if available
      if (!is.null(deriv2.bloq)) {
        d2xdp2 <- deriv2.bloq
        d2sdp2 <- if (!is.null(deriv2.err.bloq)) deriv2.err.bloq else array(0, dim = dim(deriv2.bloq))
        
        n_i <- dim(d2xdp2)[1]
        
        d2xdp2_flat <- matrix(d2xdp2, nrow = n_i, ncol = n_pars * n_pars)
        d2sdp2_flat <- matrix(d2sdp2, nrow = n_i, ncol = n_pars * n_pars)
        
        # contrib_wr: A2 * ((1/s) * d2xdp2 - (wr/s) * d2sdp2)
        w1_wr <- A2 * (1/s)
        w2_wr <- A2 * (wr/s)
        contrib_wr_flat <- crossprod(w1_wr, d2xdp2_flat) - crossprod(w2_wr, d2sdp2_flat)
        contrib_wr <- matrix(contrib_wr_flat, nrow = n_pars, ncol = n_pars)
        
        # contrib_w0: (A4 + A6) * ((1/s) * d2xdp2 - (w0/s) * d2sdp2)
        w1_w0 <- (A4 + A6) * (1/s)
        w2_w0 <- (A4 + A6) * (w0/s)
        contrib_w0_flat <- crossprod(w1_w0, d2xdp2_flat) - crossprod(w2_w0, d2sdp2_flat)
        contrib_w0 <- matrix(contrib_w0_flat, nrow = n_pars, ncol = n_pars)
        
        hessian.bloq <- hessian.bloq + 2 * (contrib_wr - contrib_w0)
      }
    }
  }
  
  out <- objlist(value = obj.bloq, gradient = grad.bloq, hessian = hessian.bloq)
  return(out)
}



## Methods for class objlist ------------------------------------------------

#' Add two lists element by element
#' 
#' @param out1 List of numerics or matrices
#' @param out2 List with the same structure as out1 (there will be no warning when mismatching)
#' @details If out1 has names, out2 is assumed to share these names. Each element of the list out1
#' is inspected. If it has a \code{names} attributed, it is used to do a matching between out1 and out2.
#' The same holds for the attributed \code{dimnames}. In all other cases, the "+" operator is applied
#' the corresponding elements of out1 and out2 as they are.
#' @return List of length of out1. 
#' @aliases sumobjlist
#' @export
#' 
"+.objlist" <- function(out1, out2) {
  
  if (is.null(out1)) return(out2)
  if (is.null(out2)) return(out1)
  
  what <- intersect(c("value", "gradient", "hessian"), c(names(out1), names(out2)))
  
  add_vector <- function(a,b) {
    # add vector b to a by names
    i <- intersect(names(a), names(b))
    a[i] <- a[i] + b[i]
    a}
  add_matrix <- function(a,b) {
    i <- intersect(rownames(a), rownames(b))
    a[i,i] <- a[i,i] + b[i,i]
    a}
  
  gn1 <- names(out1$gradient)
  gn2 <- names(out2$gradient)
  
  one_includes_two <- all(gn2 %in% gn1) 
  two_includes_one <- all(gn1 %in% gn2)
  neither_included <- !(one_includes_two | two_includes_one)
  
  out12 <- lapply(what, function(w) {
    v1 <- out1[[w]]
    v2 <- out2[[w]]
    if (w == "value") 
      return(v1 + v2)
    if (w == "gradient"){
      if (neither_included) return(add_vector(add_vector(setNames(rep(0, length(union(gn1, gn2))), union(gn1, gn2)),v1),v2))
      if (one_includes_two) return(add_vector(v1,v2))
      if (two_includes_one) return(add_vector(v2,v1))
    }
    if (w == "hessian") {
      if (neither_included) return(add_matrix(add_matrix(matrix(0, length(union(gn1,gn2)),length(union(gn1,gn2)),
                                                                dimnames = list(union(gn1,gn2), union(gn1,gn2))
      ),v1),v2))
      if (one_includes_two) return(add_matrix(v1,v2))
      if (two_includes_one) return(add_matrix(v2,v1))
    }
  })
  names(out12) <- what
  
  # Summation of numeric attributes 
  out1.attributes <- attributes(out1)[sapply(attributes(out1), is.numeric)]
  out2.attributes <- attributes(out2)[sapply(attributes(out2), is.numeric)]
  attr.names <- union(names(out1.attributes), names(out2.attributes))
  out12.attributes <- lapply(attr.names, function(n) {
    x1 <- ifelse(is.null(out1.attributes[[n]]), 0, out1.attributes[[n]])
    x2 <- ifelse(is.null(out2.attributes[[n]]), 0, out2.attributes[[n]])
    x1 + x2
  })
  attributes(out12)[attr.names] <- out12.attributes
  
  class(out12) <- "objlist"
  return(out12)
}


#' @export
print.objlist <- function(x, n1 = 20, n2 = 6, ...) {
  n1 <- min(n1,length(x$gradient))
  n2 <- min(n2,length(x$gradient))
  cat("value\n", "==================\n",x$value, "\n")
  cat("gradient[1:",n1,"] (full length = ",length(x$gradient),")\n", "==================\n", sep = "")
  print(x$gradient[1:n1])
  cat("\n")
  cat("hessian[1:",n2,",1:",n2,"]","\n", "==================\n", sep = "")
  print(x$hessian[1:n2,1:n2])
  cat("\n\n")
  cat("attributes\n", "==================\n")
  cat(capture.output(str(attributes(x), max.level = 1)), sep = "\n")
  
}



#' @export
print.objfn <- function(x, ...) {
  
  parameters <- attr(x, "parameters")
  
  cat("Objective function:\n")
  str(args(x))
  cat("\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  
}
