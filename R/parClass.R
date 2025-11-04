## Methods for the class parlist -----------------------------------------------

#' Parameter list
#' 
#' @param x list of lists, as returned by \code{trust}
#' @rdname parlist
#' @export
as.parlist <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  } else {
    class(x) <- c("parlist", "list")
    return(x)
  }
}

#' @export
#' @param object a parlist
#' @rdname parlist
summary.parlist <- function(object, ...) {
  
  x <- object
  
  # Statistics
  m_stat <- stat.parlist(x)
  m_error <- sum(m_stat == "error")
  m_converged <- sum(m_stat == "converged")
  m_notConverged <- sum(m_stat == "notconverged")
  m_sumStatus <- sum(m_error, m_converged, m_notConverged)
  m_total <- length(m_stat)
  
  # Best and worst fit
  m_parframe <- as.parframe(x)
  m_order <- order(m_parframe$value)
  m_bestWorst <- m_parframe[c(m_order[1], tail(m_order, 1)),]
  rownames(m_bestWorst) <- c("best", "worst")
  cat("Results of the best and worst fit\n")
  print(m_bestWorst)
  
  cat("\nStatistics of fit outcome",
      "\nFits aborted:       ", m_error,
      "\nFits not converged: ", m_notConverged,
      "\nFits converged:     ", m_converged,
      "\nFits total:         ", m_sumStatus, " [", m_total, "]", sep = "")
}



#' Gather statistics of a fitlist
#' @param x The fitlist
stat.parlist <- function(x) {
  status <- do.call(rbind, lapply(x, function(fit) {
    if (inherits(fit, "try-error") || any(names(fit) == "error") || any(is.null(fit))) {
      return("error")
    } else {
      if (fit$converged) {
        return("converged")
      } else {
        return("notconverged")
      }
    }
  }))
  
  rownames(status) <- 1:length(status)
  colnames(status) <- "fit status"
  
  return(status)
}


#' Plot a parameter list.
#' 
#' @param x fitlist obtained from mstrust
#' @param ... additional arguments
#' @param path print path of parameters from initials to convergence. For this
#'   option to be TRUE \code{\link{mstrust}} must have had the option
#'   \option{blather}.
#' 
#' @details If path=TRUE:        
#' @author Malenka Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
plot.parlist <- function(x, path = FALSE, ...) {
  
  pl <- x
  
  index <- do.call(rbind, lapply(pl, function(l) l$converged))
  fl <- pl[index]
  if (!path) {
    initPar <- do.call(rbind, lapply(fl, function(l) l$parinit))
    convPar <- do.call(rbind, lapply(fl, function(l) l$argument))
    
    ddata <- data.frame(cbind(matrix(initPar, ncol = 1), matrix(convPar, ncol = 1) ))
    ddata <- cbind(rep(colnames(initPar), each = nrow(initPar)), ddata, 1)
    names(ddata) <- c("parameter","x","y","run")
    
    #plot initial vs converged parameter values
    ggplot(data=ddata)+facet_wrap(~ parameter)+geom_point(aes(x=x,y=y))
  } else {
    if (!any (names(fl[[1]]) == "argpath")){
      stop("No path information in the output of mstrust. Restart mstrust with option blather.")
    }
    parNames <- names(fl[[1]]$parinit)
    
    pathPar <- do.call(rbind, mapply(function(l, idx) {
      mParPath <- as.data.frame(matrix(l$argpath, ncol = 1))
      mParPath <- cbind(rep(parNames,each = nrow(l$argpath), times = 1), rep(1:nrow(l$argpath), length(parNames)), mParPath, as.character(idx))
    }, l = fl, idx = 1:length(fl), SIMPLIFY = FALSE))
    names(pathPar) <- c("parameter", "iteration", "path", "idx")
    ggplot(data=pathPar)+geom_line(aes(x=iteration,y=path,colour=idx))+facet_wrap(~ parameter)
  }
}




#' @export
#' @importFrom data.table as.data.table rbindlist
#' @rdname as.parframe
#' @param sort.by character indicating by which colum the returned parameter frame
#' should be sorted. Defaults to \code{"value"}.
as.parframe.parlist <- function(x, sort.by = "value", ...) {
  m_stat <- stat.parlist(x)
  m_metanames <- c("index", "value", "converged", "iterations")
  m_idx <- which("error" != m_stat)
  m_parframe <- data.frame(index = m_idx, 
                           value = vapply(x[m_idx], function(.x) .x$value, 1.0),
                           converged = vapply(x[m_idx], function(.x) .x$converged, TRUE),
                           iterations = vapply(x[m_idx], function(.x) as.integer(.x$iterations), 1L))
  
  parameters <- lapply(x[m_idx], function(x) data.table::as.data.table(as.list(x$argument)))
  parameters <- data.table::rbindlist(parameters, use.names = TRUE)
  m_parframe <- cbind(m_parframe, parameters)
  
  # Sort by value
  m_parframe <- m_parframe[order(m_parframe[[sort.by]]),]
  
  parframe(m_parframe, parameters = names(parameters), metanames = m_metanames)
}



#' Concatenate parameter lists
#'
#' @description Fitlists carry an fit index which must be held unique on merging
#' multiple fitlists.
#'
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'
#' @rdname parlist
#' @export
#' @export c.parlist
c.parlist <- function(...) {
  m_fits <- lapply(list(...), unclass)
  m_fits <- do.call(c, m_fits)
  m_parlist <- mapply(function(fit, idx) {
    if (is.list(fit)) fit$index <- idx
    return(fit)
  }, fit = m_fits, idx = seq_along(m_fits), SIMPLIFY = FALSE)
  
  return(as.parlist(m_parlist))
}





## Methods for the class parframe ----


#' Coerce object to a parameter frame
#' 
#' @param x object to be coerced
#' @param ... other arguments
#' @return object of class \link{parframe}.
#' @example inst/examples/parlist.R
#' @export
as.parframe <- function(x, ...) {
  UseMethod("as.parframe", x)
}


#' Select a parameter vector from a parameter frame.
#' 
#' @description Obtain a parameter vector from a parameter frame.
#' 
#' @param x A parameter frame, e.g., the output of
#'   \code{\link{as.parframe}}.
#' @param index Integer, the parameter vector with the \code{index}-th lowest
#'   objective value.
#' @param ... not used right now
#'   
#' @details With this command, additional information included in the parameter
#'   frame as the objective value and the convergence state are removed and a
#'   parameter vector is returned. This parameter vector can be used to e.g.,
#'   evaluate an objective function.
#'   
#'   On selection, the parameters in the parameter frame are ordered such, that
#'   the parameter vector with the lowest objective value is at \option{index}
#'   1. Thus, the parameter vector with the \option{index}-th lowest objective
#'   value is easily obtained.
#'   
#' @return The parameter vector with the \option{index}-th lowest objective
#'   value.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#'   
#' @export
as.parvec.parframe <- function(x, index = 1, ...) {
  parframe <- x
  m_order <- 1:nrow(x)
  metanames <- attr(parframe, "metanames")
  if ("value" %in% metanames) m_order <- order(parframe$value)
  best <- as.parvec(unlist(as.data.frame(parframe)[m_order[index], attr(parframe, "parameters"), drop = FALSE]))
  if ("converged" %in% metanames && !parframe[m_order[index],]$converged) {
    warning("Parameter vector of an unconverged fit is selected.", call. = FALSE)
  }
  return(best)
}




#' @export
#' @rdname plotPars
plotPars.parframe <- function(x, tol = 1, ...){
  
  if (!missing(...)) x <- subset(x, ...)
  
  jumps <- stepDetect(x$value, tol)
  jump.index <- approx(jumps, jumps, xout = 1:length(x$value), method = "constant", rule = 2)$y
  
  #values <- round(x$value/tol)
  #unique.values <- unique(values)
  #jumps <- which(!duplicated(values))
  #jump.index <- jumps[match(values, unique.values)]
  x$index <- as.factor(jump.index)
  
  myparframe <- x
  parNames <- attr(myparframe,"parameters")
  parOut <- wide2long.data.frame(out = ((myparframe[, c("index", "value", parNames)])) , keep = 1:2)
  names(parOut) <- c("index", "value", "name", "parvalue")
  plot <- ggplot2::ggplot(parOut, aes(x = name, y = parvalue, color = index)) + geom_boxplot(outlier.alpha = 0) + theme_dMod() + scale_color_dMod() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  attr(plot, "data") <- parOut
  
  return(plot)
  
}


#' @export
#' @rdname plotValues
plotValues.parframe <- function(x, tol = 1, ...) {
  
  if (!missing(...)) x <- subset(x, ...)
  
  jumps <- stepDetect(x$value, tol)
  y.range <- c(min(x$value), max(max(x$value), min(x$value) + tol))
  y.jumps <- seq(y.range[2], y.range[1], length.out = length(jumps))
  
  
  pars <- x
  pars <- pars[order(pars$value),]
  pars[["index"]] <-  1:nrow(pars)
  
  
  
  P <- ggplot2::ggplot(pars, aes(x = index, y = value, pch = converged, color = iterations)) + 
    geom_vline(xintercept = jumps, lty = 2) +
    geom_point() + 
    annotate("text", x = jumps + 1, y = y.jumps, label = jumps, hjust = 0, color = "firebrick", size = 3) +
    xlab("index") + ylab("value") + 
    scale_color_gradient(low = "dodgerblue", high = "orange") +
    coord_cartesian(ylim = y.range) +
    theme_dMod()
  
  attr(P, "data") <- pars
  attr(P, "jumps") <- jumps
  
  return(P)
  
}



#' @export
#' @rdname plotProfile
plotProfile.parframe <- function(profs, ..., maxvalue = 5, parlist = NULL) {
  
  if("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])
    
    
    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }
    
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue
      
      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)
      
      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }
      
      return(sub)
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))

  
  data.zero <- subset(data, is.zero)
  
  threshold <- c(1, 2.7, 3.84)
  
  data <- droplevels.data.frame(subset(data, ...))

  
  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") + 
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")
  
  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){ 
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]  
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)
    
    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)
  
}


#' @export
#' @rdname plotProfile
plotProfile.list <- function(profs, ..., maxvalue = 5, parlist = NULL) {
  
  if("parframe" %in% class(profs)) 
    arglist <- list(profs)
  else
    arglist <- as.list(profs)
  
  
  if (is.null(names(arglist))) {
    profnames <- 1:length(arglist)
  } else {
    profnames <- names(arglist)
  }
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- as.data.frame(arglist[[i]])
    obj.attributes <- attr(arglist[[i]], "obj.attributes")
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])
    
    
    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }
    
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue
      
      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = profnames[i], mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)
      
      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = profnames[i], mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }
      
      return(sub)
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  data <- droplevels(subset(data, ...))

  
  data.zero <- subset(data, is.zero)
  
  threshold <- c(1, 2.7, 3.84)
  
  data <- droplevels.data.frame(subset(data, ...))
  
  data$y <- as.numeric(data$y)
  data$x <- as.numeric(data$x)
  
  
  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") + 
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), limits = c(NA, maxvalue)) +
    xlab("parameter value")
  
  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){ 
        origin <- which.min(arglist[[i]][["constraint"]])
        zerovalue <- arglist[[i]][origin, 1]  
      })))
      values <- parlist[, "value", drop = TRUE]
      parlist <- parlist[,!(colnames(parlist) %in% c("index", "value", "converged", "iterations"))]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)
    
    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta), color = "black", inherit.aes = FALSE)
  }
  attr(p, "data") <- data
  return(p)
  
}



#' @export
#' @rdname parframe
is.parframe <- function(x) {
  "parframe" %in% class(x)
}

#' @export
#' @param i row index in any format
#' @param j column index in any format
#' @param drop logical. If TRUE the result is coerced to the lowest possible dimension
#' @rdname parframe
"[.parframe" <- function(x, i = NULL, j = NULL, drop = FALSE){
  
  metanames <- attr(x, "metanames")
  obj.attributes <- attr(x, "obj.attributes")
  parameters <- attr(x, "parameters")
  
  out <- as.data.frame(x)
  #out <- as.data.frame(unclass(x))
  if (!is.null(i)) out <- out[i, ]
  if (!is.null(j)) out <- out[, j, drop = drop]
  
  if (drop) return(out)
  
  metanames <- intersect(metanames, colnames(out))
  obj.attributes <- intersect(obj.attributes, colnames(out))
  parameters <- intersect(parameters, colnames(out))
  
  parframe(out, parameters = parameters, metanames = metanames, obj.attributes = obj.attributes)
  
}


#' @export
#' @param ... additional arguments
#' @rdname parframe
subset.parframe <- function(x, ...) {
  
  x[with(as.list(x), ...), ]
  
}

#' Extract those lines of a parameter frame with unique elements in the value column
#' @param x parameter frame
#' @param incomparables not used. Argument exists for compatibility with S3 generic.
#' @param tol tolerance to decide when values are assumed to be equal, see \code{\link{plotValues}()}.
#' @param ... additional arguments being passed to \code{\link{plotValues}()}, e.g. for subsetting.
#' @return A subset of the parameter frame \code{x}.
#' @export
unique.parframe <- function(x, incomparables = FALSE, tol = 1, ...) {
  
  
  jumps <- attr(plotValues(x = x, tol = tol, ...), "jumps")
  x[jumps, ]
  
  
}



## Methods for the class parvec ------------------------------------------------

#' Dispatch as.parvec.
#'
#' @export
#' @rdname parvec
as.parvec <- function(x, ...) {
  UseMethod("as.parvec", x)
}


#' Parameter vector (augmented with first- and second-order derivatives)
#'
#' Creates an object of class `parvec`, storing parameter values together with
#' their first (`deriv`) and second (`deriv2`) derivatives.
#'
#' @param x Numeric vector of parameter values.
#' @param names Optional character vector of parameter names.
#' @param deriv Either a logical flag, a Jacobian matrix, or `NULL`.
#'   * If `FALSE`, no Jacobian is stored (attribute set to `NULL`).
#'   * If `NULL`, an identity matrix is assumed.
#'   * If a matrix, it is used (and missing parameters are filled with identity rows/cols).
#' @param deriv2 Either a logical flag, a 3D array, or `NULL`.
#'   * If `FALSE`, no Hessian is stored (attribute set to `NULL`).
#'   * If `NULL`, a zero tensor is assumed.
#'   * If an array, it is used (and missing parameters are filled with zeros).
#'
#' @return An object of class `"parvec"`, i.e. a named numeric vector with attributes:
#' \itemize{
#'   \item \code{attr(x, "deriv")} — Jacobian matrix (or `NULL`)
#'   \item \code{attr(x, "deriv2")} — Hessian tensor (or `NULL`)
#' }
#' @export
as.parvec.numeric <- function(x, names = NULL, deriv = NULL, deriv2 = NULL, ...) {
  
  # --- Base vector and naming ---
  p <- as.numeric(x)
  if (is.null(names)) names(p) <- names(x) else names(p) <- names
  n <- length(p)
  pnames <- names(p)
  
  # === Handle Jacobian (first derivatives) ===
  if (isFALSE(deriv)) {
    full_deriv <- NULL
  } else {
    if (is.null(deriv)) deriv <- attr(x, "deriv", exact = TRUE)
    if (is.null(deriv)) {
      full_deriv <- diag(n)
      dimnames(full_deriv) <- list(pnames, pnames)
    } else {
      full_deriv <- matrix(0, n, n, dimnames = list(pnames, pnames))
      rn <- rownames(deriv); cn <- colnames(deriv)
      if (!is.null(rn) && !is.null(cn)) {
        row.idx <- match(rn, pnames, nomatch = 0)
        col.idx <- match(cn, pnames, nomatch = 0)
        valid.r <- row.idx > 0; valid.c <- col.idx > 0
        if (any(valid.r) && any(valid.c))
          full_deriv[row.idx[valid.r], col.idx[valid.c]] <- deriv[valid.r, valid.c, drop = FALSE]
      }
      missing <- setdiff(pnames, rn)
      if (length(missing) > 0)
        full_deriv[missing, missing] <- diag(length(missing))
    }
  }
  
  # === Handle Hessian (second derivatives) ===
  if (isFALSE(deriv2)) {
    full_deriv2 <- NULL
  } else {
    if (is.null(deriv2)) deriv2 <- attr(x, "deriv2", exact = TRUE)
    if (is.null(deriv2)) {
      full_deriv2 <- array(0, dim = c(n, n, n),
                           dimnames = list(pnames, pnames, pnames))
    } else {
      full_deriv2 <- array(0, dim = c(n, n, n),
                           dimnames = list(pnames, pnames, pnames))
      dn <- dimnames(deriv2)
      if (!is.null(dn[[1]]) && !is.null(dn[[2]]) && !is.null(dn[[3]])) {
        idx1 <- match(dn[[1]], pnames, nomatch = 0)
        idx2 <- match(dn[[2]], pnames, nomatch = 0)
        idx3 <- match(dn[[3]], pnames, nomatch = 0)
        valid1 <- idx1 > 0; valid2 <- idx2 > 0; valid3 <- idx3 > 0
        if (any(valid1) && any(valid2) && any(valid3))
          full_deriv2[idx1[valid1], idx2[valid2], idx3[valid3]] <-
          deriv2[valid1, valid2, valid3, drop = FALSE]
      }
    }
  }
  
  # === Assemble final object ===
  attr(p, "deriv")  <- full_deriv
  attr(p, "deriv2") <- full_deriv2
  class(p) <- c("parvec", "numeric")
  p
}


#' Pretty printing for parvec objects (supports second derivatives)
#'
#' Prints a parameter vector along with information about
#' its attached derivative matrices (`deriv`, `deriv2`).
#'
#' @param x parvec object
#' @export
print.parvec <- function(x) {
  
  # --- Extract parameter values ---
  par <- unclass(x)
  nms <- names(par)
  n_width <- max(nchar(nms))
  
  # --- Header ---
  cat("Parameter vector:\n")
  
  # --- Print parameters neatly aligned ---
  for (i in seq_along(par)) {
    val <- formatC(par[i], digits = 6, format = "g")
    if (par[i] >= 0) val <- paste0(" ", val)
    cat(sprintf("  %s : %s\n", format(nms[i], width = n_width, justify = "right"), val))
  }
  
  # --- Print derivative info ---
  deriv  <- attr(x, "deriv")
  deriv2 <- attr(x, "deriv2")
  
  cat("\nAttributes:\n")
  if (!is.null(deriv)) {
    d <- dim(deriv)
    cat(sprintf("  deriv  : %d × %d matrix\n", d[1], d[2]))
  } else {
    cat("  deriv  : <none>\n")
  }
  
  if (!is.null(deriv2)) {
    d2 <- dim(deriv2)
    cat(sprintf("  deriv2 : %d × %d × %d array\n", d2[1], d2[2], d2[3]))
  } else {
    cat("  deriv2 : <none>\n")
  }
  
  invisible(x)
}

#' Subset method for parvec objects (supports second derivatives)
#'
#' This method subsets a parameter vector (`parvec`) while keeping its
#' first (`deriv`) and second (`deriv2`) derivative attributes consistent.
#'
#' If `drop = TRUE`, all-zero columns are removed from the Jacobian,
#' and the corresponding dimensions are also removed from the Hessian.
#'
#' @param x parvec object
#' @param ... indices or names to subset
#' @param drop logical; if TRUE, remove empty Jacobian columns
#'   and matching Hessian dimensions
#' @export
"[.parvec" <- function(x, ..., drop = FALSE) {
  
  # --- Extract numeric parameter values ---
  out <- unclass(x)[...]
  nms <- names(out)
  
  # --- Subset first derivative (Jacobian) ---
  deriv <- submatrix(attr(x, "deriv"), rows = ..., cols = TRUE)
  
  # --- Subset second derivative (Hessian tensor) ---
  deriv2 <- attr(x, "deriv2")
  if (!is.null(deriv2)) {
    if (!is.null(names(x))) {
      sel <- nms
      deriv2 <- deriv2[sel, sel, sel, drop = FALSE]
    } else {
      idx <- seq_along(out)
      deriv2 <- deriv2[idx, idx, idx, drop = FALSE]
    }
  }
  
  # --- Optionally drop empty Jacobian columns and Hessian dims ---
  if (drop) {
    empty.cols <- apply(deriv, 2, function(v) all(v == 0))
    keep.cols <- !empty.cols
    deriv <- submatrix(deriv, cols = keep.cols)
    
    if (!is.null(deriv2)) {
      keep.names <- colnames(deriv)
      deriv2 <- deriv2[keep.names, keep.names, keep.names, drop = FALSE]
    }
  }
  
  # --- Reconstruct parvec object with consistent attributes ---
  as.parvec(out, deriv = deriv, deriv2 = deriv2)
}



#' Concatenate parvec objects (supports second derivatives)
#'
#' Combines multiple parameter vectors (`parvec`) into a single object,
#' preserving their first (`deriv`) and second (`deriv2`) derivatives.
#'
#' The Jacobians are combined block-diagonally, and the Hessians are
#' combined as block-diagonal tensors (no cross-terms between distinct blocks).
#'
#' @param ... parvec objects to concatenate
#' @export
c.parvec <- function(...) {
  parlist <- list(...)
  nms <- unlist(lapply(parlist, names))
  vals <- unlist(lapply(parlist, as.numeric))
  
  if (any(duplicated(nms)))
    stop("Duplicated parameter names detected — cannot concatenate parvecs.")
  
  # --- Combine first derivatives (Jacobian) ---
  derivs <- lapply(parlist, function(p) attr(p, "deriv"))
  deriv <- Reduce(combine, derivs)
  rownames(deriv) <- nms
  
  # --- Combine second derivatives (Hessian / deriv2) ---
  deriv2s <- lapply(parlist, function(p) attr(p, "deriv2"))
  if (any(!vapply(deriv2s, is.null, TRUE))) {
    # Filter out NULL Hessians, create empty zero blocks where missing
    for (i in seq_along(deriv2s)) {
      if (is.null(deriv2s[[i]])) {
        n <- length(parlist[[i]])
        nm <- names(parlist[[i]])
        deriv2s[[i]] <- array(0, dim = c(n, n, n),
                              dimnames = list(nm, nm, nm))
      }
    }
    
    # Construct block-diagonal 3D tensor
    all_names <- nms
    n_total <- length(all_names)
    deriv2 <- array(0, dim = c(n_total, n_total, n_total),
                    dimnames = list(all_names, all_names, all_names))
    
    offset <- 0
    for (i in seq_along(deriv2s)) {
      n_i <- length(parlist[[i]])
      nm_i <- names(parlist[[i]])
      deriv2[nm_i, nm_i, nm_i] <- deriv2s[[i]][nm_i, nm_i, nm_i]
      offset <- offset + n_i
    }
  } else {
    deriv2 <- NULL
  }
  
  # --- Assemble final parvec ---
  as.parvec(vals, names = nms, deriv = deriv, deriv2 = deriv2)
}




## Methods for the class parfn--------------------------------------------------

#' Pretty printing parameter transformations
#' 
#' @param x prediction function
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
print.parfn <- function(x, ...) {
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Parameter transformation:\n")
  str(args(x))
  cat("\n")
  cat("... conditions:", paste0(conditions, collapse = ", "), "\n")
  cat("... parameters:", paste0(parameters, collapse = ", "), "\n")
  
}

#' @export
summary.parfn <- function(object, ...) {
  
  x <- object
  
  conditions <- attr(x, "conditions")
  parameters <- attr(x, "parameters")
  mappings <- attr(x, "mappings")
  
  cat("Details:\n")
  if (!inherits(x, "composed")) {
    
    
    output <- lapply(1:length(mappings), function(C) {
      
      list(
        equations = attr(mappings[[C]], "equations"),
        parameters = attr(mappings[[C]], "parameters")
      )
      
    })
    names(output) <- conditions
    
    #print(output, ...)
    output
    
  } else {
    
    cat("\nObject is composed. See original objects for more details.\n")
    
  }
}


