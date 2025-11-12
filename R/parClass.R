## Methods for the class parlist -----------------------------------------------

#' Parameter list
#' 
#' @param x list of lists, as returned by `trust`
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
#'   option to be TRUE [mstrust()] must have had the option
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
#' should be sorted. Defaults to `"value"`.
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
#' @return object of class [parframe].
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
#'   [as.parframe()].
#' @param index Integer, the parameter vector with the `index`-th lowest
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
#' @param tol tolerance to decide when values are assumed to be equal, see [plotValues()].
#' @param ... additional arguments being passed to [plotValues()], e.g. for subsetting.
#' @return A subset of the parameter frame `x`.
#' @export
unique.parframe <- function(x, incomparables = FALSE, tol = 1, ...) {
  
  
  jumps <- attr(plotValues(x = x, tol = tol, ...), "jumps")
  x[jumps, ]
  
  
}



## Methods for the class parvec ------------------------------------------------

#' Dispatch as.parvec.
#'
#' Creates an object of class `parvec`, optionally carrying first- and
#' second-order derivatives. No automatic identity or zero derivatives are
#' created.
#'
#' @param x Numeric vector of parameter values.
#' @param names Optional character vector of parameter names.
#' @param deriv Optional Jacobian matrix (first derivatives), `FALSE` to drop,
#'   or `NULL` to inherit if present attribute 'deriv'.
#' @param deriv2 Optional Hessian array (second derivatives), `FALSE` to drop,
#'   or `NULL` to inherit if present atribute 'deriv2'.
#'
#' @return
#' A numeric vector of class `"parvec"` with optional attributes:
#' \itemize{
#'   \item `attr(x, "deriv")`  — Jacobian matrix, or `NULL`
#'   \item `attr(x, "deriv2")` — Hessian tensor, or `NULL`
#' }
#'
#' @details
#' If `deriv = FALSE` or `deriv2 = FALSE`, the respective attributes are
#' forcibly removed, even if they exist in the input `x`.
#' If `deriv = NULL` or `deriv2 = NULL`, existing attributes are inherited
#' (if present), but nothing is created automatically.
#'
#' @export
#' @rdname parvec
as.parvec <- function(x, ...) {
  UseMethod("as.parvec", x)
}


#' @export
#' @rdname parvec
as.parvec.numeric <- function(x, names = NULL, deriv = NULL, deriv2 = NULL, ...) {
  
  # --- Basic setup ---
  p <- as.numeric(x)
  if (is.null(names)) names(p) <- names(x) else names(p) <- names
  pnames <- names(p)
  
  # --- First derivatives ---
  if (isFALSE(deriv)) {
    full_deriv <- NULL
  } else if (is.matrix(deriv)) {
    full_deriv <- deriv
  } else { # deriv == NULL
    full_deriv <- attr(x, "deriv", exact = TRUE)
  }
  
  # --- Second derivatives ---
  if (isFALSE(deriv2)) {
    full_deriv2 <- NULL
  } else if (is.array(deriv2)) {
    full_deriv2 <- deriv2
  } else { # deriv2 == NULL
    full_deriv2 <- attr(x, "deriv2", exact = TRUE)
  }
  
  # --- Assemble object ---
  attr(p, "deriv")  <- full_deriv
  attr(p, "deriv2") <- full_deriv2
  class(p) <- c("parvec", "numeric")
  p
}





#' Pretty printing for parvec objects
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


#' Subset method for parvec objects
#'
#' Safely subsets a `parvec` object while keeping its Jacobian (`deriv`)
#' and Hessian (`deriv2`) attributes consistent.
#'
#' This version is robust against missing derivative entries:
#' if parameters being subsetted are not present in the Jacobian or Hessian,
#' zero rows (and matching dimensions) are added automatically.
#'
#' @param x A `parvec` object.
#' @param ... Parameter indices or names to subset.
#' @param drop Logical; if `TRUE`, drop zero Jacobian columns and
#'   matching Hessian dimensions.
#'
#' @return A `parvec` object with updated numeric values and derivative
#'   attributes.
#' @export
"[.parvec" <- function(x, ..., drop = FALSE) {
  # --- Extract numeric parameter values ---
  out <- unclass(x)[...]
  nms <- names(out)
  
  # --- Subset Jacobian (first derivatives) ---
  deriv <- attr(x, "deriv")
  if (!is.null(deriv)) {
    available <- intersect(nms, rownames(deriv))
    missing <- setdiff(nms, rownames(deriv))
    
    # Subset existing rows
    deriv_sub <- deriv[available, , drop = FALSE]
    
    # Add zero rows for missing parameters
    if (length(missing)) {
      add_rows <- matrix(0, nrow = length(missing), ncol = ncol(deriv),
                         dimnames = list(missing, colnames(deriv)))
      deriv_sub <- rbind(deriv_sub, add_rows)
    }
    
    # Reorder rows to match the subset order
    deriv <- deriv_sub[nms, , drop = FALSE]
  }
  
  # --- Subset Hessian (second derivatives) ---
  deriv2 <- attr(x, "deriv2")
  if (!is.null(deriv2)) {
    basis_names <- dimnames(deriv2)[[2]]
    available1 <- intersect(nms, dimnames(deriv2)[[1]])
    missing1 <- setdiff(nms, dimnames(deriv2)[[1]])
    
    # Subset existing first dimension
    deriv2_sub <- deriv2[available1, , , drop = FALSE]
    
    # Add zero slices for missing parameters (first dimension)
    if (length(missing1)) {
      add_slices <- array(0, dim = c(length(missing1), dim(deriv2)[2], dim(deriv2)[3]),
                          dimnames = list(missing1, dimnames(deriv2)[[2]], dimnames(deriv2)[[3]]))
      deriv2_sub <- abind::abind(deriv2_sub, add_slices, along = 1)
    }
    
    # Reorder first dimension to match the subset order
    deriv2 <- deriv2_sub[nms, , , drop = FALSE]
  }
  
  # --- Drop zero columns and matching Hessian dimensions if requested ---
  if (drop && !is.null(deriv)) {
    keep.cols <- colSums(abs(deriv)) > 0
    deriv <- deriv[, keep.cols, drop = FALSE]
    if (!is.null(deriv2)) {
      deriv2 <- deriv2[, keep.cols, keep.cols, drop = FALSE]
    }
  }
  
  # --- Rebuild the parvec object ---
  as.parvec(out, deriv = deriv, deriv2 = deriv2)
}


#' Concatenate parvec objects
#'
#' Combines multiple `parvec` objects into one, preserving first (`deriv`)
#' and second (`deriv2`) derivatives when available.
#'
#' @param ... `parvec` objects to concatenate.
#' @return A combined `parvec` object.
#' @export
c.parvec <- function(...) {
  parlist <- list(...)
  stopifnot(length(parlist) > 0)
  
  # Concatenate numeric values and names
  nms  <- unlist(lapply(parlist, names), use.names = FALSE)
  vals <- unlist(lapply(parlist, as.numeric), use.names = FALSE)
  if (anyDuplicated(nms))
    stop("Duplicated parameter names detected.")
  
  # Check available derivative attributes
  derivs  <- lapply(parlist, function(p) attr(p, "deriv"))
  deriv2s <- lapply(parlist, function(p) attr(p, "deriv2"))
  has_deriv  <- any(vapply(derivs,  function(d) !is.null(d),  TRUE))
  has_deriv2 <- any(vapply(deriv2s, function(d) !is.null(d), TRUE))
  
  # Determine derivative basis only if consistent
  basis <- NULL
  if (has_deriv) {
    basis <- colnames(derivs[[which(!vapply(derivs, is.null, TRUE))[1]]])
  } else if (has_deriv2) {
    basis <- dimnames(deriv2s[[which(!vapply(deriv2s, is.null, TRUE))[1]]])[[2]]
  }
  
  nb <- if (!is.null(basis)) length(basis) else 0
  nt <- length(nms)
  
  # Preallocate only if consistent across all parvecs
  deriv  <- if (has_deriv)  matrix(0, nt, nb, dimnames = list(nms, basis)) else NULL
  deriv2 <- NULL
  
  if (has_deriv2) {
    # Ensure all non-null deriv2 blocks have same basis dims
    valid_dims <- unique(Filter(Negate(is.null),
                                lapply(deriv2s, function(h) dim(h)[2:3])))
    if (length(valid_dims) == 1L) {
      deriv2 <- array(0, dim = c(nt, nb, nb),
                      dimnames = list(nms, basis, basis))
    } else {
      warning("Inconsistent deriv2 dimensions across parvecs – dropping second derivatives.")
      has_deriv2 <- FALSE
    }
  }
  
  # Fill derivative blocks
  row_start <- 1
  for (i in seq_along(parlist)) {
    np <- length(parlist[[i]])
    
    d1 <- derivs[[i]]
    if (!is.null(deriv) && !is.null(d1))
      deriv[row_start:(row_start + np - 1), ] <- d1
    
    d2 <- deriv2s[[i]]
    if (!is.null(deriv2) && !is.null(d2) &&
        all(dim(d2)[2:3] == dim(deriv2)[2:3]))
      deriv2[row_start:(row_start + np - 1), , ] <- d2
    
    row_start <- row_start + np
  }
  
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


