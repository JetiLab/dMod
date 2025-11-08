

#' Model prediction function for ODE models. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function \code{x(times, pars, deriv = TRUE)} returning ODE output and sensitivities.
#' @param odemodel object of class 'odemodel' or 'odemodel++', see \link{odemodel}
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings. Not (yet) implemented for boost::odeint::rosenbrock4
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See \link[deSolve]{events}.
#' ATTENTION: Sensitivities for event states will only be correctly computed if defined within
#' \code{\link{odemodel}()}. Specify events within \code{Xs()} only for forward simulation.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See \link[stats]{approxfun} for possible arguments.
#' @return Object of class \link{prdfn}. If the function is called with parameters that
#' result from a parameter transformation (see \link{P}), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @export
Xs <- function(odemodel, ...) {
  UseMethod("Xs", odemodel)
}

#' @export
#' @import deSolve einsum
Xs.deSolve <- function(odemodel, 
                       forcings=NULL, 
                       events=NULL, 
                       names = NULL, 
                       condition = NULL, 
                       optionsOde=list(method = "lsoda"), 
                       optionsSens=list(method = "lsodes"), 
                       fcontrol = NULL) {
  
  func <- odemodel$func
  extended <- odemodel$extended
  if (is.null(extended)){
    warning("Element 'extended' empty. ODE model does not contain sensitivities.")
    deriv = FALSE
  }
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  if (!is.null(attr(func, "events")) & !is.null(myevents))
    warning("Events already defined in odemodel. Additional events in Xs() will be ignored. Events need to be defined in either odemodel() or Xs().")
  if (is.null(attr(func, "events")) & !is.null(myevents))
    message("Events should be definend in odemodel(). If defined in Xs(), events will be applied, but sensitivities will not be reset accordingly.")
  
  
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  forcnames <- attr(func, "forcings")
  
  # Variable and parameter names of sensitivities
  sensvar <- attr(extended, "variables")[!attr(extended, "variables")%in%variables]
  senssplit <- strsplit(sensvar, ".", fixed=TRUE)
  senssplit.1 <- unlist(lapply(senssplit, function(v) v[1]))
  senssplit.2 <- unlist(lapply(senssplit, function(v) paste(v[-1], collapse = ".")))
  svariables <- intersect(senssplit.2, variables)
  sparameters <- setdiff(senssplit.2, variables)
  
  # Initial values for sensitivities
  yiniSens <- as.numeric(senssplit.1 == senssplit.2)
  names(yiniSens) <- sensvar
  
  # Names for deriv output
  sensGrid <- expand.grid(variables, c(svariables, sparameters), stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  # Only a subset of all variables/forcings is returned
  if (is.null(names)) names <- c(variables, forcnames)
  
  # Update sensNames when names are set
  select <- sensGrid[, 1] %in% names
  sensNames <- paste(sensGrid[,1][select], sensGrid[,2][select], sep = ".")  
  
  
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    fcontrol = myfcontrol
  )
  
  P2X <- function(times, pars, deriv=TRUE, deriv2 = FALSE, env = parent.frame()){
    
    if (deriv2) {
      stop(
        "Second-order sensitivities are not available with the 'deSolve' solver.\n",
        "Consider using solver = 'boost' in odemodel()."
      )
    }
    
    yini <- unclass(pars)[variables]
    mypars <- unclass(pars)[parameters]
    
    forcings <- controls$forcings
    events <- controls$events
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    fcontrol <- controls$fcontrol
    names <- controls$names
    
    # Add event time points (required by integrator) 
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))
    
    # Sort event time points
    if (!is.null(events)) events <- events[order(events$time),]

    
    dX <- NULL
    mysensitivities <- NULL
    if (!deriv) {
      
      # Evaluate model without sensitivities
      # loadDLL(func)
      if (!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
      out <- suppressWarnings(do.call(odeC, c(list(y = unclass(yini), times = times, func = func, parms = mypars, forcings = forc, events = list(data = events), fcontrol = fcontrol), optionsOde)))
      out <- submatrix(out, cols = c("time", names))
      
      
    } else {
      
      if (!is.null(forcings)) forc <- setForcings(extended, forcings) else forc <- NULL
      
      outSens <- suppressWarnings(do.call(odeC, c(list(y = c(unclass(yini), yiniSens), times = times, func = extended, parms = mypars, 
                                      forcings = forc, fcontrol = fcontrol,
                                      events = list(data = events)), optionsSens)))
      
      out <- submatrix(outSens, cols = c("time", names))
      mysensitivities <- reshapeSens(submatrix(outSens, cols = !colnames(outSens) %in% c(variables, forcnames, "time")), variables, c(svariables, sparameters))
      
      # --- Apply parameter transformation to sensitivities (chain rule) ---
      variables <- intersect(variables, names)
      dP <- attr(pars, "deriv")
      
      if (!is.null(dP)) {
        dPsub <- submatrix(dP, rows = c(svariables, sparameters))
        dX <- einsum::einsum("aik,kj->aij", mysensitivities, dPsub)
        dimnames(dX) <- list(NULL, dimnames(mysensitivities)[[2]], colnames(dPsub))
      } else {
        dX <- mysensitivities
      }
      
    }
    
    return(prdframe(out, deriv = dX, sensitivities = mysensitivities, parameters = pars))
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
}

#' Reshape ODE sensitivities from wide matrix (var.par) to 3D array [time, variable, parameter]
#'
#' @description
#' Converts a flat sensitivity matrix (as produced by deSolve-based ODE integrations)
#' with column names of the form "variable.parameter" into a structured 3D array
#' with dimensions [n_time, n_variables, n_parameters].
#'
#' @param sensMatrix A numeric matrix of sensitivities with column names
#'   formatted as "variable.parameter". Rows correspond to time points.
#' @param variables Character vector of state variable names.
#' @param parameters Character vector of parameter names.
#'
#' @return A numeric 3D array with dimensions [time, variable, parameter],
#'   and proper `dimnames` for variables and parameters.
#'
#' @examples
#' # Example column names: "A.k1", "A.k2", "B.k1", "B.k2"
#' mysens <- matrix(runif(40), nrow = 10)
#' colnames(mysens) <- c("A.k1", "A.k2", "B.k1", "B.k2")
#' reshapeSens(mysens, variables = c("A", "B"), parameters = c("k1", "k2"))
#'
#' @keywords internal
reshapeSens <- function(sensMatrix, variables, parameters) {
  n_times <- nrow(sensMatrix)
  n_vars <- length(variables)
  n_pars <- length(parameters)
  
  # Expected column order: var1.par1, var2.par1, ..., varN.par1, var1.par2, var2.par2, ...
  expected_cols <- as.vector(outer(variables, parameters, paste, sep = "."))
  
  # Reorder columns if necessary
  sensMatrix_ordered <- sensMatrix[, expected_cols, drop = FALSE]
  
  # Convert directly to array - array() fills column-wise
  # This matches the data structure perfectly: all vars for par1, then all vars for par2, etc.
  sensArray <- array(
    data = as.matrix(sensMatrix_ordered),
    dim = c(n_times, n_vars, n_pars),
    dimnames = list(NULL, variables, parameters)
  )
  
  return(sensArray)
}

#' @export
#' @importFrom einsum einsum
Xs.boost <- function(odemodel, forcings = NULL, events = NULL, names = NULL, condition = NULL, 
                     optionsOde = list(), optionsSens = list(), optionsSens2 = list()) {
  
  if (!is.null(forcings)) {
    stop("Forcings are not yet supported for boost solver")
  }
  
  if (!is.null(events)) {
    stop("Events should be passed to odemodel()")
  }
  
  optionsDefault  <- list(atol = 1e-6, rtol = 1e-6, maxattemps = 500, maxsteps = 1e6, roottol = 1e-8, maxroot = 1)
  
  ## --- Warn about unknown options
  warn_unknown <- function(user, defaults, label) {
    bad <- setdiff(names(user), names(defaults))
    if (length(bad) > 0)
      warning(sprintf("%s: Ignoring unknown option(s): %s", label, paste(bad, collapse=", ")))
  }
  warn_unknown(optionsOde,  optionsDefault, "optionsOde")
  warn_unknown(optionsSens, optionsDefault, "optionsSens")
  warn_unknown(optionsSens2, optionsDefault, "optionsSens2")
  
  ## --- Merge user-supplied options with defaults
  optionsOde   <- modifyList(optionsDefault, optionsOde)
  optionsSens  <- modifyList(optionsDefault, optionsSens)
  optionsSens2 <- modifyList(optionsDefault, optionsSens2)
  
  func <- odemodel$func
  extended <- odemodel$extended
  extended2 <- odemodel$extended2
  
  if (is.null(extended) || is.null(extended2)) {
    warning(
      sprintf(
        "ODE model does not contain %s-order sensitivities.",
        paste(c("first", "second")[c(is.null(extended), is.null(extended2))],
              collapse = " and ")
      ),
      call. = FALSE
    )
  }
  
  # Variable and parameter names
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  
  dimnames <- attr(func, "dim_names")
  dimnames_sens <- attr(extended, "dim_names")
  
  # Only a subset of all variables is returned
  if (is.null(names)) names <- variables
  
  # Controls to be modified from outside
  controls <- list(
    names = names,
    optionsOde = optionsOde,
    optionsSens = optionsSens,
    optionsSens2 = optionsSens2,
    dimnames = dimnames,
    dimnames_sens = dimnames_sens
  )
  
  P2X <- function(times, pars, deriv=TRUE, deriv2=FALSE, env = parent.frame()) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
      deriv <- TRUE
    }
    
    paramnames <- c(variables, parameters)
    # check for missing parameters
    missing <- setdiff(paramnames, names(pars))
    if (length(missing) > 0) stop(sprintf("Missing parameters: %s", paste(missing, collapse = ", ")))
    
    mypars <- unclass(pars)[paramnames]
    
    names <- controls$names
    optionsOde <- controls$optionsOde
    optionsSens <- controls$optionsSens
    optionsSens2 <- controls$optionsSens2
    
    dX <- NULL
    mysensitivities <- NULL
    dX2 <- NULL
    mysensitivities2 <- NULL
    
    if (!deriv) {
      
      # Evaluate model without sensitivities
      out <- suppressWarnings(
        .Call(paste0("solve_", as.character(func)),
              as.numeric(times),
              as.numeric(mypars),
              as.numeric(optionsOde$atol),
              as.numeric(optionsOde$rtol),
              as.integer(optionsOde$maxattemps),
              as.integer(optionsOde$maxsteps),
              as.numeric(optionsOde$roottol),
              as.integer(optionsOde$maxroot))
      )
      
      colnames(out$variable) <- dimnames$variable
      
      out <- cbind(out$time, submatrix(out$variable, cols = names))
      colnames(out)[1] <- "time"
      
    } else if (deriv & !deriv2) {
      
      outSens <- suppressWarnings(
        .Call(paste0("solve_", as.character(extended)),
              as.numeric(times),
              as.numeric(mypars),
              as.numeric(optionsSens$atol),
              as.numeric(optionsSens$rtol),
              as.integer(optionsSens$maxattemps),
              as.integer(optionsSens$maxsteps),
              as.numeric(optionsSens$roottol),
              as.integer(optionsSens$maxroot))
      )
      
      colnames(outSens$variable) <- controls$dimnames_sens$variable
      dimnames(outSens$sens1) <- list(NULL, controls$dimnames_sens$variable, controls$dimnames_sens$sens)
      
      out <- cbind(outSens$time, submatrix(outSens$variable, cols = names))
      colnames(out)[1] <- "time"
      
      # Apply parameter transformation to the derivatives (chain rule)
      variables <- intersect(variables, names)
      
      mysensitivities <- outSens$sens1[ , variables, ]
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        dPsub <- dP[controls$dimnames_sens$sens, ]
        dX <- einsum::einsum("aik,kj->aij", mysensitivities, dPsub)
        dimnames(dX) <- list(NULL, variables, colnames(dPsub))
      } else {
        dX <- mysensitivities
      }
      
    } else {
      
      outSens2 <- suppressWarnings(
        .Call(paste0("solve_", as.character(extended2)),
              as.numeric(times),
              as.numeric(mypars),
              as.numeric(optionsSens2$atol),
              as.numeric(optionsSens2$rtol),
              as.integer(optionsSens2$maxattemps),
              as.integer(optionsSens2$maxsteps),
              as.numeric(optionsSens2$roottol),
              as.integer(optionsSens2$maxroot))
      )
      
      colnames(outSens2$variable) <- controls$dimnames_sens$variable
      
      dimnames(outSens2$sens1) <- list(NULL, controls$dimnames_sens$variable, controls$dimnames_sens$sens)
      dimnames(outSens2$sens2) <- list(NULL, controls$dimnames_sens$variable, controls$dimnames_sens$sens, controls$dimnames_sens$sens)
      
      out <- cbind(outSens2$time, submatrix(outSens2$variable, cols = names))
      colnames(out)[1] <- "time"
      
      # Apply parameter transformation to the derivatives (chain rule)
      variables <- intersect(variables, names)
      
      mysensitivities  <- outSens2$sens1[, variables, ]
      mysensitivities2 <- outSens2$sens2[, variables, , ]
      
      # Extract parameter derivatives (first and second order)
      dP  <- attr(pars, "deriv",  exact = TRUE)
      dP2 <- attr(pars, "deriv2", exact = TRUE)
      
      if (!is.null(dP) && !is.null(dP2)) {
        # Subset relevant inner parameters (theta)
        dPsub  <- dP [controls$dimnames_sens$sens, , drop = FALSE]
        dP2sub <- dP2[controls$dimnames_sens$sens, , , drop = FALSE]
        
        # --- First-order chain rule: dX/dθ = dX/dP * dP/dθ
        dX <- einsum::einsum("aik,kj->aij", mysensitivities, dPsub)
        
        # --- Second-order chain rule:
        # d²X/dθ² = d²X/dP² * dP/dθ * dP/dθ + dX/dθ * dP²/d²θ
        # term1: d²X/dP² * dP/dθ * dP/dθ
        # term2: dX/dP * dP²/d²θ
        term1 <- einsum::einsum("aikl,kj,lm->aijm", mysensitivities2, dPsub, dPsub)
        term2 <- einsum::einsum("aik,kmj->aijm",  mysensitivities,  dP2sub)
        dX2 <- term1 + term2
        
        # Assign dimension names
        dimnames(dX)  <- list(NULL, variables, colnames(dPsub))
        dimnames(dX2) <- list(NULL, variables, colnames(dPsub), colnames(dPsub))
        
      } else {
        # No parameter transformation provided
        dX  <- mysensitivities
        dX2 <- mysensitivities2
      }
    }
    prdframe(out, deriv = dX, sensitivities = mysensitivities, deriv2 = dX2, sensitivities2 = mysensitivities2, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- NULL
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
}


#' Model prediction function for ODE models without sensitivities. 
#' @description Interface to get an ODE 
#' into a model function \code{x(times, pars, forcings, events)} returning ODE output.
#' It is a reduced version of \link{Xs}, missing the sensitivities. 
#' @param odemodel Object of class \link{odemodel}.
#' @param forcings, see \link{Xs}
#' @param events, see \link{Xs}
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See \link[stats]{approxfun} for possible arguments.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to \code{f}. 
#' All quantities that are not initialised by pars 
#' in \code{x(..., forcings, events)} are initialized with 0. For more details and
#' the return value see \link{Xs}.
#' @export
Xf <- function(odemodel, forcings = NULL, events = NULL, condition = NULL, optionsOde=list(method = "lsoda"), fcontrol = NULL) {
  
  func <- odemodel$func
  
  myforcings <- forcings
  myevents <- events
  myfcontrol <- fcontrol
  
  variables <- attr(func, "variables")
  parameters <- attr(func, "parameters")
  yini <- rep(0,length(variables))
  names(yini) <- variables
  
  # Controls to be modified from outside
  controls <- list(
    forcings = myforcings,
    events = myevents,
    optionsOde = optionsOde,
    fctonrol = myfcontrol
  )
  
  P2X <- function(times, pars, deriv = TRUE){
    
    events <- controls$events
    forcings <- controls$forcings
    optionsOde <- controls$optionsOde
    
    # Add event time points (required by integrator) 
    event.times <- unique(events$time)
    times <- sort(union(event.times, times))
    
    
    yini[names(pars[names(pars) %in% variables])] <- pars[names(pars) %in% variables]
    mypars <- pars[parameters]
    #alltimes <- unique(sort(c(times, forctimes)))
    
    # loadDLL(func)
    if(!is.null(forcings)) forc <- setForcings(func, forcings) else forc <- NULL
    out <- suppressWarnings(do.call(odeC, c(list(y=yini, times=times, func=func, parms=mypars, forcings=forc,events = list(data = events), fcontrol = fcontrol), optionsOde)))
    #out <- cbind(out, out.inputs)      
    
    prdframe(out, deriv = NULL, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
  
  
  
  
  
  
  
}


#' Model prediction function from data.frame
#' 
#' @param data data.frame with columns "name", "time", and row names that 
#' are taken as parameter names. The data frame can contain a column "value"
#' to initialize the parameters.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class \link{prdfn}, i.e. 
#' a function \code{x(times pars, deriv = TRUE, conditions = NULL)}, 
#' see also \link{Xs}. Attributes are "parameters", the parameter names (row names of
#' the data frame), and possibly "pouter", a named numeric vector which is generated
#' from \code{data$value}.
#' @examples
#' # Generate a data.frame and corresponding prediction function
#' timesD <- seq(0, 2*pi, 0.5)
#' mydata <- data.frame(name = "A", time = timesD, value = sin(timesD), 
#'                      row.names = paste0("par", 1:length(timesD)))
#' x <- Xd(mydata)
#' 
#' # Evaluate the prediction function at different time points
#' times <- seq(0, 2*pi, 0.01)
#' pouter <- structure(mydata$value, names = rownames(mydata))
#' prediction <- x(times, pouter)
#' plot(prediction)
#' 
#' @export
Xd <- function(data, condition = NULL) {
  
  states <- unique(as.character(data$name))
  
  
  # List of prediction functions with sensitivities
  predL <- lapply(states, function(s) {
    subdata <- subset(data, as.character(name) == s)
    
    M <- diag(1, nrow(subdata), nrow(subdata))
    parameters.specific <- rownames(subdata)
    if(is.null(parameters.specific)) parameters.specific <- paste("par", s, 1:nrow(subdata), sep = "_")
    sensnames <- paste(s, parameters.specific, sep = ".")
    
    # return function
    out <- function(times, pars) {
      value <- approx(x = subdata$time, y = pars[parameters.specific], xout = times, rule = 2)$y
      grad <- do.call(cbind, lapply(1:nrow(subdata), function(i) {
        approx(x = subdata$time, y = M[, i], xout = times, rule = 2)$y
      }))
      colnames(grad) <- sensnames
      attr(value, "sensitivities") <- grad
      attr(value, "sensnames") <- sensnames
      return(value)
    }
    
    attr(out, "parameters") <- parameters.specific
    
    return(out)
    
  }); names(predL) <- states
  
  # Collect parameters
  parameters <- unlist(lapply(predL, function(p) attr(p, "parameters")))
  
  # Initialize parameters if available
  pouter <- NULL
  if(any(colnames(data) == "value")) 
    pouter <- structure(data$value[match(parameters, rownames(data))], names = parameters)
  
  sensGrid <- expand.grid(states, parameters, stringsAsFactors=FALSE)
  sensNames <- paste(sensGrid[,1], sensGrid[,2], sep=".")  
  
  
  controls <- list()  
  
  P2X <- function(times, pars, deriv=TRUE){
    
    
    predictions <- lapply(states, function(s) predL[[s]](times, pars)); names(predictions) <- states
    
    out <- cbind(times, do.call(cbind, predictions))
    colnames(out) <- c("time", states)
    
    mysensitivities <- NULL
    myderivs <- NULL
    if(deriv) {
      
      # Fill in sensitivities
      outSens <- matrix(0, nrow = length(times), ncol = length(sensNames), dimnames = list(NULL, c(sensNames)))
      for(s in states) {
        mysens <- attr(predictions[[s]], "sensitivities")
        mynames <- attr(predictions[[s]], "sensnames")
        outSens[, mynames] <- mysens
      }
      
      mysensitivities <- cbind(time = times, outSens)
      
      # Apply parameter transformation to the derivatives
      sensLong <- matrix(outSens, nrow = nrow(outSens)*length(states))
      dP <- attr(pars, "deriv")
      if (!is.null(dP)) {
        sensLong <- sensLong %*% submatrix(dP, rows = parameters)
        sensGrid <- expand.grid.alt(states, colnames(dP))
        sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      }
      outSens <- cbind(times, matrix(sensLong, nrow = dim(outSens)[1]))
      colnames(outSens) <- c("time", sensNames)
      
      myderivs <- outSens
      #attr(out, "deriv") <- outSens
    }
    
    #attr(out, "parameters") <- unique(sensGrid[,2])
    
    prdframe(out, deriv = myderivs, sensitivities = mysensitivities, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions. 
#' 
#' @description Creates an object of type \link{obsfn} that evaluates an observation function
#' and its derivatives based on the output of a model prediction function, see \link{prdfn}, 
#' as e.g. produced by \link{Xs}.
#' @param g Named character vector or equation vector defining the observation function
#' @param f Named character of equations or object that can be converted to eqnvec or object of class fn.
#' If f is provided, states and parameters are guessed from f.
#' @param states character vector, alternative definition of "states", usually the names of \code{f}. If both,
#' f and states are provided, the states argument overwrites the states derived from f.
#' @param parameters character vector, alternative definition of the "parameters",
#' usually the symbols contained in "g" and "f" except for \code{states} and the code word \code{time}. If both,
#' f and parameters are provided, the parameters argument overwrites the parameters derived from f and g.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param attach.input logical, indiating whether the original input should be
#' returned with the output.
#' @param deriv logical, generate function to evaluate derivatives of observables. Necessary for parameter estimation.
#' @param compile Logical, compile the function (see \link{funC0})
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to R command line.
#' @return Object of class \link{obsfn}, i.e.
#' a function \code{y(..., deriv = TRUE, conditions = NULL)} representing the evaluation of the 
#' observation function. Arguments \code{out} (model prediction) and \code{pars} (parameter values)
#' shoudl be passed by the \code{...} argument.
#' If \code{out} has the attribute  "sensitivities", the result of
#' \code{y(out, pars)}, will have an attributed "deriv" which reflecs the sensitivities of 
#' the observation with respect to the parameters.
#' If \code{pars} is the result of a parameter transformation \code{p(pars)} (see \link{P}), 
#' the Jacobian 
#' of the parameter transformation and the sensitivities of the observation function
#' are multiplied according to the chain rule for differentiation.
#' @details For \link{odemodel}s with forcings, it is best, to pass the prediction function \code{x} to the "f"-argument 
#' instead of the equations themselves. If an eqnvec is passed to "f" in this case, the forcings and states
#' have to be specified manually via the "states"-argument.
#' @example inst/examples/prediction.R
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL, condition = NULL, attach.input = TRUE, deriv = TRUE, compile = FALSE, modelname = NULL, verbose = FALSE) {
 
  
  # Idea: 
  # If replicate scaling is undispensible and different 
  # observable names for different replicates is not an option, then
  # g could be a list of observables. For this case, the observation
  # function has to return a list of observations for each condition.
  # Not yet clear how this works with the "+" operator.
  myattach.input <- attach.input
  
  warnings <- FALSE
  modelname_deriv <- NULL
  
  if (is.null(f) && is.null(states) && is.null(parameters)) 
    stop("Not all three arguments f, states and parameters can be NULL")
  
  # Modify modelname by condition
  if (!is.null(modelname) && !is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Then add suffix(es) for derivative function
  if (!is.null(modelname)) modelname_deriv <- paste(modelname, "deriv", sep = "_")
  
  # Get potential paramters from g, forcings are treated as parameters because
  # sensitivities dx/dp with respect to forcings are zero.
  # Distinguish between
  # symbols = any symbol that occurrs in g
  # states = states of the underlying prediction function
  # parameters = parameters of the underlying prediction function
  # estimate = parameters p for which derivatives should be returned and states x assumed to provide derivatives, dx/dp
  symbols <- getSymbols(unclass(g))
  if (is.null(f)) {
    states <- union(states, "time")
    estimate <- union(states, parameters)
    parameters <- union(parameters, setdiff(symbols, c(states, "time")))
  } else if (inherits(f, "fn")) {
    myforcings <- Reduce(union, lapply(lapply(attr(f, "mappings"), 
                                              function(mymapping) {attr(mymapping, "forcings")}), 
                                       function(myforcing) {as.character(myforcing$name)}))
    mystates <- unique(c(do.call(c, lapply(getEquations(f), names)), "time"))
    if(length(intersect(myforcings, mystates)) > 0)
      stop("Forcings and states overlap in different conditions. Please run Y for each condition by supplying only the condition specific f.")
    
    mystates <- c(mystates, myforcings)
    myparameters <- setdiff(union(getParameters(f), getSymbols(unclass(g))), c(mystates, myforcings))
    
    estimate <- c(states, parameters)
    if (is.null(states)) estimate <- c(estimate, setdiff(mystates, myforcings))
    if (is.null(parameters)) estimate <- c(estimate, myparameters)
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  } else {
    # Get all states and parameters from f
    f <- as.eqnvec(f)
    mystates <- union(names(f), "time")
    myparameters <- getSymbols(c(unclass(g), unclass(f)), exclude = mystates)
    # Set states and parameters to be estimated according to arguments, and
    # take values from mystates and myparameters, if NULL
    estimate <- c(states, parameters)
    if (is.null(states)) estimate <- c(estimate, mystates)
    if (is.null(parameters)) estimate <- c(estimate, myparameters)
    # Return states and parameters according to what is found in the equations and what is supplied by the user (probably not needed)
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  }

  # cat("States:\n")
  # print(states)
  # cat("Parameters:\n")
  # print(parameters)
  # cat("Estimate:\n")
  # print(estimate)
  
  # Observables defined by g
  observables <- names(g)
  
  gEval <- funC0(g, variables = states, parameters = parameters, compile = compile, modelname = modelname, 
                 verbose = verbose, convenient = FALSE, warnings = FALSE)
  
  # Produce everything that is needed for derivatives
  if (deriv) {
    # Character matrices of derivatives
    dxdp <- dgdx <- dgdp <- NULL    
    states.est <- intersect(states, estimate)
    pars.est <- intersect(parameters, estimate)
    
    variables.deriv <- c(
      states, 
      as.vector(outer(states.est, c(states.est, pars.est), paste, sep = "."))
    )
    
    if (length(states.est) > 0 & length(pars.est) > 0) {
      dxdp <- apply(expand.grid.alt(states.est, c(states.est, pars.est)), 1, paste, collapse = ".")
      dxdp <- matrix(dxdp, nrow = length(states.est))
    }
    if (length(states.est) > 0)
      dgdx <- matrix(jacobianSymb(g, states.est), nrow = length(g))
    if (length(pars.est) > 0) {
      dgdp <- cbind(
        matrix("0", nrow = length(g), ncol = length(states.est)), 
        matrix(jacobianSymb(g, pars.est), nrow = length(g))
      )
    }
    
    # Sensitivities of the observables
    derivs <- as.vector(sumSymb(prodSymb(dgdx, dxdp), dgdp))
    if (length(derivs) == 0) stop("Neither states nor parameters involved. Use Y() with argument 'deriv = FALSE' instead.")
    names(derivs) <- apply(expand.grid.alt(observables, c(states.est, pars.est)), 1, paste, collapse = ".")
    
    derivsEval <- funC0(derivs, variables = variables.deriv, parameters = parameters, compile = compile, modelname = modelname_deriv,
                        verbose = verbose, convenient = FALSE, warnings = FALSE)
    
  }
  
  # Vector with zeros for possibly missing derivatives
  # zeros <- rep(0, length(dxdp))
  # names(zeros) <- dxdp
  # Redundant -> missing values have been implemented in funC0
  
  controls <- list(attach.input = attach.input) 
  
  X2Y <- function(out, pars) {
    
    attach.input <- controls$attach.input
    
    # Prepare list for with()
    nOut <- ncol(out)
    values <- gEval(M = out, p = pars)
    
    sensitivities.export <- NULL
    myderivs <- NULL
    
    dout <- attr(out, "sensitivities")
    if (!is.null(dout) & deriv) {
      dvalues <- derivsEval(M = cbind(out, dout), p = pars)
      sensitivities.export <- cbind(time = out[, 1], dvalues)
    }
    
    
    # Parameter transformation
    dP <- attr(pars, "deriv")
    if (!is.null(dP) & !is.null(dout) & deriv) {
      
      parameters.all <- c(states.est, pars.est)
      parameters.missing <- parameters.all[!parameters.all %in% rownames(dP)]
      
      if (length(parameters.missing) > 0 & warnings)
        warning("Parameters ", paste(parameters.missing, collapse = ", ", "are missing in the Jacobian of the parameter transformation. Zeros are introduced."))
      
      dP.full <- matrix(0, nrow = length(parameters.all), ncol = ncol(dP), dimnames = list(parameters.all, colnames(dP)))
      dP.full[intersect(rownames(dP), parameters.all),] <- dP[intersect(rownames(dP), parameters.all),]

      # Multiplication with tangent map
      sensLong <- matrix(dvalues, nrow = nrow(out)*length(observables))
      sensLong <- sensLong %*% dP.full
      dvalues <- matrix(sensLong, nrow = dim(out)[1])
      
      # Naming
      sensGrid <- expand.grid.alt(observables, colnames(dP.full))
      sensNames <- paste(sensGrid[,1], sensGrid[,2], sep = ".")
      colnames(dvalues) <- sensNames
      
    }
    
    
    
    # Format output
    values <- cbind(time = out[,"time"], values)
    if (attach.input)
      values <- cbind(values, submatrix(out, cols = -1))
    
    
    myderivs <- myparameters <- NULL
    if (!is.null(dout) & deriv & !attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues)
      if (is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    if (!is.null(dout) & deriv & attach.input) {
      myderivs <- cbind(time = out[,"time"], dvalues, submatrix(attr(out, "deriv"), cols = -1))
      if (is.null(dP)) myparameters <- names(pars) else myparameters <- colnames(dP)
    }
    
    
    # Output 
    prdframe(prediction = values, deriv = myderivs, sensitivities = sensitivities.export, parameters = pars) 
    
    
    
  }
  
  attr(X2Y, "equations") <- g
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states") <- states
  attr(X2Y, "modelname") <- modelname
  
  obsfn(X2Y, parameters, condition)
  
}

#' Generate a prediction function that returns times
#' 
#' Function to deal with non-ODE models within the framework of dMod. See example.
#' 
#' @param condition  either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @return Object of class \link{prdfn}.
#' @examples 
#' x <- Xt()
#' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"))
#' 
#' times <- seq(-1, 1, by = .05)
#' pars <- c(a = .1, b = 1)
#' 
#' plot((g*x)(times, pars))
#' @export
Xt <- function(condition = NULL) {
  
  
  
  # Controls to be modified from outside
  controls <- list()
  
  P2X <- function(times, pars, deriv=TRUE){
    
    out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))
    sens <- deriv <- out
    
    prdframe(out, deriv = deriv, sensitivities = sens, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- NULL
  attr(P2X, "equations") <- NULL
  attr(P2X, "forcings") <- NULL
  attr(P2X, "events") <- NULL
  
  
  prdfn(P2X, NULL, condition) 
  
  
  
  
}



#' An identity function which vanishes upon concatenation of fns
#'
#' @return fn of class idfn
#' @export
#'
#' @examples
#' x <- Xt()
#' id <- Id()
#'
#' (id*x)(1:10, pars = c(a = 1))
#' (x*id)(1:10, pars = c(a = 1))
#' str(id*x)
#' str(x*id)
Id <- function() {
  outfn <- function() return(NULL)
  class(outfn) <- c("idfn", "fn")
  return(outfn)
}
