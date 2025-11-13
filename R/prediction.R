

#' Model prediction function for ODE models. 
#' @description Interface to combine an ODE and its sensitivity equations
#' into one model function `x(times, pars, deriv = TRUE)` returning ODE output and sensitivities.
#' @param odemodel object of class 'odemodel' or 'odemodel++', see [odemodel]
#' @param forcings data.frame with columns name (factor), time (numeric) and value (numeric).
#' The ODE forcings. Not (yet) implemented for boost::odeint::rosenbrock4
#' @param events data.frame of events with columns "var" (character, the name of the state to be
#' affected), "time" (numeric, time point), "value" (numeric, value), "method" (character, either
#' "replace", "add" or "multiply"). See [events][deSolve::events].
#' ATTENTION: Sensitivities for event states will only be correctly computed if defined within
#' [odemodel()]. Specify events within `Xs()` only for forward simulation.
#' @param names character vector with the states to be returned. If NULL, all states are returned.
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param optionsSens list with arguments to be passed to odeC() for integration of the extended system
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See [approxfun][stats::approxfun] for possible arguments.
#' @return Object of class [prdfn]. If the function is called with parameters that
#' result from a parameter transformation (see [P]), the Jacobian of the parameter transformation
#' and the sensitivities of the ODE are multiplied according to the chain rule for
#' differentiation. The result is saved in the attributed "deriv", 
#' i.e. in this case the attibutes "deriv" and "sensitivities" do not coincide. 
#' @export
Xs <- function(odemodel, ...) {
  UseMethod("Xs", odemodel)
}

#' @export
#' @importFrom einsum einsum
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
    
    return(prdframe(out, deriv = dX, parameters = pars))
  }
  
  attr(P2X, "parameters") <- c(variables, parameters)
  attr(P2X, "equations") <- as.eqnvec(attr(func, "equations"))
  attr(P2X, "forcings") <- forcings
  attr(P2X, "events") <- events
  attr(P2X, "modelname") <- func[1]
  
  
  prdfn(P2X, c(variables, parameters), condition) 
}

#' Reshape ODE sensitivities from wide matrix (var.par) to 3D array
#'
#' @description
#' Converts a flat sensitivity matrix (as produced by deSolve-based ODE integrations)
#' with column names of the form "variable.parameter" into a structured 3D array.
#'
#' @param sensMatrix A numeric matrix of sensitivities with column names
#'   formatted as "variable.parameter". Rows correspond to time points.
#' @param variables Character vector of state variable names.
#' @param parameters Character vector of parameter names.
#'
#' @return A numeric 3D array with dimensions: time, variable, parameter,
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
        
        # --- First-order chain rule: ∂x_i/∂θ_j = ∂x_i/∂p_k * ∂p_k/∂θ_j
        dX <- einsum::einsum("aik,kj->aij", mysensitivities, dPsub)
        
        # --- Second-order chain rule:
        # ∂²x_i/∂θ_k∂θ_j =
        #   (∂²x_i/∂p_a∂p_b) * (∂p_b/∂θ_k) * (∂p_a/∂θ_j) # term1
        # + (∂x_i/∂p_a) * (∂²p_a/∂θ_k∂θ_j)               # term2
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
    prdframe(out, deriv = dX, deriv2 = dX2, parameters = pars)
    
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
#' into a model function `x(times, pars, forcings, events)` returning ODE output.
#' It is a reduced version of [Xs], missing the sensitivities. 
#' @param odemodel Object of class [odemodel].
#' @param forcings, see [Xs]
#' @param events, see [Xs]
#' @param condition either NULL (generic prediction for any condition) or a character, denoting
#' the condition for which the function makes a prediction.
#' @param optionsOde list with arguments to be passed to odeC() for the ODE integration.
#' @param fcontrol list with additional fine-tuning arguments for the forcing interpolation. 
#' See [approxfun][stats::approxfun] for possible arguments.
#' @details Can be used to integrate additional quantities, e.g. fluxes, by adding them to `f`. 
#' All quantities that are not initialised by pars 
#' in `x(..., forcings, events)` are initialized with 0. For more details and
#' the return value see [Xs].
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
#' @return Object of class [prdfn], i.e. 
#' a function `x(times pars, deriv = TRUE, conditions = NULL)`, 
#' see also [Xs]. Attributes are "parameters", the parameter names (row names of
#' the data frame), and possibly "pouter", a named numeric vector which is generated
#' from `data$value`.
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
    
    prdframe(out, deriv = myderivs, parameters = pars)
    
  }
  
  attr(P2X, "parameters") <- structure(parameters, names = NULL)
  attr(P2X, "pouter") <- pouter
  
  prdfn(P2X, attr(P2X, "parameters"), condition)
  
}


#' Observation functions. 
#' 
#' @description 
#' Creates an object of type [obsfn] that evaluates an observation function
#' and, if requested, its first and second derivatives based on the output of a model 
#' prediction function, see [prdfn], as e.g. produced by [Xs].
#' 
#' @param g Named character vector or [eqnvec] defining the observation function.
#' @param f Named character vector of equations or an object that can be converted 
#' to [eqnvec], or an object of class 'fn'. If `f` is provided, states and parameters 
#' are automatically inferred from `f`.
#' @param states Character vector, alternative definition of state variables, usually 
#' the names of `f`. If both `f` and `states` are provided, the `states` argument 
#' overrides those derived from `f`.
#' @param parameters Character vector, alternative definition of parameters, usually 
#' the symbols contained in `g` and `f` except for `states` and the keyword `time`. 
#' If both `f` and `parameters` are provided, the `parameters` argument overrides those 
#' derived from `f` and `g`.
#' @param condition Either `NULL` (generic prediction for any condition) or a character 
#' string specifying the condition for which the function generates predictions.
#' @param attach.input Logical, indicating whether the original model input should be 
#' included in the output.
#' @param deriv Logical, if `TRUE`, the function evaluates first-order derivatives
#' of observables with respect to parameters.
#' @param deriv2 Logical, if `TRUE`, the function also evaluates second derivatives 
#' of observables with respect to parameters.
#' @param compile Logical, if `TRUE`, the function is compiled (see [CppODE::funCpp]).
#' @param modelname Character, used if `compile = TRUE`, specifies a fixed filename 
#' for the generated C file.
#' @param verbose Logical, print compiler output to the R console.
#' 
#' @return 
#' An object of class [obsfn], i.e. a function 
#' `g(..., deriv = TRUE, deriv2 = FALSE, condition = NULL, verbose = F)` representing the evaluation of the 
#' observation function. The function returns observable values and, if requested, 
#' their first- and second-order derivatives with respect to the parameters.
#' 
#' @example inst/examples/prediction.R
#' 
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @importFrom abind abind
#' @export
Y <- function(g, f = NULL, states = NULL, parameters = NULL, condition = NULL,
              attach.input = TRUE, deriv = TRUE, deriv2 = TRUE,
              compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  if (deriv2 && !deriv) {
    warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
    deriv <- TRUE
  }
  
  if (is.null(f) && is.null(states) && is.null(parameters)) 
    stop("Not all three arguments f, states and parameters can be NULL")
  
  # --- Define model name with condition suffix (to avoid name collisions) ---
  if (is.null(modelname)) modelname <- "obsfn"
  if (!is.null(condition)) modelname <- paste("obsfn",modelname, sanitizeConditions(condition), sep = "_")
  
  # --- Identify symbols in g ---
  symbols <- getSymbols(unclass(g))
  
  # --- Infer states and parameters ---
  if (is.null(f)) {
    states <- union(states, "time")
    parameters <- union(parameters, setdiff(symbols, states))
  } else if (inherits(f, "fn")) {
    myforcings <- Reduce(union, lapply(lapply(attr(f, "mappings"), 
                                              function(m) attr(m, "forcings")), 
                                       function(ff) as.character(ff$name)))
    mystates <- unique(c(do.call(c, lapply(getEquations(f), names)), "time"))
    if (length(intersect(myforcings, mystates)) > 0)
      stop("Forcings and states overlap in different conditions.")
    
    mystates <- c(mystates, myforcings)
    myparameters <- setdiff(union(getParameters(f), getSymbols(unclass(g))), 
                            c(mystates, myforcings))
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  } else {
    f <- as.eqnvec(f)
    mystates <- union(names(f), "time")
    myparameters <- getSymbols(c(unclass(g), unclass(f)), exclude = mystates)
    states <- union(mystates, states)
    parameters <- union(myparameters, parameters)
  }
  
  observables <- names(g)
  obsParams <- intersect(symbols, parameters)
  obsStates <- setdiff(symbols, parameters)
  
  # --- Compile evaluator for g (value, Jacobian, Hessian) ---
  gEval <- CppODE::funCpp(
    g,
    variables  = obsStates,
    parameters = obsParams,
    compile    = compile,
    modelname  = modelname,
    verbose    = verbose,
    convenient = FALSE,
    warnings   = FALSE,
    deriv      = deriv,
    deriv2     = deriv2
  )
  
  controls <- list(attach.input = attach.input)
  
  # --- Core observation mapping function ---
  X2Y <- function(out, pars, deriv = TRUE, deriv2 = FALSE, env = parent.frame()) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE. Setting deriv = TRUE automatically.")
      deriv <- TRUE
    }
    
    attach.input <- controls$attach.input
    
    outEval <- gEval(out[, obsStates], 
                     pars[obsParams], 
                     deriv = deriv, 
                     deriv2 = deriv2)
    
    # --- Observable values ---
    values <- cbind(time = out[,"time"], outEval$out)
    if (attach.input) values <- cbind(values, submatrix(out, cols = -1))
    
    # --- Compute first and second derivatives ---
    myderivs <- myderivs2 <- NULL
    if (deriv && !deriv2) {
      
      # --- Direct first derivatives from g ---
      # dGdX[j,a,i] = ∂g_j/∂x_a
      # dGdP[j,b,i] = ∂g_j/∂p_b
      # dX[i,a,k]   = ∂x_a/∂θ_k
      # dP[b,k]     = ∂p_b/∂θ_k
      dGdX <- outEval$jacobian[, obsStates, , drop = FALSE]
      dGdP <- outEval$jacobian[, obsParams, , drop = FALSE]
      dX   <- attr(out,  "deriv")
      dP   <- attr(pars, "deriv")
      
      if (!is.null(dX)) dXsub <- dX[, obsStates, , drop = FALSE]
      if (!is.null(dP)) dPsub <- dP[obsParams, , drop = FALSE]
      
      outer_pars <- character(0)
      
      # ---------------------------------------------------------------------
      # CASE 1: dX ≠ NULL, dP ≠ NULL → full chain rule
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂x_a)(∂x_a/∂θ_k) + (∂g_j/∂p_b)(∂p_b/∂θ_k)
      if (!is.null(dX) && !is.null(dP)) {
        term11 <- einsum::einsum("jai,iak->ijk", dGdX, dXsub)
        term12 <- einsum::einsum("jbi,bk->ijk",  dGdP, dPsub)
        myderivs <- term11 + term12
        outer_pars <- colnames(dP)
      }
      
      # ---------------------------------------------------------------------
      # CASE 2: dX ≠ NULL, dP = NULL → dynamic + local observation parameters
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂x_a)(∂x_a/∂θ_k) + ∂g_j/∂p_local
      if (!is.null(dX) && is.null(dP)) {
        term11 <- einsum::einsum("jai,iak->ijk", dGdX, dXsub)
        dyn_params   <- intersect(obsParams, dimnames(dX)[[3]])
        local_params <- setdiff(obsParams, dyn_params)
        
        if (length(local_params) > 0) {
          term12 <- aperm(dGdP[, local_params, , drop = FALSE], c(3, 1, 2))
          myderivs <- abind::abind(term11, term12, along = 3)
          outer_pars <- c(dimnames(dX)[[3]], local_params)
        } else {
          myderivs <- term11
          outer_pars <- dimnames(dX)[[3]]
        }
      }
      
      # ---------------------------------------------------------------------
      # CASE 3: dX = NULL, dP ≠ NULL → parameter-only derivatives
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂p_b)(∂p_b/∂θ_k)
      if (is.null(dX) && !is.null(dP)) {
        myderivs <- einsum::einsum("jbi,bk->ijk", dGdP, dPsub)
        outer_pars <- colnames(dP)
      }
      
      # ---------------------------------------------------------------------
      # CASE 4: dX = NULL, dP = NULL → purely algebraic observables
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = ∂g_j/∂p_b  (direct dependence on local parameters)
      if (is.null(dX) && is.null(dP)) {
        myderivs <- aperm(dGdP, c(3, 1, 2))
        outer_pars <- obsParams
      }
      
      # --- Assign dimension names and optionally attach inputs ---
      if (!is.null(myderivs))
        dimnames(myderivs) <- list(NULL, observables, outer_pars)
      
      if (attach.input && !is.null(myderivs) && !is.null(dX)) {
        # Extract and align parameter names
        dyn_params <- dimnames(dX)[[3]]
        all_params <- outer_pars
        
        # Pad dX with zeros if myderivs has extra (local) parameters
        if (length(extra <- setdiff(all_params, dyn_params)) > 0) {
          pad <- array(0, dim = c(dim(dX)[1:2], length(extra)),
                       dimnames = list(NULL, NULL, extra))
          dX <- abind::abind(dX, pad, along = 3)
        }
        
        # Reorder parameters to match myderivs
        dX <- dX[, , all_params, drop = FALSE]
        
        # Combine along observable/state dimension
        myderivs <- abind::abind(myderivs, dX, along = 2)
      }
      
    } else if (deriv && deriv2) {
      
      # --- Direct first and second derivatives from g ---
      # dGdX[j,a,i] = ∂g_j/∂x_a,   dGdP[j,b,i] = ∂g_j/∂p_b
      # dG2dX2[j,a,b,i] = ∂²g_j/∂x_a∂x_b,   dG2dXdP[j,a,b,i] = ∂²g_j/∂x_a∂p_b
      # dG2dPdX[j,b,a,i] = ∂²g_j/∂p_b∂x_a,  dG2dP2[j,b,c,i] = ∂²g_j/∂p_b∂p_c
      dGdX <- outEval$jacobian[, obsStates, , drop = FALSE]
      dGdP <- outEval$jacobian[, obsParams, , drop = FALSE]
      dG2dX2  <- outEval$hessian[, obsStates,  obsStates,  , drop = FALSE]
      dG2dXdP <- outEval$hessian[, obsStates,  obsParams, , drop = FALSE]
      dG2dPdX <- outEval$hessian[, obsParams, obsStates,  , drop = FALSE]
      dG2dP2  <- outEval$hessian[, obsParams, obsParams, , drop = FALSE]
      
      # dX[i,a,k] = ∂x_a/∂θ_k,   dP[b,k] = ∂p_b/∂θ_k
      # dX2[i,a,k,l] = ∂²x_a/∂θ_k∂θ_l,   dP2[b,k,l] = ∂²p_b/∂θ_k∂θ_l
      dX  <- attr(out,  "deriv")
      dP  <- attr(pars, "deriv")
      dX2 <- attr(out,  "deriv2")
      dP2 <- attr(pars, "deriv2")
      
      if (!is.null(dX))  dXsub  <- dX[,  obsStates, , drop = FALSE]
      if (!is.null(dX2)) dX2sub <- dX2[, obsStates, , , drop = FALSE]
      if (!is.null(dP))  dPsub  <- dP[obsParams, , drop = FALSE]
      if (!is.null(dP2)) dP2sub <- dP2[obsParams, , , drop = FALSE]
      
      outer_pars <- character(0)
      
      # ---------------------------------------------------------------------
      # CASE 1: dX ≠ NULL, dP ≠ NULL → full chain rule
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂x_a)(∂x_a/∂θ_k) + (∂g_j/∂p_b)(∂p_b/∂θ_k)
      # ∂²g_j/∂θ_k∂θ_l = (1)+(2)+(3)+(4)+(5)+(6)
      if (!is.null(dX) && !is.null(dP)) {
        term11 <- einsum::einsum("jai,iak->ijk", dGdX, dXsub)
        term12 <- einsum::einsum("jbi,bk->ijk",  dGdP, dPsub)
        myderivs <- term11 + term12
        outer_pars <- colnames(dP)
        
        term21 <- einsum::einsum("jabi,iak,ibl->ijkl", dG2dX2, dXsub, dXsub)
        term22 <- einsum::einsum("jabi,iak,bl->ijkl",  dG2dXdP, dXsub, dPsub)
        term23 <- einsum::einsum("jbai,bk,ial->ijkl",  dG2dPdX, dPsub, dXsub)
        term24 <- einsum::einsum("jbci,bk,cl->ijkl",   dG2dP2,  dPsub, dPsub)
        term25 <- einsum::einsum("jai,iakl->ijkl",     dGdX,    dX2sub)
        term26 <- einsum::einsum("jbi,bkl->ijkl",      dGdP,    dP2sub)
        
        myderivs2 <- Reduce(`+`, list(term21, term22, term23, term24, term25, term26))
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 2: dX ≠ NULL, dP = NULL → dynamic + local observation parameters
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂x_a)(∂x_a/∂θ_k) + ∂g_j/∂p_local
      # ∂²g_j/∂θ_k∂θ_l = block Hessian [dyn θ | local obs-θ]
      if (!is.null(dX) && is.null(dP)) {
        term11 <- einsum::einsum("jai,iak->ijk", dGdX, dXsub)
        dyn_params   <- dimnames(dX)[[3]]
        local_params <- setdiff(obsParams, dyn_params)
        
        if (length(local_params) > 0) {
          term12 <- aperm(dGdP[, local_params, , drop = FALSE], c(3,1,2))
          myderivs <- abind::abind(term11, term12, along = 3)
          outer_pars <- c(dyn_params, local_params)
        } else {
          myderivs <- term11
          outer_pars <- dyn_params
        }
        
        n_i <- dim(dGdX)[3]; n_j <- length(observables)
        Kx <- length(dyn_params); Kl <- length(local_params); Ktot <- Kx + Kl
        myderivs2 <- array(0, dim = c(n_i, n_j, Ktot, Ktot))
        
        # xx-block: (1)+(5)
        term21 <- einsum::einsum("jabi,iak,ibl->ijkl", dG2dX2, dXsub, dXsub)
        term25 <- einsum::einsum("jai,iakl->ijkl",     dGdX,    dX2sub)
        myderivs2[, , 1:Kx, 1:Kx] <- term21 + term25
        
        # x–p_local, p_local–x, p_local–p_local
        if (Kl > 0) {
          myderivs2[, , 1:Kx, Kx+seq_len(Kl)] <- einsum::einsum("jabi,iak->ijkb",
                                                                dG2dXdP[, obsStates, local_params, , drop = FALSE], dXsub)
          myderivs2[, , Kx+seq_len(Kl), 1:Kx] <- einsum::einsum("jbai,ial->ijbl",
                                                                dG2dPdX[, local_params, obsStates, , drop = FALSE], dXsub)
          myderivs2[, , Kx+seq_len(Kl), Kx+seq_len(Kl)] <-
            aperm(dG2dP2[, local_params, local_params, , drop = FALSE], c(4,1,2,3))
        }
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 3: dX = NULL, dP ≠ NULL → parameter-only derivatives
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = (∂g_j/∂p_b)(∂p_b/∂θ_k)
      # ∂²g_j/∂θ_k∂θ_l = (4)+(6)
      if (is.null(dX) && !is.null(dP)) {
        myderivs <- einsum::einsum("jbi,bk->ijk", dGdP, dPsub)
        outer_pars <- colnames(dP)
        term24 <- einsum::einsum("jbci,bk,cl->ijkl", dG2dP2, dPsub, dPsub)
        term26 <- einsum::einsum("jbi,bkl->ijkl",    dGdP,   dP2sub)
        myderivs2 <- term24 + term26
        outer_pars2 <- outer_pars
      }
      
      # ---------------------------------------------------------------------
      # CASE 4: dX = NULL, dP = NULL → purely algebraic observables
      # ---------------------------------------------------------------------
      # ∂g_j/∂θ_k = ∂g_j/∂p_local,  ∂²g_j/∂θ_k∂θ_l = ∂²g_j/∂p_b∂p_c
      if (is.null(dX) && is.null(dP)) {
        myderivs  <- aperm(dGdP, c(3,1,2))
        myderivs2 <- aperm(dG2dP2, c(4,1,2,3))
        outer_pars <- obsParams
        outer_pars2 <- obsParams
      }
      
      # --- Assign dimension names and optionally attach inputs ---
      if (!is.null(myderivs))
        dimnames(myderivs) <- list(NULL, observables, outer_pars)
      if (!is.null(myderivs2))
        dimnames(myderivs2) <- list(NULL, observables, outer_pars2, outer_pars2)
      
      if (attach.input && !is.null(dX) && !is.null(myderivs)) {
        # ----- First-order: pad dX over θ if local obs-params exist -----
        dyn_params <- dimnames(dX)[[3]]
        all_params <- outer_pars
        
        if (length(extra <- setdiff(all_params, dyn_params)) > 0L) {
          pad <- array(0, dim = c(dim(dX)[1:2], length(extra)),
                       dimnames = list(NULL, NULL, extra))
          dX <- abind::abind(dX, pad, along = 3)
        }
        dX <- dX[, , all_params, drop = FALSE]
        myderivs <- abind::abind(myderivs, dX, along = 2)
        
        # ----- Second-order: build full-sized zero tensor and place dX2 block -----
        if (!is.null(myderivs2) && !is.null(dX2)) {
          dyn_params2 <- dimnames(dX2)[[3]]
          all_params2 <- outer_pars2
          
          n_i <- dim(dX2)[1]
          n_a <- dim(dX2)[2]
          Ktot <- length(all_params2)
          
          # Build full tensor with proper dimnames
          dX2_full <- array(0, dim = c(n_i, n_a, Ktot, Ktot),
                            dimnames = list(NULL,
                                            dimnames(dX2)[[2]],  # <- keep proper state names
                                            all_params2, all_params2))
          
          idx <- match(dyn_params2, all_params2)
          dX2_full[, , idx, idx] <- dX2
          
          # Concatenate along observable/state axis
          myderivs2 <- abind::abind(myderivs2, dX2_full, along = 2)
        }
      }
    }
    
    # --- Return predictions and sensitivities ---------------------------
    prdframe(prediction = values, deriv = myderivs, deriv2 = myderivs2, parameters = pars)
  }
  
  # --- Attach metadata --------------------------------------------------
  attr(X2Y, "equations")  <- as.eqnvec(g)
  attr(X2Y, "parameters") <- parameters
  attr(X2Y, "states")     <- states
  attr(X2Y, "modelname")  <- modelname
  
  obsfn(X2Y, parameters, condition)
}


 
# #' Generate a prediction function that returns times
# #' 
# #' Function to deal with non-ODE models within the framework of dMod. See example.
# #' 
# #' @param condition  either NULL (generic prediction for any condition) or a character, denoting
# #' the condition for which the function makes a prediction.
# #' @return Object of class [prdfn].
# #' @examples 
# #' x <- Xt()
# #' g <- Y(c(y = "a*time^2+b"), f = NULL, parameters = c("a", "b"))
# #' 
# #' times <- seq(-1, 1, by = .05)
# #' pars <- c(a = .1, b = 1)
# #' 
# #' plot((g*x)(times, pars))
# #' @export
# Xt <- function(condition = NULL) {
# 
# 
# 
#   # Controls to be modified from outside
#   controls <- list()
# 
#   P2X <- function(times, pars, deriv=TRUE){
# 
#     out <- matrix(times, ncol = 1, dimnames = list(NULL, "time"))
#     sens <- deriv <- out
# 
#     prdframe(out, deriv = deriv, sensitivities = sens, parameters = pars)
# 
#   }
# 
#   attr(P2X, "parameters") <- NULL
#   attr(P2X, "equations") <- NULL
#   attr(P2X, "forcings") <- NULL
#   attr(P2X, "events") <- NULL
# 
# 
#   prdfn(P2X, NULL, condition)
# 
# 
# 
# 
# }



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
