#' Compare data and model prediction by computing residuals
#'
#' When sigma is given in the data, it always has priority over the sigma in the error model
#'
#' @param data data.frame with name (factor), time (numeric), value (numeric) and sigma (numeric)
#' @param out output of ode(), optionally augmented with attributes
#' "deriv" (output of ode() for the sensitivity equations) and
#' "parameters" (character vector of parameter names, a subsest of those
#' contained in the sensitivity equations). If "deriv" is given, also "parameters"
#' needs to be given.
#' @param err output of the error model function
#' @return data.frame with the original data augmented by columns "prediction" (
#' numeric, the model prediction), "residual" (numeric, difference between
#' prediction and data value), "weighted.residual" (numeric, residual devided
#' by sigma). If "deriv" was given, the returned data.frame has an
#' attribute "deriv" (data.frame with the derivatives of the residuals with
#' respect to the parameters).
#' @export
#' @family dMod interface
#' @importFrom stats setNames
res <- function(data, out, err = NULL) {
  
  # .. 1 Preparations to match prediction values with data values ----#
  data$name <- as.character(data$name)
  # Unique times, names and parameter names
  times <- sort(unique(data$time))
  names <- unique(data$name)
  # Match data times/names in unique times/names
  data.time <- match.num(data$time, times)
  data.name <- match(data$name, names)
  # Match unique times/names in out times/names
  out.time <- match.num(times, out[,1])
  out.name <- match(names, colnames(out))
  if (any(is.na(out.name))) stop("The following observable in data does not have a prediction: ", paste0(setdiff(names, colnames(out)), collapse = ","))
  # Match data times/names in out times/names
  timeIndex <- out.time[data.time]
  nameIndex <- out.name[data.name]
  prediction <- sapply(1:nrow(data), function(i) out[timeIndex[i], nameIndex[i]])
  
  
  # .. Propagate derivatives if available ----#
  deriv <- attr(out, "deriv")
  deriv.data <- NULL
  if (!is.null(deriv)) {
    pars <- unique(unlist(lapply(strsplit(colnames(deriv)[-1], split = ".", fixed = TRUE), function(i) i[2])))
    sensnames <- as.vector(outer(names, pars, paste, sep = "."))
    # Match names to the corresponding sensitivities in sensnames
    names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)))
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv))
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
    # Derivatives of the prediction
    deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) submatrix(deriv, timeIndex[i], derivnameIndex[, i])))
    colnames(deriv.prediction) <- pars
    
    deriv.data <- data.frame(time = data$time, name = data$name, deriv.prediction)
  }
  
  
  # .. Modifications if error model is available ----#
  # * There are six cases to consider
  #   * 1 all(!is.na(data$sigma)) &  is.null(err)  --> only data$sigma counts
  #   * 2 all(!is.na(data$sigma)) & !is.null(err)  --> only data$sigma counts
  #   * 3 any_not_all(!is.na(data$sigma)) &  is.null(err)  --> invalid
  #   * 4 any_not_all(!is.na(data$sigma)) & !is.null(err)  --> complicated
  #   * 5 all(is.na(data$sigma)) &  is.null(err)  --> invalid
  #   * 6 all(is.na(data$sigma)) & !is.null(err)  --> only err counts
  
  # .. Informative error message for cases 3 and 5 ----#
  if (any(is.na(data$sigma)) & is.null(err))
    stop("In data, some sigmas are NA and no errmodel exists for the respective condition. Please fix data$sigma or supply errmodel.")
  
  # [] validate that this does exactly what I want
  # .. This index vector is used to identify candidates for replacing sigmas by errmodel-predictions ----#
  sNAIndex <- is.na(data$sigma)
  if (!any(sNAIndex))
    err <- NULL # if all sigmas are present, err can be thrown away, since no results of it are used anyway
  
  if (!is.null(err)) {
    time.err <- match.num(times, err[,1])
    name.err <- match(names, colnames(err))
    timeIndex <- time.err[data.time]
    nameIndex <- name.err[data.name]
    errprediction <- sapply(1:nrow(data), function(i) err[timeIndex[i], nameIndex[i]])
    if (any(sNAIndex & is.na(errprediction)))
      stop("errmodel predicts NA for some observables with is.na(data$sigma).")
    data$sigma[sNAIndex] <- errprediction[sNAIndex]
  }
  
  # .. Propagate derivatives of err model if available ----#
  deriv.err <- attr(err, "deriv")
  deriv.err.data <- NULL
  if (!is.null(err) && !is.null(deriv.err)) {
    
    pars <- unique(unlist(lapply(strsplit(colnames(deriv.err)[-1], split = ".", fixed = TRUE), function(i) i[2])))
    sensnames <- as.vector(outer(names, pars, paste, sep = "."))
    # Match names to the corresponding sensitivities in sensnames
    names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names), ncol = length(pars)))
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv.err))
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name]], ncol = length(data.name))
    # Derivatives of the prediction
    deriv.prediction <- do.call(rbind, lapply(1:nrow(data), function(i) submatrix(deriv.err, timeIndex[i], derivnameIndex[, i])))
    colnames(deriv.prediction) <- pars
    deriv.prediction[is.na(deriv.prediction)] <- 0
    
    deriv.prediction[!sNAIndex, ] <- 0 # set the derivatives of unused error-predictions to zero.
    
    deriv.err.data <- data.frame(time = data$time, name = data$name, deriv.prediction)
    
  }
  
  # .. Set value to loq if below loq ----#
  data$value <- pmax(data$value, data$lloq)
  is.bloq <- data$value <= data$lloq
  
  # .. Compute residuals ----#
  residuals <- prediction - data$value
  weighted.residuals <- (prediction - data$value)/data$sigma
  weighted.0 <- prediction/data$sigma
  
  data[["prediction"]] <- prediction
  data[["residual"]] <- residuals
  data[["weighted.residual"]] <- weighted.residuals
  data[["weighted.0"]] <- weighted.0
  data[["bloq"]] <- is.bloq
  
  # .. output ----#
  objframe(data, deriv = deriv.data, deriv.err = deriv.err.data)
  
}


#' Compute residuals and Jacobian with BLOQ handling (C++ backend)
#'
#' @description
#' `resCpp()` is the high-performance C++ equivalent of [res()].  
#' It computes residuals, weighted residuals, and the combined Jacobian
#' (prediction + error model sensitivities), with optional BLOQ handling
#' via methods M1, M3, M4NM, and M4BEAL.
#'
#' @details
#' The function calls directly into C++ via Rcpp, so it is much faster
#' than a pure R implementation.  
#' If `optBLOQ = "none"`, the results are equivalent to [res()].  
#' Otherwise, BLOQ-specific transformations of residuals and Jacobian
#' entries are applied.
#'
#' @param data A data.frame with columns:
#'   * `time` (numeric)
#'   * `name` (factor or character, observable name)
#'   * `value` (numeric, observed value)
#'   * `sigma` (numeric, residual standard deviation)
#'   * `lloq` (numeric, lower limit of quantification)
#' @param out Numeric matrix of model predictions (rows = times, cols = observables).
#' @param deriv Optional numeric matrix of sensitivities of predictions
#'   w.r.t. parameters (first two columns = time and name).
#' @param deriv_err Optional numeric matrix of sensitivities of the error
#'   model w.r.t. parameters (first two columns = time and name).
#' @param optBLOQ Character string: `"none"`, `"M1"`, `"M3"`, `"M4NM"`, or `"M4BEAL"`.
#'
#' @return A data.frame of class `objframe` with columns:
#'   * `time`, `name`, `value`, `prediction`, `sigma`,
#'     `residual`, `weighted.residual`, `weighted.0`, `bloq`
#'
#'   and with an attribute `jacobian` (numeric matrix of residual sensitivities).
#'
#' @seealso [res()] for the pure-R implementation without BLOQ handling.
#'
#' @examples
#' df <- data.frame(time = 1:3,
#'                  name = c("A","A","A"),
#'                  value = c(1,2,3),
#'                  sigma = c(1,1,1),
#'                  lloq = c(0,0,0))
#' out <- matrix(c(1.1, 2.1, 3.1), ncol = 1)
#'
#' resCpp(df, out, optBLOQ = "M3")
#'
#' @export
resCpp <- function(data, out, deriv = NULL, deriv_err = NULL, optBLOQ = "none") {
  .Call(`_dMod_resCpp`, data, out, deriv, deriv_err, optBLOQ)
}




#' Time-course data for the JAK-STAT cell signaling pathway
#'
#' Phosphorylated Epo receptor (pEpoR), phosphorylated STAT in the
#' cytoplasm (tpSTAT) and total STAT (tSTAT) in the cytoplasmhave been 
#' measured at times 0, ..., 60.
#'
#' @name jakstat
#' @docType data
#' @keywords data
NULL




# Match with numeric tolerance 
match.num <- function(x, y, tol = 1e-8) {
  
  digits <- -log10(tol)
  match(round(x, digits), round(y, digits))
  
} 

