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


#' Compute residuals and weighted residuals with BLOQ handling (C++ backend)
#'
#' @description
#' `resCpp()` is a high-performance variant of [res()] that computes residuals
#' and weighted residuals. It supports optional BLOQ (Below Limit of Quantification)
#' handling using various methods (`"M1"`, `"M3"`, `"M4NM"`, `"M4BEAL"`).
#'
#' @details
#' ## Input requirements
#' 
#' `data` must contain columns:
#' - `time` (numeric): Time points
#' - `name` (character/factor): Observable names
#' - `value` (numeric): Measured values
#' - `sigma` (numeric): Standard errors (can contain NA if `err` is provided)
#' - `lloq` (numeric, optional): Lower limit of quantification
#'
#' `out` can be:
#' - A **data.frame** in wide format (typical ODE output): First column is `time`,
#'   remaining columns are named observables
#' - A **data.frame** in long format: Must have columns `time`, `name`, and 
#'   either `prediction` or `value`
#' - A **numeric matrix** in wide format: First column is `time`, remaining columns
#'   are named observables (column names must match data$name values)
#'
#' `err` (optional) provides error model predictions, same format as `out`.
#' Used to fill NA values in `data$sigma`.
#'
#' ## BLOQ handling
#' 
#' When `value <= lloq`, the observation is flagged as BLOQ (`bloq = TRUE`).
#' The behavior depends on `optBLOQ`:
#' 
#' - `"none"` (default): Standard residuals, equivalent to [res()]
#' - `"M1"`: BLOQ observations contribute 0 to the objective (weighted.residual = 0)
#' - `"M3"`: Uses `f^2 = -2*log(Phi(-wr))` where `wr` is the weighted residual
#' - `"M4NM"` or `"M4BEAL"`: Uses `f^2 = -2*log(1 - Phi(wr)/Phi(w0))`
#'   where `w0 = prediction/sigma`. Requires `lloq >= 0`.
#'
#' Note: For `optBLOQ = "none"`, values below LLOQ are set to LLOQ before
#' computing residuals (same as [res()]).
#'
#' @param data data.frame with columns `time`, `name`, `value`, `sigma`, 
#'   optionally `lloq`.
#' @param out Predictions as data.frame or matrix. See Details.
#' @param err Optional error model predictions (data.frame or matrix). 
#'   Used to fill NA values in `data$sigma`.
#' @param optBLOQ Character string specifying BLOQ handling method:
#'   `"none"` (default), `"M1"`, `"M3"`, `"M4NM"`, or `"M4BEAL"`.
#'
#' @return An `objframe` (data.frame with class `c("objframe", "data.frame")`) 
#'   containing:
#'   - `time`: Time points from data
#'   - `name`: Observable names from data
#'   - `value`: Measured values (set to lloq if originally below lloq)
#'   - `prediction`: Model predictions
#'   - `sigma`: Standard errors (filled from err if originally NA)
#'   - `residual`: `prediction - value`
#'   - `weighted.residual`: Weighted residuals, potentially transformed by BLOQ method
#'   - `weighted.0`: `prediction / sigma` (used internally for M4 methods)
#'   - `bloq`: Logical, TRUE if original value was below lloq
#'
#' @export
#' @family dMod interface
#' @seealso [res()] for the R implementation, [objframe()] for the output class
#' 
#' @examples
#' \dontrun{
#' # Typical usage with ODE output
#' data <- data.frame(
#'   time = c(0, 1, 2, 0, 1, 2),
#'   name = c("A", "A", "A", "B", "B", "B"),
#'   value = c(1.0, 0.8, 0.6, 0.0, 0.2, 0.3),
#'   sigma = c(0.1, 0.1, 0.1, 0.05, 0.05, 0.05),
#'   lloq = c(0, 0, 0, 0.1, 0.1, 0.1)
#' )
#' 
#' # ODE output (wide format)
#' out <- data.frame(
#'   time = c(0, 1, 2),
#'   A = c(1.0, 0.85, 0.65),
#'   B = c(0.0, 0.15, 0.35)
#' )
#' 
#' # Standard residuals (like res())
#' result <- resCpp(data, out, optBLOQ = "none")
#' 
#' # With M3 BLOQ handling
#' result_m3 <- resCpp(data, out, optBLOQ = "M3")
#' 
#' # Check which observations are BLOQ
#' subset(result_m3, bloq == TRUE)
#' }
resCpp <- function(data, out, err = NULL, optBLOQ = "none") {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  
  required_cols <- c("time", "name", "value", "sigma")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("'data' is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate optBLOQ
  valid_bloq <- c("none", "M1", "M3", "M4NM", "M4BEAL")
  if (!optBLOQ %in% valid_bloq) {
    stop("'optBLOQ' must be one of: ", paste(valid_bloq, collapse = ", "))
  }
  
  # Convert name to character if factor
  if (is.factor(data$name)) {
    data$name <- as.character(data$name)
  }
  
  # Ensure lloq exists (fill with -Inf if missing)
  if (!"lloq" %in% names(data)) {
    data$lloq <- -Inf
  }
  # Call C++ function
  result <- .Call(`_dMod_resCpp`, data, out, err, optBLOQ)
  
  # --- Transform deriv attribute to match res() output format ---
  deriv_out <- attr(out, "deriv")
  if (!is.null(deriv_out)) {
    # Extract parameter names from column names (format: observable.parameter)
    deriv_cols <- colnames(deriv_out)[-1]  # skip time column
    pars <- unique(unlist(lapply(strsplit(deriv_cols, split = ".", fixed = TRUE), function(i) i[2])))
    
    # Get unique observables from data
    names_unique <- unique(as.character(result$name))
    
    # Build sensitivity column names: observable.parameter
    sensnames <- as.vector(outer(names_unique, pars, paste, sep = "."))
    
    # Match names to corresponding sensitivities
    names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names_unique), ncol = length(pars)))
    
    # Get positions of sensnames in colnames of deriv
    sensnames.deriv <- match(sensnames, colnames(deriv_out))
    
    # Match data names to unique names
    data.name.idx <- match(as.character(result$name), names_unique)
    
    # Get the columns in deriv corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name.idx]], ncol = length(data.name.idx))
    
    # Match times
    times_unique <- sort(unique(result$time))
    data.time.idx <- match(result$time, times_unique)
    out.time.idx <- match.num(times_unique, deriv_out[, 1])
    timeIndex <- out.time.idx[data.time.idx]
    
    # Extract derivatives for each data row
    deriv.prediction <- do.call(rbind, lapply(1:nrow(result), function(i) {
      submatrix(deriv_out, timeIndex[i], derivnameIndex[, i])
    }))
    colnames(deriv.prediction) <- pars
    
    # Create deriv.data in res() format
    deriv.data <- data.frame(time = result$time, name = result$name, deriv.prediction)
    attr(result, "deriv") <- deriv.data
  }
  
  # --- Transform deriv.err attribute to match res() output format ---
  deriv_err <- attr(err, "deriv")
  if (!is.null(err) && !is.null(deriv_err)) {
    # Extract parameter names from column names
    deriv_err_cols <- colnames(deriv_err)[-1]  # skip time column
    pars <- unique(unlist(lapply(strsplit(deriv_err_cols, split = ".", fixed = TRUE), function(i) i[2])))
    
    # Get unique observables from data
    names_unique <- unique(as.character(result$name))
    
    # Build sensitivity column names
    sensnames <- as.vector(outer(names_unique, pars, paste, sep = "."))
    
    # Match names to corresponding sensitivities
    names.sensnames <- t(matrix(1:length(sensnames), nrow = length(names_unique), ncol = length(pars)))
    
    # Get positions of sensnames in colnames of deriv.err
    sensnames.deriv <- match(sensnames, colnames(deriv_err))
    
    # Match data names to unique names
    data.name.idx <- match(as.character(result$name), names_unique)
    
    # Get the columns in deriv.err corresponding to data names
    derivnameIndex <- matrix(sensnames.deriv[names.sensnames[, data.name.idx]], ncol = length(data.name.idx))
    
    # Match times
    times_unique <- sort(unique(result$time))
    data.time.idx <- match(result$time, times_unique)
    err.time.idx <- match.num(times_unique, deriv_err[, 1])
    timeIndex <- err.time.idx[data.time.idx]
    
    # Extract derivatives for each data row
    deriv.prediction <- do.call(rbind, lapply(1:nrow(result), function(i) {
      submatrix(deriv_err, timeIndex[i], derivnameIndex[, i])
    }))
    colnames(deriv.prediction) <- pars
    
    # Set NA to 0
    deriv.prediction[is.na(deriv.prediction)] <- 0
    
    # Set derivatives to 0 for rows where sigma was not NA originally
    sNAIndex <- is.na(data$sigma)
    deriv.prediction[!sNAIndex, ] <- 0
    
    # Create deriv.err.data in res() format
    deriv.err.data <- data.frame(time = result$time, name = result$name, deriv.prediction)
    attr(result, "deriv.err") <- deriv.err.data
  }
  
  return(result)
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

