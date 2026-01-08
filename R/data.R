#' Compare data and model prediction by computing residuals
#'
#' When sigma is given in the data, it always has priority over the sigma in the error model
#'
#' @param data data.frame with columns: time, name, value, sigma, lloq
#' @param out matrix with predictions (first column = time)
#' @param err Nullable matrix with error model predictions
#' @export
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
  out.time <- match.num(times, out[, 1])
  out.name <- match(names, colnames(out))
  if (any(is.na(out.name))) {
    stop("The following observable in data does not have a prediction: ", 
         paste0(setdiff(names, colnames(out)), collapse = ","))
  }
  # Match data times/names in out times/names
  timeIndex <- out.time[data.time]
  nameIndex <- out.name[data.name]
  prediction <- sapply(1:nrow(data), function(i) out[timeIndex[i], nameIndex[i]])
  
  
  # .. Propagate first derivatives if available ----#
  # deriv is 3D array: [time, observable, parameter]
  deriv <- attr(out, "deriv")
  deriv.data <- NULL
  if (!is.null(deriv)) {
    dims <- dim(deriv)
    n_pars <- dims[3]
    pars <- dimnames(deriv)[[3]]
    if (is.null(pars)) {
      pars <- paste0("p", seq_len(n_pars))
    }
    
    # Map observable names to indices in deriv
    obs_names_deriv <- dimnames(deriv)[[2]]
    obs_idx_map <- match(names, obs_names_deriv)
    
    # Extract derivatives: result is matrix [n_data, n_pars]
    deriv.prediction <- matrix(0, nrow = nrow(data), ncol = n_pars)
    colnames(deriv.prediction) <- pars
    
    for (i in seq_len(nrow(data))) {
      time_idx <- timeIndex[i]
      obs_idx <- obs_idx_map[data.name[i]]
      deriv.prediction[i, ] <- deriv[time_idx, obs_idx, ]
    }
    
    deriv.data <- deriv.prediction
  }
  
  
  # .. Propagate second derivatives if available ----#
  # deriv2 is 4D array: [time, observable, par1, par2]
  deriv2 <- attr(out, "deriv2")
  deriv2.data <- NULL
  if (!is.null(deriv2)) {
    dims <- dim(deriv2)
    n_pars <- dims[3]
    pars <- dimnames(deriv2)[[3]]
    if (is.null(pars)) {
      pars <- paste0("p", seq_len(n_pars))
    }
    
    # Map observable names to indices in deriv2
    obs_names_deriv2 <- dimnames(deriv2)[[2]]
    obs_idx_map <- match(names, obs_names_deriv2)
    
    # Extract second derivatives: result is array [n_data, n_pars, n_pars]
    deriv2.array <- array(0, dim = c(nrow(data), n_pars, n_pars))
    dimnames(deriv2.array) <- list(NULL, pars, pars)
    
    for (i in seq_len(nrow(data))) {
      time_idx <- timeIndex[i]
      obs_idx <- obs_idx_map[data.name[i]]
      deriv2.array[i, , ] <- deriv2[time_idx, obs_idx, , ]
    }
    
    deriv2.data <- deriv2.array
  }
  
  
  # .. Modifications if error model is available ----#
  # .. Informative error message ----#
  if (any(is.na(data$sigma)) & is.null(err)) {
    stop("In data, some sigmas are NA and no errmodel exists. ",
         "Please fix data$sigma or supply errmodel.")
  }
  
  # .. Index vector for replacing sigmas by errmodel-predictions ----#
  sNAIndex <- is.na(data$sigma)
  if (!any(sNAIndex)) {
    err <- NULL # if all sigmas are present, err can be thrown away
  }
  
  if (!is.null(err)) {
    time.err <- match.num(times, err[, 1])
    name.err <- match(names, colnames(err))
    timeIndex.err <- time.err[data.time]
    nameIndex.err <- name.err[data.name]
    errprediction <- sapply(1:nrow(data), function(i) err[timeIndex.err[i], nameIndex.err[i]])
    if (any(sNAIndex & is.na(errprediction))) {
      stop("errmodel predicts NA for some observables with is.na(data$sigma).")
    }
    data$sigma[sNAIndex] <- errprediction[sNAIndex]
  }
  
  
  # .. Propagate first derivatives of err model if available ----#
  # deriv.err is 3D array: [time, observable, parameter]
  deriv.err <- attr(err, "deriv")
  deriv.err.data <- NULL
  if (!is.null(err) && !is.null(deriv.err)) {
    dims <- dim(deriv.err)
    n_pars <- dims[3]
    pars <- dimnames(deriv.err)[[3]]
    if (is.null(pars)) {
      pars <- paste0("p", seq_len(n_pars))
    }
    
    # Map observable names
    obs_names_deriv_err <- dimnames(deriv.err)[[2]]
    obs_idx_map <- match(names, obs_names_deriv_err)
    
    # Extract derivatives: matrix [n_data, n_pars]
    deriv.prediction <- matrix(0, nrow = nrow(data), ncol = n_pars)
    colnames(deriv.prediction) <- pars
    
    for (i in seq_len(nrow(data))) {
      time_idx <- timeIndex.err[i]
      obs_idx <- obs_idx_map[data.name[i]]
      deriv.prediction[i, ] <- deriv.err[time_idx, obs_idx, ]
    }
    
    # Set NA derivatives to 0
    deriv.prediction[is.na(deriv.prediction)] <- 0
    # Set derivatives to 0 where sigma was not NA (error model not used)
    deriv.prediction[!sNAIndex, ] <- 0
    
    deriv.err.data <- deriv.prediction
  }
  
  
  # .. Propagate second derivatives of err model if available ----#
  # deriv2.err is 4D array: [time, observable, par1, par2]
  deriv2.err <- attr(err, "deriv2")
  deriv2.err.data <- NULL
  if (!is.null(err) && !is.null(deriv2.err)) {
    dims <- dim(deriv2.err)
    n_pars <- dims[3]
    pars <- dimnames(deriv2.err)[[3]]
    if (is.null(pars)) {
      pars <- paste0("p", seq_len(n_pars))
    }
    
    # Map observable names
    obs_names_deriv2_err <- dimnames(deriv2.err)[[2]]
    obs_idx_map <- match(names, obs_names_deriv2_err)
    
    # Extract second derivatives: array [n_data, n_pars, n_pars]
    deriv2.array <- array(0, dim = c(nrow(data), n_pars, n_pars))
    dimnames(deriv2.array) <- list(NULL, pars, pars)
    
    for (i in seq_len(nrow(data))) {
      time_idx <- timeIndex.err[i]
      obs_idx <- obs_idx_map[data.name[i]]
      deriv2.array[i, , ] <- deriv2.err[time_idx, obs_idx, , ]
    }
    
    # Set NA derivatives to 0
    deriv2.array[is.na(deriv2.array)] <- 0
    # Set derivatives to 0 where sigma was not NA
    for (i in which(!sNAIndex)) {
      deriv2.array[i, , ] <- 0
    }
    
    deriv2.err.data <- deriv2.array
  }
  
  
  # .. Set value to lloq if below lloq ----#
  data$value <- pmax(data$value, data$lloq)
  is.bloq <- data$value <= data$lloq
  
  # .. Compute residuals (vectorized!) ----#
  residuals <- prediction - data$value
  weighted.residuals <- (prediction - data$value) / data$sigma
  weighted.0 <- prediction / data$sigma
  
  data[["prediction"]] <- prediction
  data[["residual"]] <- residuals
  data[["weighted.residual"]] <- weighted.residuals
  data[["weighted.0"]] <- weighted.0
  data[["bloq"]] <- is.bloq
  
  # .. Output ----#
  out_list <- list(
    data = data,
    deriv = deriv.data,
    deriv2 = deriv2.data,
    deriv.err = deriv.err.data,
    deriv2.err = deriv2.err.data
  )
  
  # Set class if objframe exists
  if (exists("objframe", mode = "function")) {
    return(objframe(data, deriv = deriv.data, deriv2 = deriv2.data,
                    deriv.err = deriv.err.data, deriv2.err = deriv2.err.data))
  } else {
    return(out_list)
  }
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


#' Time-course data for the Bile-Acid demonstration model
#'
#' @name badata
#' @docType data
#' @keywords data
NULL


# Match with numeric tolerance 
match.num <- function(x, y, tol = 1e-8) {
  
  digits <- -log10(tol)
  match(round(x, digits), round(y, digits))
  
} 

