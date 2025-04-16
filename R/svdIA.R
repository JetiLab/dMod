#' SVD-based Identifiability Analysis
#' 
#' This function performs a Singular Value Decomposition (SVD) based identifiability analysis
#' using Sobol sampling to explore parameter sets for a given prediction function (usually 'g*x*p').
#' It calculates the sensitivity matrix for each set of parameters, performs SVD,
#' and returns the numerical results.
#'
#' @param prdfun The model function for which the sensitivity matrix is calculated.
#' @param times Vector of time points for evaluation.
#' @param HC A list with the parameter bounds. The list can be in two forms:
#'            - `HC = list(parlower = lower_value, parupper = upper_value)` for the same bounds for all parameters.
#'            - `HC = list(parlower = named_numeric_vector, parupper = named_numeric_vector)` for specific bounds for each parameter.
#'            Defaults to `list(parlower = -5, parupper = 3)`.
#' @param nsamples The number of random samples to generate (default is 10).
#' @param neglect Vector of parameter names to exclude from the analysis (default is NULL).
#' @param cores Number of cores to use for parallel computation (default is 1).
#' 
#' @return A list containing:
#'   - `tilde_S`: The concatenated sensitivity matrix.
#'   - `svd_result`: The result of the Singular Value Decomposition.
#'   - `samples`: The randomly sampled parameter sets.
#'   - `singular_values`: The singular values from SVD.
#'   - `right_singular_vectors`: The right singular vectors from SVD.
#'   - `param_names`: Names of parameters included in the analysis.
#' 
#' @examples
#' # Assuming you have a model function `prdfun` and parameter bounds `HC`
#' HC <- list(parlower = -5, parupper = 5)  # Default bounds for all parameters
#' result <- SVDIA(prdfun, times, HC, nsamples = 100)
#' plotSVDs(result)
#' plotRightVector(result, 1)  # Plot the first right singular vector
#' 
#' @author Simon Beyer, \email{simon.beyer@@fdm.uni-freiburg.de}
#' @importFrom parallel mclapply
#' @export
SVDIA <- function(prdfun, times, HC = list(parlower = -5, parupper = 3), nsamples = 10, neglect = NULL, cores = 1) {
  
  # Get parameter names using the model function
  param_names <- getParameters(prdfun)
  npar <- length(param_names)
  
  # Check the structure of HC and extract bounds for parameters
  if (is.numeric(HC$parlower) && is.numeric(HC$parupper)) {
    # If HC has a single value for all parameters, expand it
    parlower <- structure(rep(HC$parlower, npar), names = param_names)
    parupper <- structure(rep(HC$parupper, npar), names = param_names)
  } else if (length(HC$parlower) == npar && length(HC$parupper) == npar) {
    # If HC provides specific bounds for each parameter, use them
    parlower <- HC$parlower
    parupper <- HC$parupper
  } else {
    stop("HC must either provide a single value for all parameters or named vectors for each parameter.")
  }
  
  # Generate random samples using runif instead of Sobol sampling
  pars_matrix <- matrix(NA, nrow = nsamples, ncol = npar)
  for (i in 1:npar) {
    pars_matrix[, i] <- runif(nsamples, min = parlower[i], max = parupper[i])
  }
  colnames(pars_matrix) <- param_names
  
  # Calculate the sensitivity matrix for each parameter set and store them
  if (cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' not available. Falling back to single core processing.")
    cores <- 1
  }
  
  tryCatch({
    if (cores > 1) {
      S_list <- parallel::mclapply(1:nsamples, function(i) {
        pars <- pars_matrix[i, ]
        tryCatch({
          S_i <- do.call(rbind, getDerivs(prdfun(times, pars, deriv = T)))
          return(as.matrix(S_i))
        }, error = function(e) {
          warning(sprintf("Failed to calculate sensitivity for parameter set %s: %s", 
                          paste(pars, collapse = ", "), e$message))
          return(NULL)
        })
      }, mc.cores = cores)
    } else {
      S_list <- lapply(1:nsamples, function(i) {
        pars <- pars_matrix[i, ]
        tryCatch({
          S_i <- do.call(rbind, getDerivs(prdfun(times, pars, deriv = T)))
          return(as.matrix(S_i))
        }, error = function(e) {
          warning(sprintf("Failed to calculate sensitivity for parameter set %s: %s", 
                          paste(pars, collapse = ", "), e$message))
          return(NULL)
        })
      })
    }
  }, error = function(e) {
    stop(paste("Error in sensitivity calculation:", e$message))
  })
  
  # Remove NULL entries (failed calculations) from S_list
  S_list <- Filter(Negate(is.null), S_list)
  
  if (length(S_list) == 0) {
    stop("All sensitivity calculations failed. Please check your model and parameters.")
  } else if (length(S_list) < nsamples) {
    warning(sprintf("Only %d out of %d sensitivity calculations succeeded.", length(S_list), nsamples))
  }
  
  # Concatenate all sensitivity matrices vertically to form tilde_S
  tilde_S <- do.call(rbind, S_list)
  tilde_S <- tilde_S[, !(colnames(tilde_S) == "time")]
  
  tilde_S <- as.data.frame(tilde_S) %>%
    tibble::rownames_to_column("row") %>%
    tidyr::pivot_longer(cols = -row, names_to = "colname", values_to = "value") %>%
    dplyr::mutate(
      prefix = sub("\\..*$", "", colname),
      parameter = sub("^[^.]+\\.", "", colname)
    ) %>%
    dplyr::select(-colname) %>%
    tidyr::pivot_wider(
      names_from = parameter,
      values_from = value
    ) %>%
    dplyr::select(-row, -prefix) %>% 
    as.matrix()
  

  
  # Remove neglected parameters
  if (!is.null(neglect)) {
    tilde_S <- tilde_S[, !(param_names %in% neglect)]
    param_names <- setdiff(param_names, neglect)
  }
  
  # Perform SVD
  svd_result <- svd(tilde_S)
  singular_values <- svd_result$d
  right_singular_vectors <- svd_result$v
  
  # Return only the numerical results
  return(list(
    samples = pars_matrix,
    singular_values = singular_values,
    right_singular_vectors = right_singular_vectors
  ))
}

#' Plot Singular Values from SVDIA Result
#'
#' This function plots the singular values from an SVDIA analysis result,
#' optionally filtering by upper and lower bounds.
#'
#' @param SVDIAResult The result from the SVDIA function.
#' @param sUpper Upper bound for singular values to plot (default is Inf).
#' @param sLower Lower bound for singular values to plot (default is 0).
#'
#' @return A ggplot object showing the singular values.
#'
#' @examples
#' result <- SVDIA(prdfun, times)
#' plotSVDs(result)
#' plotSVDs(result, sUpper = 10, sLower = 1e-3)
#' 
#' @author Simon Beyer, \email{simon.beyer@@fdm.uni-freiburg.de}
#' @import ggplot2
#' @export
plotSVDs <- function(SVDIAResult, sUpper = Inf, sLower = 0) {
  # Extract singular values
  singular_values <- SVDIAResult$singular_values
  
  # Filter singular values based on bounds
  sv_df <- data.frame(
    Index = 1:length(singular_values),
    SingularValue = singular_values
  )
  
  # Apply filters
  sv_df_filtered <- sv_df[sv_df$SingularValue >= sLower & sv_df$SingularValue <= sUpper, ]
  
  if (nrow(sv_df_filtered) == 0) {
    warning("No singular values within the specified bounds.")
    return(NULL)
  }
  
  # Create the plot
  p_singularvalues <- ggplot(sv_df_filtered, aes(x = Index, y = SingularValue)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_line(alpha = 0.4) +
    scale_y_log10() +
    labs(title = "Singular Values of the Sensitivity Matrix",
         x = "Index",
         y = "Singular Value (log-scale)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p_singularvalues)
}

#' Plot Right Singular Vector from SVDIA Result
#'
#' This function plots a specific right singular vector from an SVDIA analysis result.
#'
#' @param SVDIAResult The result from the SVDIA function.
#' @param index The index of the right singular vector to plot.
#'
#' @return A ggplot object showing the right singular vector contributions.
#'
#' @examples
#' result <- SVDIA(prdfun, times)
#' plotRightVector(result, 1)  # Plot the first right singular vector
#' 
#' @author Simon Beyer, \email{simon.beyer@@fdm.uni-freiburg.de}
#' @import ggplot2
#' @export
plotRightVector <- function(SVDIAResult, index) {
  # Check if the index is valid
  if (index < 1 || index > length(SVDIAResult$singular_values)) {
    stop(paste("Invalid index. Must be between 1 and", length(SVDIAResult$singular_values)))
  }
  
  # Extract the right singular vector and corresponding singular value
  vec <- SVDIAResult$right_singular_vectors[, index]
  singular_value <- SVDIAResult$singular_values[index]
  
  # Get parameter names from the right singular vectors matrix
  # This ensures we only get names for parameters that weren't neglected
  param_names <- rownames(SVDIAResult$right_singular_vectors)
  
  # If rownames are NULL, try to determine parameter names from dimensions
  if (is.null(param_names)) {
    # Use column names from samples, but only take as many as we have in the vector
    all_param_names <- colnames(SVDIAResult$samples)
    if (length(all_param_names) >= length(vec)) {
      param_names <- all_param_names[1:length(vec)]
    } else {
      # If we can't determine proper names, use generic ones
      param_names <- paste0("Param", 1:length(vec))
    }
  }
  
  # Create data frame for plotting
  vec_df <- data.frame(
    Parameter = param_names,
    Contribution = vec
  )
  
  # Create the plot
  p_vector <- ggplot(vec_df, aes(x = Parameter, y = Contribution)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    labs(title = paste("Right Singular Vector", index),
         subtitle = paste("Singular Value:", format(singular_value, digits = 4)),
         x = "Parameter",
         y = "Contribution") +
    theme_minimal() +
    coord_flip()
  
  return(p_vector)
}
