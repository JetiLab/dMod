#' Compute residuals
#' 
#' @param data Data frame with columns time, name, value, sigma, lloq
#' @param out Numeric matrix with predictions (first column = time)
#' @param err Optional error model predictions
#' @return Object of class objframe with data and derivatives
#' @export
res <- function(data, out, err = NULL) {
  
  # Call Rcpp implementation
  result <- res_cpp(data, out, err)
  
  # Convert to objframe with 3D/4D arrays
  objframe(
    result$data,
    deriv = result$deriv,
    deriv2 = result$deriv2,
    deriv.err = result$deriv.err,
    deriv2.err = result$deriv2.err
  )
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

