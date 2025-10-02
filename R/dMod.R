#' dMod: Dynamic Modeling in R
#'
#' @docType package
#' @name dMod
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib dMod, .registration = TRUE
## usethis namespace: end
NULL


utils::globalVariables(c("value", "sigma", "condition", "x", "y", "name", "proflist", "delta", "loq", "bloq",
                  "sigmaLS", "cbLower95", "cbUpper95", "cbLower68", "cbUpper68",
                  "iteration", "idx", "is.zero", "index", "converged", "iterations",
                  ".runbgOutput", "terminal", "parvalue", "dataErrorModel",
                  "Rate", "weighted.residual", "i", "whichIndex"))
