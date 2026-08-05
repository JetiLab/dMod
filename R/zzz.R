#' Package initialization
#'
#' Declares the Python dependencies needed by SBML / PEtab import/export and
#' by `symmetryDetection()`. `reticulate::py_require()` does not start Python
#' here -- it just records the requirement so the first call into Python (e.g.
#' via `importSbml()` or `symmetryDetection()`) provisions an env that has
#' the listed packages available. Users with an existing libsbml install can
#' bypass reticulate entirely by setting the `DMOD_LIBSBML_PYTHON` environment
#' variable; see `.dmod_libsbml_python()` in `R/SBMLinterface.R`.
#'
#' @keywords internal
#' @importFrom stats setNames predict density median dnorm
#' @importFrom utils globalVariables
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("reticulate", quietly = TRUE))
    reticulate::py_require(c("python-libsbml", "sympy", "scipy", "numpy", "symengine"))
}

# Guard for optionally-Suggested feature dependencies: fail with an actionable
# message if the package a feature needs is not installed. Lets heavy/optional
# deps (Python bridge, SBML/PEtab I/O, LP solver) live in Suggests, not Imports.
.require_ns <- function(pkg, feature) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for ", feature,
         "; install it with install.packages(\"", pkg, "\").", call. = FALSE)
  invisible(TRUE)
}

# cppDE's install-time SUNDIALS/KLU flags, read from its (unexported)
# cvodeConfig environment. All-empty when unavailable.
.cppDE_config <- function() {
  empty <- list(available = FALSE, cflags = "", libs = "",
                klu_available = FALSE, klu_cflags = "", klu_libs = "")
  if (!requireNamespace("cppDE", quietly = TRUE)) return(empty)
  cfg <- get0("cvodeConfig", envir = asNamespace("cppDE"), inherits = FALSE)
  if (!is.environment(cfg)) return(empty)
  utils::modifyList(empty, as.list(cfg)[names(empty)][
    !vapply(as.list(cfg)[names(empty)], is.null, logical(1))])
}
