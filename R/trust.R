#' Non-linear optimisation via a trust-region method
#'
#' \code{trust} minimises (or maximises) a smooth objective function for which
#' value, gradient and Hessian are available. \code{trustL1} additionally adds
#' an L1 (lasso-style) penalty \code{lambda * sum(|p - mu|)} on a user-selected
#' subset of parameters, with kink-aware step handling. Both routines use a
#' Moré-Sorensen trust-region subproblem on a reduced working space that drops
#' coordinates pinned at a bound where the gradient pushes further into the
#' bound (active-set treatment). \code{trustL1} additionally drops parameters
#' pinned at their L1 kink.
#'
#' @param objfun R function whose first argument is a numeric vector of
#'   parameters. Must return a list with components \code{value},
#'   \code{gradient}, \code{hessian}. Extra arguments accepted by
#'   \code{objfun} can be supplied via \code{...}.
#' @param parinit Named numeric starting vector. Must be finite. Values
#'   outside \code{[parlower, parupper]} are clipped with a warning.
#' @param rinit Initial trust-region radius.
#' @param rmax Maximum allowed trust-region radius.
#' @param parscale Optional named or unnamed numeric of length
#'   \code{length(parinit)} for parameter rescaling. The subproblem
#'   operates on \code{g / parscale} and
#'   \code{H / outer(parscale, parscale)}.
#' @param iterlim Maximum number of outer trust-region iterations.
#' @param fterm Convergence threshold on \code{|f - f_try|}.
#' @param mterm Convergence threshold on
#'   \code{|g^T p + 0.5 * p^T H p|}.
#' @param minimize If \code{TRUE} (default) minimise; if \code{FALSE}
#'   maximise.
#' @param blather If \code{TRUE} return the per-iteration trace
#'   (\code{argpath}, \code{argtry}, \code{steptype}, \code{accept},
#'   \code{r}, \code{rho}, \code{valpath}, \code{valtry},
#'   \code{preddiff}, \code{stepnorm}).
#' @param parupper,parlower Named or scalar numeric bounds. If unnamed,
#'   the first element broadcasts to all parameters; if named, the
#'   entries slot by name into a length-K vector defaulting to
#'   \code{+/- Inf}.
#' @param printIter If \code{TRUE} print iteration count and objective
#'   value to the console at each function evaluation.
#' @param traceFile Optional path. If non-\code{NULL}, CSV-log per
#'   evaluation \code{iter, value, p1, p2, ...}.
#' @param ... Additional named arguments forwarded to \code{objfun}.
#'
#' @return A list with components \code{argument}, \code{value},
#'   \code{gradient}, \code{hessian}, \code{iterations},
#'   \code{converged}. When \code{blather = TRUE} the list also
#'   contains \code{argpath}, \code{argtry}, \code{steptype},
#'   \code{accept}, \code{r}, \code{rho}, \code{valpath},
#'   \code{valtry}, \code{preddiff}, \code{stepnorm}.
#'
#' @export
trust <- function(objfun, parinit, rinit, rmax,
                  parscale  = NULL,
                  iterlim   = 100L,
                  fterm     = 1e-6,
                  mterm     = 1e-6,
                  minimize  = TRUE,
                  blather   = FALSE,
                  parupper  = NULL,
                  parlower  = NULL,
                  printIter = FALSE,
                  traceFile = NULL,
                  ...) {
  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }
  trust_impl(fn, parinit, rinit, rmax, parscale, as.integer(iterlim),
             fterm, mterm, minimize, blather,
             parupper, parlower, printIter, traceFile)
}


#' @export
#' @rdname trust
#' @param mu Named numeric vector of reference values for the L1-penalised
#'   parameters. Names must be a subset of \code{names(parinit)}; only the
#'   named parameters receive a penalty. Defaults to a zero vector covering
#'   all of \code{parinit}.
#' @param one.sided Logical. If \code{TRUE}, the penalty is one-sided and
#'   acts as a lower wall at \code{mu}: \code{lambda * max(0, mu - p)}.
#'   Otherwise it is the symmetric \code{lambda * |p - mu|}.
#' @param lambda Strength of the L1 penalty. Either a scalar (broadcast to
#'   all entries of \code{mu}) or a named numeric aligned with \code{mu}.
trustL1 <- function(objfun, parinit, mu = 0 * parinit, one.sided = FALSE, lambda = 1,
                    rinit, rmax,
                    parscale  = NULL,
                    iterlim   = 100L,
                    fterm     = 1e-6,
                    mterm     = 1e-6,
                    minimize  = TRUE,
                    blather   = FALSE,
                    parupper  = NULL,
                    parlower  = NULL,
                    printIter = FALSE,
                    traceFile = NULL,
                    ...) {
  sanePars <- sanitizePars(parinit, list(...)$fixed)
  parinit  <- sanePars$pars
  
  if (is.null(names(parinit)))
    stop("trustL1: parinit must be a named numeric vector")
  if (length(mu) > 0L && is.null(names(mu)))
    stop("trustL1: mu must be a named numeric vector")
  
  unknown <- setdiff(names(mu), names(parinit))
  if (length(unknown) > 0L)
    stop("trustL1: mu has names not present in parinit: ",
         paste(unknown, collapse = ", "))
  
  if (length(lambda) == 1L) {
    lambda <- structure(rep(as.numeric(lambda), length(mu)),
                        names = names(mu))
  } else {
    if (is.null(names(lambda)))
      stop("trustL1: lambda must be scalar or a named numeric vector")
    if (!setequal(names(lambda), names(mu)))
      stop("trustL1: names(lambda) must equal names(mu)")
    lambda <- lambda[names(mu)]
  }
  
  dots <- list(...)
  fn <- if (length(dots) > 0L) {
    function(x) do.call(objfun, c(list(x), dots))
  } else {
    objfun
  }
  
  mu <- structure(as.numeric(mu), names = names(mu))
  lambda <- structure(as.numeric(lambda), names = names(lambda))
  
  trustL1_impl(fn, parinit, mu, lambda,
               as.logical(one.sided)[1L], rinit, rmax,
               parscale, as.integer(iterlim), fterm, mterm,
               minimize, blather, parupper, parlower, printIter, traceFile)
}

