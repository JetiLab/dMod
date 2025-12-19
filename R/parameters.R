## Functions to generate parameter transformation ----


#' Generate a parameter transformation function
#' 
#' @description
#' This function provides a unified interface for generating condition-specific
#' parameter transformations, as commonly required in ODE-based modeling workflows.
#' 
#' `P()` can operate in two modes:
#' 
#' - **Explicit mode** (`method = "explicit"`, see [Pexpl]):
#'   Inner parameters are directly computed from symbolic expressions,
#'   for example
#'   \deqn{p_{\text{inner}} = \mathrm{parfn}(p_{\text{outer}})}
#'   A common application of the explicit mode is the log-transformation
#'   \deqn{p_{\text{outer}} \mapsto \exp(p_{\text{outer}})}
#'   which ensures positive parameters.
#'
#' - **Implicit mode** (`method = "implicit"`, see [Pimpl]):
#'   Typically used to infer initial values \eqn{p_{\text{ini}}} satisfying the
#'   steady-state condition
#'   \eqn{f(p_{\text{ini}}, p_{\text{dyn}}) = 0}.
#'   This yields an overall **partially implicit mapping**
#'   \deqn{p_{\text{dyn}} \mapsto (p_{\text{ini}}, p_{\text{dyn}})}
#'   where \eqn{f} usually represents the right-hand side (RHS) of an ODE model.
#' 
#' Both transformation types can be combined with other mappings via arithmetic
#' operators (`+` and `*`) thanks to the [parfn] interface.
#'
#' @param x Either an `eqnlist` or an object of class `eqnvec` or a named character vector,
#' or a list thereof. If a list is provided, `P()` is called on each
#' element and conditions are taken from the list names.
#' @param method Character, either `"explicit"` or `"implicit"`.
#' Determines whether to use [Pexpl] or [Pimpl].
#' @param parameters Character vector of outer parameters.
#' @param deriv Logical. If `TRUE`, compute and attach the Jacobian of the
#' transformation as attribute `"deriv"`.
#' @param deriv2 Logical. If `TRUE`, compute and attach the Hessian as
#' attribute `"deriv2"`. Implies `deriv = TRUE`.
#' @param fixed Character vector of parameter names treated as fixed (no
#' derivatives returned with respect to them).
#' @param keep.root Logical. Applies to `method = "implicit"` only.
#' If `TRUE`, reuse the root from the previous call as a warm-start guess
#' for faster convergence.
#' @param positive Logical or named logical vector. Applies to `method = "implicit"` only.
#' If `TRUE`, the steady-state solver is performed in log-space
#' (\eqn{p_{\text{ini}} = \exp(z_{\text{ini}})}) to enforce positive solutions.
#' If a named logical vector, only the states with `TRUE` values are constrained
#' to be positive.
#' @param parlower Named numeric vector of lower bounds for dependent states.
#'   Used for multistart initialization on cold starts. Names must match
#'   dependent state names. States not mentioned use default bounds
#'   (1e-6 for positive states, -1e6 otherwise).
#' @param parupper Named numeric vector of upper bounds for dependent states.
#'   Used for multistart initialization on cold starts. Names must match
#'   dependent state names. States not mentioned use default bound of 1e6.
#' @param nstart Integer. Number of random starting points for multistart
#'   root finding (default 100). Only used on cold starts when `parlower`
#'   and `parupper` are provided.
#' @param optionsSolver List. Applies to `method = "implicit"` only.
#'   Options passed to [nleqslv::nleqslv]. Recognized entries include:
#'   - `method`: solver method, either `"Newton"` or `"Broyden"` (default `"Newton"`)
#'   - `global`: globalization strategy, one of `"dbldog"` (double dogleg, default),
#'     `"pwldog"` (Powell dogleg), `"qline"` (quadratic line search), 
#'     `"gline"` (geometric line search), `"none"` (pure local method)
#'   - `control`: list with `ftol` (function tolerance, default `1e-10`),
#'     `xtol` (step tolerance, default `1e-10`), `maxit` (max iterations, default `200`)
#'   - Top-level `ftol`, `xtol`, `maxit` are also accepted for convenience.
#' @param attach.input Logical. Attach input parameters to the output if they
#' are not overwritten by the transformation (identity mapping).
#' @param condition Character. Condition label for which the transformation is
#' generated. If `trafo` is a list, this is inferred from list names.
#' @param compile Logical. If `TRUE`, compile the transformation via
#' [funCpp] for improved performance.
#' @param modelname Character. Used when `compile = TRUE` to define the
#' base name of the generated C code file.
#' @param verbose Logical. Print compiler output and diagnostic messages to the
#' R console.
#'
#' @return
#' An object of class [parfn], representing the parameter transformation.
#' The returned function
#' `p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE)`
#' computes inner parameters and attaches derivatives as attributes
#' `"deriv"` and `"deriv2"` when requested.
#'
#' @seealso
#' [Pexpl], [Pimpl], [parfn]
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @export
P <- function(x = NULL,
              method = c("explicit", "implicit"), 
              parameters=NULL, 
              deriv = TRUE,
              deriv2 = FALSE,
              fixed = NULL,
              keep.root = TRUE, 
              positive = TRUE,
              parlower     = NULL,
              parupper     = NULL,
              nstart       = 100L,
              optionsSolver = list(),
              attach.input = FALSE,
              condition = NULL, 
              compile = FALSE, 
              modelname = NULL, 
              verbose = FALSE) {
  
  method <- match.arg(method)
  if (is.null(x)) return()
  if (inherits(x, "eqnlist")) {
    method <- "implicit"
    fnout <- Pimpl(x, 
                   parameters = parameters,
                   deriv = deriv,
                   deriv2 = deriv2,
                   fixed = fixed,
                   keep.root = keep.root,
                   positive = positive,
                   parlower = parlower,
                   parupper = parupper,
                   nstart = nstart,
                   optionsSolver = optionsSolver,
                   attach.input = attach.input,
                   condition = condition,
                   compile = compile, 
                   modelname = modelname, 
                   verbose = verbose)
  } else {
    
    if (!is.list(x)) {
      x <- list(x)
      names(x) <- condition
    }
    
    fnout <- Reduce("+", lapply(1:length(x), function(i) {
      switch(method, 
             explicit = Pexpl(trafo = as.eqnvec(x[[i]]), 
                              parameters = parameters,
                              deriv = deriv,
                              deriv2 = deriv2,
                              fixed = fixed,
                              attach.input = attach.input,
                              condition = names(x[i]),
                              compile = compile,
                              modelname = modelname, 
                              verbose = verbose),
             implicit = Pimpl(x = as.eqnvec(x[[i]]), 
                              parameters = parameters,
                              deriv = deriv,
                              deriv2 = deriv2,
                              fixed = fixed,
                              keep.root = keep.root,
                              positive = positive,
                              parlower = parlower,
                              parupper = parupper,
                              nstart = nstart,
                              optionsSolver = optionsSolver,
                              attach.input = attach.input,
                              condition = names(x[i]),
                              compile = compile, 
                              modelname = modelname, 
                              verbose = verbose))
    }))
  }
  return(fnout)
}


#' Parameter transformation (explicit)
#'
#' Constructs a parameter transformation function that maps **outer parameters**
#' \eqn{p_{\text{outer}}} to **inner parameters** \eqn{p_{\text{inner}}}
#' according to symbolic expressions.
#'
#' @description
#' The explicit parameter transformation defines a direct, algebraic mapping
#'
#' \deqn{p_{\text{inner}} = \mathrm{parfn}(p_{\text{outer}}),}
#'
#' where \eqn{\mathrm{parfn}} is a vector-valued function composed from symbolic
#' expressions. Each element of `trafo` defines one component of
#' \eqn{p_{\text{inner}}}.
#'
#' Derivatives are obtained by **symbolic differentiation**:
#'
#' - The **Jacobian**
#'   \deqn{J_{ij} = \dfrac{\partial p_{\text{inner},i}}{\partial p_{\text{outer},j}}}
#'   is computed from the transformation expressions.
#'
#' - If `deriv2 = TRUE`, the **Hessian tensor**
#'   \deqn{H_{ijk} = \dfrac{\partial^2 p_{\text{inner},i}}{\partial p_{\text{outer},j}\,\partial p_{\text{outer},k}}}
#'   is also precomputed symbolically.
#'
#' These derivatives are attached as attributes `"deriv"` and `"deriv2"`
#' to the resulting function output and automatically composed when
#' transformations are combined via the [parfn] interface.
#'
#' @param trafo `eqnvec` or named character vector.
#' Names correspond to **inner parameters**; each element defines how it depends
#' on **outer parameters**.
#' @param parameters Character vector of outer parameter names. If omitted,
#' all symbols in `trafo` are used.
#' @param deriv Logical. If `TRUE`, compute and attach the Jacobian of the transformation.
#' @param deriv2 Logical. If `TRUE`, compute and attach the Hessian as well.
#' Implies `deriv = TRUE`.
#' @param attach.input Logical. If `TRUE`, include unchanged input parameters
#' in the output vector (identity mapping).
#' @param compile Logical. If `TRUE`, compile the transformation via [funCpp]
#' for faster evaluation.
#' @param condition Character label for which the transformation is generated.
#' @param modelname Base name for generated C++ code if `compile = TRUE`.
#' @param verbose Logical. Print compiler messages.
#'
#' @return
#' A function
#' `p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE, verbose = FALSE)`
#' that evaluates the parameter transformation.
#' The result is an object of class [parvec], which contains
#'
#' - the transformed parameters (`numeric` vector)
#' - attribute `"deriv"`: the Jacobian matrix
#' - attribute `"deriv2"`: the Hessian tensor (if requested)
#'
#' @seealso
#' [Pimpl] for implicit (steady-state) parameter transformations,
#' [P] for automatic mode selection.
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @export
Pexpl <- function(trafo,
                  parameters = NULL,
                  deriv = TRUE,
                  deriv2 = FALSE,
                  fixed = NULL,
                  attach.input = FALSE,
                  condition = NULL,
                  compile = FALSE,
                  modelname = NULL,
                  verbose = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Determine parameter sets
  # ---------------------------------------------------------------------------
  if (is.null(parameters)) {
    parameters <- getSymbols(trafo)
  } else {
    identity <- parameters[!(parameters %in% names(trafo))]
    names(identity) <- identity
    trafo <- c(trafo, identity)
    parameters <- getSymbols(trafo)
  }
  
  # Model name with condition label
  if (is.null(modelname)) modelname <- "expl_parfn" else paste(modelname, "expl_parfn", sep = "_")
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # ---------------------------------------------------------------------------
  # Build compiled (or fallback R) evaluator for transformation
  # ---------------------------------------------------------------------------
  PEval <- suppressWarnings(
    CppODE::funCpp(
      unclass(trafo),
      variables  = NULL,
      parameters = parameters,
      fixed      = fixed,
      compile    = compile,
      modelname  = modelname,
      verbose    = verbose,
      convenient = FALSE,
      deriv      = deriv,
      deriv2     = deriv2
    )
  )
  
  # ---------------------------------------------------------------------------
  # Define returned parameter transformation function
  # ---------------------------------------------------------------------------
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, verbose = FALSE) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE; enabling deriv = TRUE.")
      deriv <- TRUE
    }
    
    pars <- c(pars, fixed)
    
    # Evaluate inner parameters (fixed handled by funCpp at runtime)
    outPEval <- PEval(NULL, pars, fixed = names(fixed), attach.input = attach.input,
                      deriv = deriv, deriv2 = deriv2, verbose = verbose)
    
    pinner <- setNames(as.numeric(outPEval$out), colnames(outPEval$out))
    
    if (any(is.nan(pinner))) {
      stop(
        paste0(
          "The following inner parameter(s) evaluate to NaN:\n\t",
          paste0(names(pinner)[is.nan(pinner)], collapse = "\n\t"),
          ".\nLikely cause: division by zero or missing inputs."
        )
      )
    }
    
    # ----------------- Apply chain rules -----------------
    myderiv  <- NULL
    myderiv2 <- NULL
    
    # Get upstream derivatives if present
    dP  <- attr(pars, "deriv",  exact = TRUE)
    dP2 <- attr(pars, "deriv2", exact = TRUE)
    has_upstream <- !is.null(dP)
    
    if (deriv && !is.null(outPEval$jacobian)) {
      J_outer <- outPEval$jacobian[, , 1, drop = FALSE]
      dim(J_outer) <- dim(J_outer)[1:2]
      dimnames(J_outer) <- dimnames(outPEval$jacobian)[1:2]
      
      if (has_upstream) {
        # Chain rule: J_inner = J_outer %*% dP
        nonfixed <- colnames(J_outer)
        dP_sub <- dP[nonfixed, , drop = FALSE]
        myderiv <- einsum::einsum("ij,jk->ik", J_outer, dP_sub)
      } else {
        myderiv <- J_outer
      }
    }
    
    if (deriv2 && !is.null(outPEval$hessian)) {
      H_outer <- outPEval$hessian[, , , 1, drop = FALSE]
      dim(H_outer) <- dim(H_outer)[1:3]
      dimnames(H_outer) <- dimnames(outPEval$hessian)[1:3]
      
      if (has_upstream && !is.null(dP2)) {
        # Chain rule for Hessian
        nonfixed <- colnames(J_outer)
        dP_sub  <- dP[nonfixed, , drop = FALSE]
        dP2_sub <- dP2[nonfixed, , , drop = FALSE]
        
        # term1: H_outer contracted with two Jacobians
        term1 <- einsum::einsum("imn,mj,nk->ijk", H_outer, dP_sub, dP_sub)
        # term2: J_outer contracted with upstream Hessian
        term2 <- einsum::einsum("im,mjk->ijk", J_outer, dP2_sub)
        myderiv2 <- term1 + term2
      } else {
        myderiv2 <- H_outer
      }
    }
    
    # -------------------------------------------------------------------------
    # Assemble result and return
    # -------------------------------------------------------------------------
    pinner <- as.parvec(pinner,
                        deriv  = if (deriv)  myderiv  else FALSE,
                        deriv2 = if (deriv2) myderiv2 else FALSE)
    
    if (attach.input && !all(names(pars) %in% names(pinner))) {
      pinner <- c(pinner,
                  as.parvec(pars[setdiff(names(pars), names(pinner))],
                            deriv  = if (deriv)  NULL else FALSE,
                            deriv2 = if (deriv2) NULL else FALSE))
    }
    
    pinner
  }
  
  # ---------------------------------------------------------------------------
  # Attach metadata
  # ---------------------------------------------------------------------------
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
  attr(p2p, "modelname")  <- modelname
  
  parfn(p2p, parameters, condition)
}


#' Parameter transformation (implicit)
#'
#' @description
#' #' Constructs an **implicit parameter transformation** that defines a subset of
#' inner parameters \eqn{p_{\text{ini}}} as steady-state solutions of a nonlinear,
#' **autonomous** system
#'
#' \deqn{f(p_{\text{ini}}, p_{\text{dyn}}) = 0,}
#'
#' where \eqn{f} typically represents the right-hand side (RHS) of an ODE model
#' \eqn{\dot{x} = f(x, p)}. The system must be *autonomous*, i.e. it should not **explicitly** 
#' depend on time. Time-dependent (non-autonomous) systems have no stationary solutions and
#' are therefore incompatible with `Pimpl()`.
#'
#' The remaining parameters \eqn{p_{\text{dyn}}} are passed through unchanged,
#' yielding an overall **partially implicit mapping**
#'
#' \deqn{p_{\text{dyn}} \mapsto (p_{\text{ini}}, p_{\text{dyn}}).}
#'
#' The steady-state is obtained numerically using [nleqslv::nleqslv] with
#' analytical Jacobians derived via symbolic differentiation.
#' If `positive = TRUE`, the system is solved in log-space,
#' \eqn{p_{\text{ini}} = \exp(z_{\text{ini}})}, to ensure positivity of steady states.
#'
#' ## Derivatives
#'
#' Sensitivities of the implicit mapping are derived analytically via the
#' **implicit function theorem** applied to \eqn{f(p_{\text{ini}}, p_{\text{dyn}})=0}.
#'
#' **First-order derivative (Jacobian)**
#'
#' The Jacobian of the mapping \eqn{p_{\text{ini},i}(p_{\text{dyn}})} satisfies
#'
#' \deqn{
#' \dfrac{\partial f_i}{\partial p_{\text{ini},k}}
#' \dfrac{\partial p_{\text{ini},k}}{\partial p_{\text{dyn},j}}
#' +
#' \dfrac{\partial f_i}{\partial p_{\text{dyn},j}}
#' = 0,
#' }
#'
#' which implies
#'
#' \deqn{
#' \dfrac{\partial p_{\text{ini},i}}{\partial p_{\text{dyn},j}}
#' =
#' -\!
#' \left(\dfrac{\partial f_k}{\partial p_{\text{ini},i}}\right)^{-1}
#' \dfrac{\partial f_k}{\partial p_{\text{dyn},j}}.
#' }
#'
#' **Second-order derivative (Hessian)**
#'
#' A second differentiation with respect to \eqn{p_{\text{dyn},k}} yields
#'
#' \deqn{
#' \dfrac{\partial^2 p_{\text{ini},i}}{\partial p_{\text{dyn},j}\,\partial p_{\text{dyn},k}}
#' =
#' -\!
#' \left(\dfrac{\partial f_l}{\partial p_{\text{ini},i}}\right)^{-1}
#' \!\!
#' \left[
#'   \dfrac{\partial^2 f_l}{\partial p_{\text{ini},m}\,\partial p_{\text{ini},n}}
#'     \dfrac{\partial p_{\text{ini},m}}{\partial p_{\text{dyn},j}}
#'     \dfrac{\partial p_{\text{ini},n}}{\partial p_{\text{dyn},k}}
#'   +
#'   \dfrac{\partial^2 f_l}{\partial p_{\text{ini},m}\,\partial p_{\text{dyn},k}}
#'     \dfrac{\partial p_{\text{ini},m}}{\partial p_{\text{dyn},j}}
#'   +
#'   \dfrac{\partial^2 f_l}{\partial p_{\text{dyn},j}\,\partial p_{\text{ini},m}}
#'     \dfrac{\partial p_{\text{ini},m}}{\partial p_{\text{dyn},k}}
#'   +
#'   \dfrac{\partial^2 f_l}{\partial p_{\text{dyn},j}\,\partial p_{\text{dyn},k}}
#' \right]\!,
#' }
#'
#' where repeated indices imply summation (Einstein convention), and
#' all quantities are evaluated at the steady state satisfying
#' \eqn{f_i(p_{\text{ini}}, p_{\text{dyn}}) = 0}.
#' The Hessian tensor \eqn{\partial^2 p_{\text{ini},a} / (\partial p_{\text{dyn},j}\,\partial p_{\text{dyn},k})}
#' is symmetric in indices \eqn{(j,k)} for sufficiently smooth \eqn{f_i}.
#' 
#' ## Conservation relations
#'
#' When the input `trafo` is an `eqnlist` (ODE model), any algebraic conservation
#' laws in the stoichiometric matrix are automatically detected via
#' [conservedQuantities()] and converted into total quantities such as
#' `totAKT = AKT + pAKT`.  
#'
#' Each total represents a **new outer parameter** (added to \eqn{p_{\text{dyn}}})
#' that defines the conserved amount of a species pool. The corresponding
#' conservation relations are incorporated into the system as additional
#' algebraic equations of the form
#'
#' \deqn{(\text{sum of species}) - \text{tot} = 0,}
#'
#' ensuring that the steady-state computation enforces the conservation laws
#' exactly.  
#'
#' All detected conservation relations and their associated total parameters
#' are printed to the console when `Pimpl()` is called.
#'
#' ## Solver options
#'
#' The nonlinear system is solved with [nleqslv::nleqslv] using analytical
#' Jacobians. See `optionsSolver` parameter for customization.
#'
#' ## Composition
#'
#' The resulting mapping composes naturally with other parameter transformations
#' created by [Pexpl] or [Pimpl]. Derivatives (`"deriv"`, `"deriv2"`) are automatically
#' propagated when chained via [parfn].
#'
#' @param x An `eqnvec` defining
#'   \eqn{f(p_{\text{ini}}, p_{\text{dyn}})=0},
#'   or an `eqnlist` specifying an ODE model in stoichiometric form.
#' @param parameters Character vector of outer parameter names
#'   \eqn{p_{\text{dyn}}}. If omitted, inferred automatically.
#' @param deriv Logical. If `TRUE`, compute and attach the Jacobian via
#'   the implicit function theorem.
#' @param deriv2 Logical. If `TRUE`, compute and attach the Hessian.
#'   Implies `deriv = TRUE`.
#' @param fixed Character vector of parameter names treated as fixed
#'   (excluded from symbolic differentiation at compile time).
#' @param keep.root Logical. Reuse the previous root as starting guess
#'   to accelerate convergence.
#' @param positive Logical or named logical vector. If `TRUE`, solve in log-space
#'   (\eqn{p_{\text{ini}} = \exp(z_{\text{ini}})}) to ensure positive steady states
#'   for all dependent states. If a named logical vector, only the states with
#'   `TRUE` values are constrained to be positive. Names must match dependent
#'   state names. States not mentioned default to `FALSE` (unconstrained).
#' @param parlower Named numeric vector of lower bounds for dependent states.
#'   Used for multistart initialization on cold starts. Names must match
#'   dependent state names. States not mentioned use default bounds
#'   (1e-6 for positive states, -1e6 otherwise).
#' @param parupper Named numeric vector of upper bounds for dependent states.
#'   Used for multistart initialization on cold starts. Names must match
#'   dependent state names. States not mentioned use default bound of 1e6.
#' @param nstart Integer. Number of random starting points for multistart
#'   root finding (default 10). Only used on cold starts when `parlower`
#'   and `parupper` are provided.
#' @param optionsSolver List of options passed to [nleqslv::nleqslv].
#'   Merged with internal defaults via [modifyList()]. Recognized entries include:
#'   - `method`: solver method, either `"Newton"` or `"Broyden"` (default `"Newton"`)
#'   - `global`: globalization strategy, one of `"dbldog"` (double dogleg, default),
#'     `"pwldog"` (Powell dogleg), `"qline"` (quadratic line search), 
#'     `"gline"` (geometric line search), `"none"` (pure local method)
#'   - `control`: list with `ftol` (function tolerance, default `1e-10`),
#'     `xtol` (step tolerance, default `1e-10`), `maxit` (max iterations, default `200`)
#'   - Top-level `ftol`, `xtol`, `maxit` are also accepted for convenience.
#' @param attach.input Logical. Include unchanged input parameters in the output.
#' @param condition Character. Condition label for which the transformation is generated.
#' @param compile Logical. Compile generated residual functions via [funCpp].
#' @param modelname Character. Base name for compiled code if `compile = TRUE`.
#' @param verbose Logical. If TRUE, print diagnostic and convergence information.
#'
#' @return
#' A function  
#' `p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE, verbose = FALSE)`
#' where `fixed` is a named numeric vector of parameters to fix at runtime
#' (excluded from derivatives but used for evaluation).
#' Returns the steady-state parameter vector with attached attributes
#' `"deriv"` (Jacobian) and `"deriv2"` (Hessian, if requested).
#'
#' @seealso
#' [Pexpl] for explicit transformations,  
#' [nleqslv::nleqslv] for numerical steady-state solving,  
#' [P] for the unified high-level interface.
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @importFrom nleqslv nleqslv
#' @importFrom abind abind
#' @export
Pimpl <- function(x,
                  parameters   = NULL,
                  deriv        = TRUE,
                  deriv2       = FALSE,
                  fixed        = NULL,
                  keep.root    = TRUE,
                  positive     = TRUE,
                  parlower     = NULL,
                  parupper     = NULL,
                  nstart       = 10L,
                  optionsSolver = list(),
                  attach.input = FALSE,
                  condition    = NULL,
                  compile      = FALSE,
                  modelname    = NULL,
                  verbose      = FALSE) {
  
  # Normalize input: accept eqnlist; insert conservation totals if needed
  totals <- character(0)
  forcing_states <- character(0)
  
  if (inherits(x, "eqnlist")) {
    eq    <- x
    trafo <- as.eqnvec(eq)
    
    # Detect forcing/event states: states with RHS = 0 are always in steady state
    # and should be treated as parameters (pass-through)
    trafo_vec <- as.character(trafo)
    names(trafo_vec) <- names(trafo)
    
    for (st in names(trafo_vec)) {
      rhs <- trafo_vec[st]
      # Check if RHS is effectively zero: "0", "1*(0)", "-1*(0)+1*(0)", etc.
      rhs_clean <- gsub("\\s+", "", rhs)
      rhs_clean <- gsub("^[+-]?[0-9]*\\*?\\(0\\)([+-][0-9]*\\*?\\(0\\))*$", "0", rhs_clean)
      
      if (rhs_clean == "0" || rhs_clean == "0.0" || rhs_clean == "") {
        forcing_states <- c(forcing_states, st)
      }
    }
    
    # Remove forcing states from trafo - they become parameters
    if (length(forcing_states) > 0L) {
      trafo <- trafo[setdiff(names(trafo), forcing_states)]
      message("Detected forcing/event states (treated as parameters): ", 
              paste(forcing_states, collapse = ", "))
    }
    
    cq_df <- conservedQuantities(eq$smatrix)
    n_cq  <- if (!is.null(attr(cq_df, "n"))) attr(cq_df, "n") else 0L
    
    states_detected <- names(trafo)
    symbols_all     <- getSymbols(trafo)
    nonstates_auto  <- setdiff(symbols_all, states_detected)
    
    if (n_cq > 0L) {
      if (!is.null(parameters)) {
        cand_totals <- setdiff(parameters, nonstates_auto)
        totals <- if (length(cand_totals) == n_cq) cand_totals else character(0)
      } else {
        totals <- character(0)
      }
      
      # Heuristic for descriptive total names (e.g. ERK + ppERK → totERK)
      totals_created <- character(n_cq)
      for (i in seq_len(n_cq)) {
        cons_expr <- cq_df[i, 1]
        symbols_i <- getSymbols(as.character(cons_expr))
        
        if (length(symbols_i) >= 2) {
          # Determine common prefix among species names
          prefix <- Reduce(function(a, b) {
            n <- nchar(a)
            while (n > 0 && substr(b, 1, n) != substr(a, 1, n)) n <- n - 1
            substr(a, 1, n)
          }, symbols_i)
          
          # If no prefix found, fall back to longest common substring
          if (nchar(prefix) == 0) {
            common_sub <- function(a, b) {
              max_sub <- ""
              for (ii in seq_len(nchar(a))) {
                for (jj in ii:nchar(a)) {
                  sub <- substr(a, ii, jj)
                  if (grepl(sub, b) && nchar(sub) > nchar(max_sub)) max_sub <- sub
                }
              }
              max_sub
            }
            prefix <- Reduce(common_sub, symbols_i)
          }
          
          prefix <- gsub("_+$", "", prefix)
          prefix <- gsub("^[^A-Za-z0-9]+", "", prefix)
          prefix <- if (nchar(prefix) == 0) paste(symbols_i, collapse = "_") else prefix
          
          totals_created[i] <- paste0("tot", prefix)
        } else {
          totals_created[i] <- paste0("total", i)
        }
      }
      
      totals <- totals_created
      
      # Replace residuals by conservation constraints "(sum) - total_i = 0"
      # For each conservation relation, replace the ODE of ONE involved species
      # with the algebraic constraint. We choose the last species in each relation.
      
      # Convert to plain character vector for reliable manipulation
      trafo_vec <- as.character(trafo)
      names(trafo_vec) <- names(trafo)
      
      replaced_states <- character(n_cq)
      for (i in seq_len(n_cq)) {
        cons_expr <- as.character(cq_df[i, 1])
        symbols_i <- getSymbols(cons_expr)
        
        # Select one species from this conservation relation to replace
        # Filter to those that are actually states, preserving order from symbols_i
        valid_species <- symbols_i[symbols_i %in% names(trafo_vec)]
        if (length(valid_species) == 0L) {
          stop("Conservation relation involves no model states: ", cons_expr)
        }
        # Take the last valid species
        repl_state <- valid_species[length(valid_species)]
        replaced_states[i] <- repl_state
        
        # Replace the ODE with the conservation constraint
        trafo_vec[repl_state] <- paste0("(", cons_expr, ") - ", totals[i])
      }
      
      # If user provided totals in parameters, match them to the correct conservation
      # relations based on overlap between species names and total names
      if (!is.null(parameters)) {
        cand_totals <- setdiff(parameters, nonstates_auto)
        if (length(cand_totals) == n_cq) {
          # Match user totals to conservation relations by finding overlaps
          # e.g., species {ERK, pERK, ppERK} should match "totERK"
          matched_totals <- character(n_cq)
          used <- logical(length(cand_totals))
          
          for (i in seq_len(n_cq)) {
            cons_expr <- as.character(cq_df[i, 1])
            species_i <- getSymbols(cons_expr)
            
            # For each candidate total, check if any species name overlaps with it
            best_match <- NULL
            best_overlap <- 0
            
            for (j in seq_along(cand_totals)) {
              if (used[j]) next
              total_name <- cand_totals[j]
              # Remove common prefixes like "tot", "total", "Total" for comparison
              total_base <- gsub("^[Tt]ot(al)?", "", total_name)
              
              # Check overlap with each species in the conservation relation
              for (sp in species_i) {
                # Check if species contains the total base or vice versa
                # e.g., "ERK" in "totERK", or "ERK" in "pERK"
                overlap <- .string_overlap(sp, total_base)
                if (overlap > best_overlap) {
                  best_overlap <- overlap
                  best_match <- j
                }
              }
            }
            
            if (!is.null(best_match) && best_overlap > 0) {
              matched_totals[i] <- cand_totals[best_match]
              used[best_match] <- TRUE
            }
          }
          
          # If all matched successfully, use user-provided names
          if (all(nchar(matched_totals) > 0)) {
            # Update trafo with matched total names
            for (i in seq_len(n_cq)) {
              cons_expr <- as.character(cq_df[i, 1])
              trafo_vec[replaced_states[i]] <- paste0("(", cons_expr, ") - ", matched_totals[i])
            }
            totals <- matched_totals
          }
        }
      }
      
      # Convert back to eqnvec
      trafo <- as.eqnvec(trafo_vec)
      
      # Recompute symbol sets and exclude totals from nonstates
      states_detected <- names(trafo)
      symbols_all     <- getSymbols(trafo)
      nonstates_auto  <- setdiff(symbols_all, c(states_detected, totals))
      
      # Report detected conservation relations
      msg <- paste(sprintf("  • %s = %s (replacing %s)", 
                           totals, cq_df[[1]], replaced_states), collapse = "\n")
      message("Detected conservation relations:\n", msg)
    }
    
    required_outer <- unique(c(nonstates_auto, totals, forcing_states))
    if (!is.null(parameters)) {
      miss  <- setdiff(required_outer, parameters)
      extra <- setdiff(parameters, required_outer)
      if (length(miss) || length(extra)) {
        stop(
          "Pimpl(eqnlist): mismatch in 'parameters'.\n",
          if (length(miss))  paste0("  Missing: ", paste(miss,  collapse = ", "), "\n") else "",
          if (length(extra)) paste0("  Unexpected: ", paste(extra, collapse = ", "), "\n") else "",
          "  Required outer parameters are: ", paste(required_outer, collapse = ", ")
        )
      }
    } else {
      parameters <- required_outer
    }
  } else {
    trafo <- as.eqnvec(x)
  }
  
  # Partition symbols: states vs outer parameters
  states    <- names(trafo)
  nonstates <- setdiff(getSymbols(trafo), states)
  # Set parameters if not provided (for eqnvec input)
  if (is.null(parameters)) {
    parameters <- nonstates
  }
  
  # Allow some states to be independent by listing them in `parameters`
  indep_st <- intersect(states, if (is.null(parameters)) character(0) else parameters)
  dep_st   <- setdiff(states, indep_st)
  if (length(dep_st) == 0L) {
    stop("Pimpl: no dependent states to solve; 'parameters' covers all states.")
  }
  
  # Validate: dep_st cannot be compile-time fixed (we need dfdx_dep for IFT)
  if (!is.null(fixed) && any(fixed %in% dep_st)) {
    stop("Dependent states cannot be compile-time fixed: ",
         paste(intersect(fixed, dep_st), collapse = ", "))
  }
  
  # Remove forcing states from 'positive' if present (they are not solved)
  if (length(forcing_states) > 0L && !is.null(names(positive))) {
    positive <- positive[setdiff(names(positive), forcing_states)]
  }
  
  # Normalize 'positive' argument to a named logical vector over dep_st
  pos_mask <- .normalize_positive_mask(positive, dep_st)
  any_positive <- any(pos_mask)
  
  # Build compiled evaluator for the dependent residuals only
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition))
    modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # Parameters that will be symbolically differentiated (excluding compile-time fixed)
  all_params <- c(states, nonstates)
  diff_params <- setdiff(all_params, fixed)
  
  FEval <- suppressWarnings(
    CppODE::funCpp(
      as.eqnvec(trafo[dep_st]),
      variables  = NULL,
      parameters = all_params,
      fixed      = fixed,
      compile    = compile,
      modelname  = modelname,
      verbose    = verbose,
      convenient = FALSE,
      deriv      = TRUE,  # Always need Jacobian for nleqslv
      deriv2     = isTRUE(deriv2)
    )
  )
  
  # nleqslv default options (merged with user list)
  defaultsSolver <- list(
    method = "Newton",
    global = "dbldog",
    control = list(
      ftol   = 1e-6,
      xtol   = 1e-6,
      maxit  = 1000
    )
  )
  
  # Merge user options
  optsSolver <- defaultsSolver
  if (!is.null(optionsSolver$method)) optsSolver$method <- optionsSolver$method
  if (!is.null(optionsSolver$global)) optsSolver$global <- optionsSolver$global
  if (!is.null(optionsSolver$control)) {
    optsSolver$control <- modifyList(optsSolver$control, optionsSolver$control)
  }
  # Also allow top-level ftol, xtol, maxit for convenience
  if (!is.null(optionsSolver$ftol)) optsSolver$control$ftol <- optionsSolver$ftol
  if (!is.null(optionsSolver$xtol)) optsSolver$control$xtol <- optionsSolver$xtol
  if (!is.null(optionsSolver$maxit)) optsSolver$control$maxit <- optionsSolver$maxit
  
  # Warm-start state
  guess_env <- new.env(parent = emptyenv())
  guess_env$guess <- NULL
  
  # Store compile-time fixed for use in p2p
  compile_fixed <- fixed
  
  # Returned transformation closure
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, verbose = FALSE) {
    
    # Cold-Start bei neuem Fit-Kontext (z.B. mstrust)
    current_token <- getOption(".dMod.fit_token", default = NULL)
    if (!is.null(current_token) && !identical(guess_env$token, current_token)) {
      if (verbose) cat("[Pimpl] New fit context detected, resetting warm-start\n")
      guess_env$guess <- NULL
      guess_env$token <- current_token
    }
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE; enabling deriv = TRUE.")
      deriv <- TRUE
    }
    
    # Upstream derivative attributes (for chained transformations)
    dP  <- attr(pars, "deriv",  exact = TRUE)
    dP2 <- attr(pars, "deriv2", exact = TRUE)
    
    # Merge runtime fixed values into p
    p <- pars
    if (!is.null(fixed)) {
      if (is.null(names(fixed))) stop("Runtime fixed must be a named numeric vector.")
      if (any(names(fixed) %in% names(p))) {
        p[names(fixed)] <- fixed
      } else {
        p <- c(p, fixed)
      }
    }
    
    # Runtime fixed names (for filtering derivatives)
    runtime_fixed_names <- if (!is.null(fixed)) names(fixed) else character(0)
    
    # Provide initial guesses for dependent states; validate others
    if (!all(dep_st %in% names(p))) p[setdiff(dep_st, names(p))] <- 1
    if (!all(indep_st %in% names(p))) {
      stop("Missing independent states: ", paste(setdiff(indep_st, names(p)), collapse = ", "))
    }
    if (!all(nonstates %in% names(p))) {
      stop("Missing outer parameters: ", paste(setdiff(nonstates, names(p)), collapse = ", "))
    }
    
    x_dep0 <- p[dep_st]
    x_ind  <- p[indep_st]
    r_ns   <- p[nonstates]
    
    # Warm-start
    if (keep.root && !is.null(guess_env$guess)) {
      x_dep0[names(x_dep0) %in% names(guess_env$guess)] <-
        guess_env$guess[names(x_dep0) %in% names(guess_env$guess)]
    }
    
    # Helper to pack state+nonstate vector in FEval's expected order
    pack_full <- function(x_dep) {
      x_full <- setNames(numeric(length(states)), states)
      x_full[dep_st]   <- x_dep
      x_full[indep_st] <- x_ind
      c(x_full, r_ns)
    }
    
    # Define residual and Jacobian functions for nleqslv
    # These need access to FEval and pack_full
    
    # Get column indices for dependent states in the Jacobian
    dep_st_idx <- match(dep_st, diff_params)
    
    if (any_positive) {
      # Transform to z-space: z[i] = log(x[i]) for positive states, z[i] = x[i] otherwise
      x_to_z <- function(x) {
        z <- x
        z[pos_mask] <- log(pmax(x[pos_mask], 1e-12))
        z
      }
      z_to_x <- function(z) {
        x <- z
        # Clamp z to avoid overflow in exp(): exp(709) ≈ 8e307, exp(710) = Inf
        z_clamped <- pmax(pmin(z[pos_mask], 700), -700)
        x[pos_mask] <- exp(z_clamped)
        x
      }
      
      # Residual in z-space
      f_z <- function(z) {
        x_dep <- z_to_x(z)
        tryCatch(
          FEval(NULL, p = pack_full(x_dep),
                attach.input = FALSE, deriv = FALSE, deriv2 = FALSE, verbose = FALSE)$out[1, ],
          error = function(e) rep(NaN, length(z))
        )
      }
      
      # Jacobian in z-space: J_z = J_x * diag(dx/dz)
      # For positive states: dx/dz = exp(z) = x
      # For unconstrained states: dx/dz = 1
      jac_z <- function(z) {
        x_dep <- z_to_x(z)
        E <- tryCatch(
          FEval(NULL, p = pack_full(x_dep),
                attach.input = FALSE, deriv = TRUE, deriv2 = FALSE, verbose = FALSE),
          error = function(e) NULL
        )
        if (is.null(E)) return(NULL)
        J_full <- E$jacobian[, , 1, drop = TRUE]
        J_x <- J_full[, dep_st_idx, drop = FALSE]
        
        # Chain rule: multiply columns by dx/dz
        dxdz <- ifelse(pos_mask, x_dep, 1)
        sweep(J_x, 2, dxdz, `*`)
      }
      
      # Multistart if cold start and bounds are provided
      if (is.null(guess_env$guess) && !is.null(parlower) && !is.null(parupper)) {
        root <- .multistart_nleqslv(f_z, jac_z, x_dep0, dep_st, pos_mask, 
                                    parlower, parupper, nstart, optsSolver,
                                    x_to_z, z_to_x, verbose)
        x_dep_star <- root$x_star
      } else {
        z0 <- x_to_z(x_dep0)
        root <- .quiet_nleqslv(f_z, jac_z, z0, optsSolver)
        x_dep_star <- z_to_x(root$x)
      }
      
      # Check convergence
      if (!root$converged || any(is.na(x_dep_star)) || any(is.nan(x_dep_star))) {
        stop(sprintf(
          "[Pimpl] Steady-state solver did not converge (termcd: %d, fval: %.2e).\n  %s\n  Consider providing better initial guesses via the input parameter vector or using multistart (parlower/parupper).",
          root$termcd, max(abs(root$fvec)), .nleqslv_message(root$termcd)
        ))
      }
      
      if (verbose) {
        cat(sprintf(
          "[Pimpl] nleqslv converged in %d iterations (max|f|: %.2e, termcd: %d)\n",
          root$iter, max(abs(root$fvec)), root$termcd
        ))
        if (sum(pos_mask) < length(pos_mask)) {
          cat(sprintf("  Positivity enforced for: %s\n", 
                      paste(names(pos_mask)[pos_mask], collapse = ", ")))
        }
      }
    } else {
      # No positivity constraints: solve directly in x-space
      f_x <- function(x_dep) {
        tryCatch(
          FEval(NULL, p = pack_full(x_dep),
                attach.input = FALSE, deriv = FALSE, deriv2 = FALSE, verbose = FALSE)$out[1, ],
          error = function(e) rep(NaN, length(x_dep))
        )
      }
      
      jac_x <- function(x_dep) {
        E <- tryCatch(
          FEval(NULL, p = pack_full(x_dep),
                attach.input = FALSE, deriv = TRUE, deriv2 = FALSE, verbose = FALSE),
          error = function(e) NULL
        )
        if (is.null(E)) return(NULL)
        J_full <- E$jacobian[, , 1, drop = TRUE]
        J_full[, dep_st_idx, drop = FALSE]
      }
      
      # Multistart if cold start and bounds are provided
      if (is.null(guess_env$guess) && !is.null(parlower) && !is.null(parupper)) {
        root <- .multistart_nleqslv(f_x, jac_x, x_dep0, dep_st, pos_mask,
                                    parlower, parupper, nstart, optsSolver,
                                    NULL, NULL, verbose)
        x_dep_star <- root$x_star
      } else {
        root <- .quiet_nleqslv(f_x, jac_x, x_dep0, optsSolver)
        x_dep_star <- root$x
      }
      
      # Check convergence
      if (!root$converged || any(is.na(x_dep_star)) || any(is.nan(x_dep_star))) {
        stop(sprintf(
          "[Pimpl] Steady-state solver did not converge (termcd: %d, fval: %.2e).\n  %s\n  Consider providing better initial guesses via the input parameter vector or using multistart (parlower/parupper).",
          root$termcd, max(abs(root$fvec)), .nleqslv_message(root$termcd)
        ))
      }
      
      if (verbose) {
        cat(sprintf(
          "[Pimpl] nleqslv converged in %d iterations (max|f|: %.2e, termcd: %d)\n",
          root$iter, max(abs(root$fvec)), root$termcd
        ))
      }
    }
    
    if (keep.root) guess_env$guess <- setNames(x_dep_star, dep_st)
    
    # Evaluate Jacobian/Hessian of the dependent residuals at the solution
    # Note: we do NOT pass runtime fixed here - we need full Jacobians for IFT
    # The runtime filtering happens later via cols_nonfixed
    E <- FEval(NULL, p = pack_full(x_dep_star),
               attach.input = FALSE, deriv = TRUE, deriv2 = deriv2, verbose = verbose)
    
    # FEval Jacobian is w.r.t. diff_params (compile-time fixed already excluded)
    J <- E$jacobian[,,1, drop = TRUE]
    all_cols <- diff_params
    
    # Block partition df/d(·) - indices relative to diff_params
    dep_st_diff   <- intersect(dep_st, diff_params)
    indep_st_diff <- intersect(indep_st, diff_params)
    nonstates_diff <- intersect(nonstates, diff_params)
    
    dfdx_dep   <- J[, match(dep_st_diff, all_cols), drop = FALSE]
    dfdx_indep <- if (length(indep_st_diff)) J[, match(indep_st_diff, all_cols), drop = FALSE]
    else matrix(0, nrow = nrow(dfdx_dep), ncol = 0, dimnames = list(dep_st, character(0)))
    dfdp_ns    <- if (length(nonstates_diff)) J[, match(nonstates_diff, all_cols), drop = FALSE]
    else matrix(0, nrow = nrow(dfdx_dep), ncol = 0, dimnames = list(dep_st, character(0)))
    
    rownames(dfdx_dep) <- dep_st
    colnames(dfdx_dep) <- dep_st_diff
    if (ncol(dfdx_indep)) {
      rownames(dfdx_indep) <- dep_st
      colnames(dfdx_indep) <- indep_st_diff
    }
    if (ncol(dfdp_ns)) {
      rownames(dfdp_ns) <- dep_st
      colnames(dfdp_ns) <- nonstates_diff
    }
    
    # First-order IFT in x-space
    # Check condition number before inverting
    rcond_val <- tryCatch(rcond(dfdx_dep), error = function(e) 0)
    
    if (rcond_val < 1e-12) {
      stop("[Pimpl] Jacobian df/dx is singular or near-singular at the solution (rcond = ",
           format(rcond_val, digits = 3), "). ",
           "The steady-state problem may be ill-posed or the solution is not a regular root. ",
           "Consider providing better initial guesses or checking model structure.")
    }
    
    inv_fxd <- solve(dfdx_dep)
    
    Dx_ns   <- if (ncol(dfdp_ns)) - inv_fxd %*% dfdp_ns 
    else matrix(0, nrow = length(dep_st), ncol = 0, dimnames = list(dep_st, character(0)))
    Dx_ind  <- if (ncol(dfdx_indep)) - inv_fxd %*% dfdx_indep
    else matrix(0, nrow = length(dep_st), ncol = 0, dimnames = list(dep_st, character(0)))
    
    # Prepare the base output vector: states + nonstates (excluding auto-generated totals)
    x_full_star <- setNames(numeric(length(states)), states)
    x_full_star[dep_st]   <- x_dep_star
    x_full_star[indep_st] <- x_ind
    
    # Add forcing states (RHS = 0, pass-through from input)
    forcing_vals <- if (length(forcing_states) > 0L) p[forcing_states] else numeric(0)
    
    # Add nonstates to output, but exclude totals and forcing_states (forcing_states added separately above)
    nonstates_output <- setdiff(nonstates, c(totals, forcing_states))
    out_vec <- c(x_full_star, forcing_vals, r_ns[nonstates_output])
    
    # Columns for derivatives: exclude both compile-time and runtime fixed
    cols_nonfixed <- setdiff(names(pars), c(compile_fixed, runtime_fixed_names))
    
    # Build Jacobian over incoming parameter basis (excluding all fixed)
    J_outer <- matrix(0, nrow = length(out_vec), ncol = length(cols_nonfixed),
                      dimnames = list(names(out_vec), cols_nonfixed))
    
    # Identity for pass-through outputs already in out_vec
    passthrough <- setdiff(names(out_vec), dep_st)
    if (length(passthrough)) {
      idx <- intersect(passthrough, cols_nonfixed)
      if (length(idx)) J_outer[idx, idx] <- diag(1, length(idx))
    }
    
    # Fill derivatives for dependent states (only for non-fixed parameters)
    active_ns  <- intersect(nonstates_diff, cols_nonfixed)
    active_ind <- intersect(indep_st_diff, cols_nonfixed)
    if (length(active_ns))  J_outer[dep_st, active_ns]  <- Dx_ns[, active_ns, drop = FALSE]
    if (length(active_ind)) J_outer[dep_st, active_ind] <- Dx_ind[, active_ind, drop = FALSE]
    
    # Second-order IFT (x-space)
    H_outer <- NULL
    if (deriv2) {
      H <- E$hessian[,,,1, drop = TRUE]
      
      idx_xdep <- match(dep_st_diff, all_cols)
      idx_xind <- if (length(indep_st_diff)) match(indep_st_diff, all_cols) else integer(0)
      idx_ns   <- if (length(nonstates_diff)) match(nonstates_diff, all_cols) else integer(0)
      
      H_xx   <- H[, idx_xdep, idx_xdep, drop = FALSE]
      H_xp   <- if (length(idx_ns)) H[, idx_xdep, idx_ns, drop = FALSE] 
      else array(0, dim = c(length(dep_st), length(dep_st), 0))
      H_pp   <- if (length(idx_ns)) H[, idx_ns, idx_ns, drop = FALSE]
      else array(0, dim = c(length(dep_st), 0, 0))
      H_xind <- if (length(idx_xind)) H[, idx_xdep, idx_xind, drop = FALSE] else NULL
      H_indp <- if (length(idx_xind) && length(idx_ns)) H[, idx_xind, idx_ns, drop = FALSE] else NULL
      
      P_cols <- c(nonstates_diff, indep_st_diff)
      D_depP <- cbind(Dx_ns, Dx_ind)
      if (length(P_cols)) colnames(D_depP) <- P_cols
      
      rhs <- array(0, dim = c(length(dep_st), length(P_cols), length(P_cols)),
                   dimnames = list(dep_st, P_cols, P_cols))
      
      if (length(P_cols) > 0) {
        rhs <- rhs + einsum::einsum("iab,aj,bk->ijk", H_xx, D_depP, D_depP)
        
        if (length(nonstates_diff)) {
          D_ns <- Dx_ns
          rhs[, nonstates_diff, nonstates_diff] <- rhs[, nonstates_diff, nonstates_diff] +
            H_pp +
            einsum::einsum("iak,aj->ijk", H_xp, D_ns) +
            einsum::einsum("iaj,ak->ijk", H_xp, D_ns)
        }
        
        if (length(indep_st_diff)) {
          D_ind <- Dx_ind
          if (length(nonstates_diff)) {
            rhs[, nonstates_diff, indep_st_diff] <- rhs[, nonstates_diff, indep_st_diff] +
              einsum::einsum("iak,aj->ijk", H_xp, D_ind) +
              einsum::einsum("iab,aj,bk->ijk", H_xx, Dx_ns, D_ind)
            
            rhs[, indep_st_diff, nonstates_diff] <- rhs[, indep_st_diff, nonstates_diff] +
              einsum::einsum("iak,aj->ijk", H_indp, Dx_ns) +
              einsum::einsum("iab,aj,bk->ijk", H_xx, D_ind, Dx_ns)
          }
          
          rhs[, indep_st_diff, indep_st_diff] <- rhs[, indep_st_diff, indep_st_diff] +
            einsum::einsum("iab,aj,bk->ijk", H_xx, D_ind, D_ind)
          if (!is.null(H_xind)) {
            rhs[, indep_st_diff, indep_st_diff] <- rhs[, indep_st_diff, indep_st_diff] +
              einsum::einsum("iak,aj->ijk", H_xind, D_ind) +
              einsum::einsum("iaj,ak->ijk", H_xind, D_ind)
          }
        }
      }
      
      if (length(P_cols) > 0) {
        rhs_mat <- matrix(rhs, nrow = nrow(dfdx_dep), ncol = length(P_cols) * length(P_cols))
        sol_mat <- - inv_fxd %*% rhs_mat
        d2x_dep <- array(sol_mat, dim = c(length(dep_st), length(P_cols), length(P_cols)),
                         dimnames = list(dep_st, P_cols, P_cols))
      } else {
        d2x_dep <- array(0, dim = c(length(dep_st), 0, 0),
                         dimnames = list(dep_st, character(0), character(0)))
      }
      
      # Embed into full Hessian over current outputs and incoming cols (filtered)
      H_outer <- array(0, dim = c(length(out_vec), length(cols_nonfixed), length(cols_nonfixed)),
                       dimnames = list(names(out_vec), cols_nonfixed, cols_nonfixed))
      active_P <- intersect(P_cols, cols_nonfixed)
      if (length(active_P)) {
        H_outer[dep_st, active_P, active_P] <- d2x_dep[, active_P, active_P, drop = FALSE]
      }
    }
    
    # If attach.input = TRUE, append any incoming parameters that were not
    # part of out_vec as additional pass-through outputs and EXPAND derivatives
    if (attach.input) {
      extras <- setdiff(names(pars), names(out_vec))
      if (length(extras)) {
        # Append to the output vector
        out_vec <- c(out_vec, pars[extras])
        
        # Expand Jacobian with identity rows for pass-through extras
        addJ <- matrix(0, nrow = length(extras), ncol = ncol(J_outer),
                       dimnames = list(extras, colnames(J_outer)))
        hit  <- intersect(extras, colnames(J_outer))
        if (length(hit)) addJ[cbind(hit, hit)] <- 1
        J_outer <- rbind(J_outer, addJ)
        
        # Expand Hessian with zero rows for the new outputs
        if (deriv2) {
          addH <- array(0, dim = c(length(extras), ncol(J_outer), ncol(J_outer)),
                        dimnames = list(extras, colnames(J_outer), colnames(J_outer)))
          H_outer <- abind::abind(H_outer, addH, along = 1)
        }
      }
    }
    
    # Compose with upstream derivatives AFTER any expansion so shapes match
    myderiv <- if (!is.null(dP)) {
      dP_use <- dP[colnames(J_outer), , drop = FALSE]
      einsum::einsum("ij,jk->ik", J_outer, dP_use)
    } else J_outer
    
    myderiv2 <- NULL
    if (deriv2) {
      if (!is.null(dP) && !is.null(dP2)) {
        dP_use  <- dP [colnames(J_outer), , drop = FALSE]
        dP2_use <- dP2[colnames(J_outer), , , drop = FALSE]
        termA   <- einsum::einsum("imn,mj,nk->ijk", H_outer, dP_use, dP_use)
        termB   <- einsum::einsum("im,mjk->ijk",   myderiv, dP2_use)
        myderiv2 <- termA + termB
      } else {
        myderiv2 <- H_outer
      }
    }
    
    # Dimension names for derivatives
    if (deriv && !is.null(myderiv)) {
      rownames(myderiv) <- names(out_vec)
      if (!is.null(dP)) {
        colnames(myderiv) <- colnames(dP)
      } else {
        colnames(myderiv) <- cols_nonfixed
      }
    }
    
    if (deriv2 && !is.null(myderiv2)) {
      if (!is.null(dP)) {
        dimnames(myderiv2) <- list(names(out_vec), colnames(dP), colnames(dP))
      } else {
        dimnames(myderiv2) <- list(names(out_vec), cols_nonfixed, cols_nonfixed)
      }
    }
    
    # Build final parvec
    res <- as.parvec(out_vec,
                     deriv  = if (deriv)  myderiv  else FALSE,
                     deriv2 = if (deriv2) myderiv2 else FALSE)
    
    res
  }
  
  # Metadata
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- unique(parameters)
  attr(p2p, "modelname")  <- modelname
  
  parfn(p2p, parameters, condition)
}


#' Normalize positivity constraint specification
#'
#' Internal helper to convert the `positive` argument into a named logical
#' vector over all dependent states.
#'
#' @param positive Logical scalar or named logical vector.
#' @param dep_st Character vector of dependent state names.
#' @return Named logical vector of length `length(dep_st)`.
#' @keywords internal
.normalize_positive_mask <- function(positive, dep_st) {
  
  
  if (is.logical(positive) && length(positive) == 1L && is.null(names(positive))) {
    # Scalar TRUE/FALSE: apply to all dependent states
    return(setNames(rep(positive, length(dep_st)), dep_st))
  }
  
  if (is.logical(positive) && !is.null(names(positive))) {
    # Named logical vector: validate names and build mask
    unknown <- setdiff(names(positive), dep_st)
    if (length(unknown) > 0L) {
      warning("Unknown state names in 'positive' (ignored): ", 
              paste(unknown, collapse = ", "))
    }
    
    # Default to FALSE for states not mentioned
    pos_mask <- setNames(rep(FALSE, length(dep_st)), dep_st)
    overlap <- intersect(names(positive), dep_st)
    pos_mask[overlap] <- positive[overlap]
    return(pos_mask)
  }
  
  stop("'positive' must be TRUE, FALSE, or a named logical vector ",
       "with names matching dependent states: ",
       paste(dep_st, collapse = ", "))
}


#' Multistart root finding using nleqslv
#'
#' Internal helper that performs multiple root-finding attempts with uniformly
#' sampled starting points within given bounds. Returns the best solution found.
#'
#' @param fn Residual function to solve.
#' @param jac Jacobian function.
#' @param x_dep0 Initial guess (used as fallback and for structure).
#' @param dep_st Character vector of dependent state names.
#' @param pos_mask Named logical vector indicating which states have positivity constraints.
#' @param parlower Named numeric vector of lower bounds.
#' @param parupper Named numeric vector of upper bounds.
#' @param nstart Number of random starting points to try.
#' @param optsSolver Options for nleqslv.
#' @param x_to_z Optional transformation function (for positive constraints).
#' @param z_to_x Optional inverse transformation function.
#' @param verbose Logical; print progress information.
#' @return List with elements: x_star (best solution), converged, termcd, fvec, iter.
#' @keywords internal
.multistart_nleqslv <- function(fn, jac, x_dep0, dep_st, pos_mask, parlower, parupper,
                                nstart, optsSolver, x_to_z, z_to_x, verbose) {
  
  # Extract bounds for dependent states only
  lower <- x_dep0
  upper <- x_dep0
  
  for (nm in dep_st) {
    if (nm %in% names(parlower)) lower[nm] <- parlower[nm]
    else lower[nm] <- if (pos_mask[nm]) 1e-6 else -1e6
    
    if (nm %in% names(parupper)) upper[nm] <- parupper[nm]
    else upper[nm] <- if (pos_mask[nm]) 1e6 else 1e6
  }
  
  best_root   <- NULL
  best_fval   <- Inf
  use_transform <- !is.null(x_to_z) && !is.null(z_to_x)
  
  if (verbose) {
    cat(sprintf("[Pimpl] Multistart with %d random initial points\n", nstart))
  }
  
  for (i in seq_len(nstart)) {
    # Sample uniformly within bounds
    x_start <- setNames(runif(length(dep_st), lower, upper), dep_st)
    
    # Transform to z-space if needed
    start_val <- if (use_transform) x_to_z(x_start) else x_start
    
    # Attempt root finding
    root <- tryCatch(
      .quiet_nleqslv(fn, jac, start_val, optsSolver),
      error = function(e) NULL
    )
    
    if (!is.null(root) && !any(is.na(root$x)) && !any(is.nan(root$x))) {
      max_fval <- max(abs(root$fvec))
      if (max_fval < best_fval) {
        best_fval <- max_fval
        best_root <- root
        
        if (verbose) {
          cat(sprintf("  [%d] max|f|: %.2e (new best, termcd: %d)\n", 
                      i, max_fval, root$termcd))
        }
        
        # Early exit if converged well
        if (best_fval < 1e-10) break
      }
    }
  }
  
  # Fallback to original guess if multistart failed completely
  if (is.null(best_root)) {
    if (verbose) cat("  Multistart failed, using original guess\n")
    start_val <- if (use_transform) x_to_z(x_dep0) else x_dep0
    best_root <- .quiet_nleqslv(fn, jac, start_val, optsSolver)
  }
  
  # Return in consistent format
  x_star <- if (use_transform) z_to_x(best_root$x) else best_root$x
  
  list(
    x_star    = setNames(x_star, dep_st),
    x         = best_root$x,
    converged = best_root$converged,
    termcd    = best_root$termcd,
    fvec      = best_root$fvec,
    iter      = best_root$iter
  )
}


#' Quietly run nleqslv suppressing all output
#'
#' Internal helper that runs nleqslv::nleqslv while suppressing
#' all warnings, messages, and stdout output.
#'
#' @param fn Residual function.
#' @param jac Jacobian function.
#' @param start Initial guess vector.
#' @param optsSolver Options for nleqslv.
#' @return Result from nleqslv with added 'converged' field.
#' @keywords internal
.quiet_nleqslv <- function(fn, jac, start, optsSolver) {
  result <- suppressWarnings(suppressMessages(
    nleqslv::nleqslv(
      x       = start,
      fn      = fn,
      jac     = jac,
      method  = optsSolver$method,
      global  = optsSolver$global,
      control = optsSolver$control
    )
  ))
  
  
  # Add convenience field for convergence check
  # termcd == 1 means convergence
  result$converged <- (result$termcd == 1)
  
  result
}


#' Get human-readable message for nleqslv termination code
#'
#' @param termcd Integer termination code from nleqslv.
#' @return Character string describing the termination reason.
#' @keywords internal
.nleqslv_message <- function(termcd) {
  switch(as.character(termcd),
         "1" = "Convergence achieved",
         "2" = "x-values within tolerance, but not converged",
         "3" = "No acceptable step found during line search",
         "4" = "Iteration limit exceeded",
         "5" = "Jacobian is singular",
         "6" = "Jacobian is ill-conditioned",
         "-10" = "User interrupt",
         "Unknown termination code"
  )
}


#' Compute string overlap between two strings
#'
#' Internal helper that finds the longest common substring between two strings.
#' Used for fuzzy matching of species names to total parameter names.
#'
#' @param a First string.
#' @param b Second string.
#' @return Integer length of the longest common substring.
#' @keywords internal
.string_overlap <- function(a, b) {
  if (nchar(a) == 0 || nchar(b) == 0) return(0L)
  
  # Case-insensitive comparison
  a <- tolower(a)
  b <- tolower(b)
  
  # Check direct containment first (most common case)
  if (grepl(a, b, fixed = TRUE)) return(nchar(a))
  if (grepl(b, a, fixed = TRUE)) return(nchar(b))
  
  
  # Find longest common substring
  max_len <- 0L
  len_a <- nchar(a)
  len_b <- nchar(b)
  
  for (i in seq_len(len_a)) {
    for (j in i:len_a) {
      sub <- substr(a, i, j)
      if (nchar(sub) > max_len && grepl(sub, b, fixed = TRUE)) {
        max_len <- nchar(sub)
      }
    }
  }
  
  max_len
}


#' Construct parameter transformations
#'
#' Helper functions to construct and modify symbolic parameter transformations
#' used by prediction functions such as [P()] and [Xs()].
#'
#' The functions [define()], [insert()] and [branch()] operate exclusively on
#' the symbolic level. They are used to build transformation objects that
#' describe how *outer parameters* are expressed in terms of *inner parameters*
#' or constants.
#'
#' No model evaluation, sensitivity calculation or parameter checking is
#' performed by these functions. The resulting transformations are interpreted
#' later when prediction or objective functions are constructed.
#'
#' \describe{
#'   \item{define}{
#'     Reset or redefine a transformation rule by explicitly specifying a new
#'     right-hand side.
#'   }
#'   \item{insert}{
#'     Insert symbolic substitutions into existing transformation rules without
#'     resetting them.
#'   }
#'   \item{branch}{
#'     Duplicate a transformation for multiple conditions and optionally apply
#'     condition-specific substitutions.
#'   }
#' }
#'
#' When transformations are branched, a condition table is stored as metadata
#' (attribute \code{"tree"}) and may be used to restrict subsequent calls to
#' [define()] or [insert()] to specific conditions.
#'
#' @param trafo
#'   A named character vector, an object of class \code{eqnvec}, or a list
#'   thereof representing parameter transformations.
#'
#' @param expr
#'   Character string of the form \code{"lhs ~ rhs"} defining a symbolic
#'   transformation or substitution.
#'
#' @param table
#'   Optional data frame specifying condition-specific substitutions. Rownames
#'   identify conditions; columns correspond to parameter names.
#'
#' @param conditions
#'   Character vector of condition names. If supplied, overrides
#'   \code{rownames(table)}.
#'
#' @param apply
#'   Character string specifying whether and how entries of \code{table} are
#'   applied when branching:
#'   \describe{
#'     \item{"nothing"}{Only duplicate the transformation (default).}
#'     \item{"insert"}{Apply entries via [insert()].}
#'     \item{"define"}{Apply entries via [define()].}
#'   }
#'
#' @param conditionMatch
#'   Optional character string (regular expression). If provided, the operation
#'   is applied only to conditions whose names match this expression.
#'
#' @param ...
#'   Named values used to substitute symbols occurring in \code{expr}.
#'
#' @return
#' An object of the same type as \code{trafo}, possibly expanded to a list if
#' branching has been applied.
#'
#' @export
#' @example inst/examples/define.R
define <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(mytrafo[[i]])
    
    with(mytable, {
      args <- c(list(expr = expr, trafo = mytrafo[[i]], reset = TRUE), eval(dots))
      do.call(repar, args)
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}


#' @export
#' @rdname define
insert <- function(trafo, expr, ..., conditionMatch = NULL) {
  
  
  if (missing(trafo)) trafo <- NULL
  lookuptable <- attr(trafo, "tree")
  
  
  if (is.list(trafo) & is.null(names(trafo)))
    stop("If trafo is a list, elements must be named.")
  
  if (is.list(trafo) & !all(names(trafo) %in% rownames(lookuptable)))
    stop("If trafo is a list and contains a lookuptable (is branched from a tree), the list names must be contained in the rownames of the tree.")
  
  if (!is.list(trafo)) {
    mytrafo <- list(trafo)
  } else {
    mytrafo <- trafo
  }
  
  dots <- substitute(alist(...))
  out <- lapply(1:length(mytrafo), function(i) {
    
    .currentTrafo <- mytrafo[[i]]
    .currentSymbols <- NULL
    if (!is.null(.currentTrafo))
      .currentSymbols <- getSymbols(.currentTrafo)
    
    if (is.list(trafo)) {
      mytable <- lookuptable[names(mytrafo)[i], , drop = FALSE]
    } else {
      mytable <- lookuptable[1, , drop = FALSE]
    }
    
    if ((!is.null(conditionMatch)))
      if((!str_detect(rownames(mytable), conditionMatch))) 
        return(.currentTrafo)
    
    
    
    with(mytable, {
      .fun <- function() {
        # subset conditions by logicals expressions supplied by the dots
        dots_eval <- eval(dots)                                            # convert from substituted to language
        if (length(dots_eval) == 0) 
          return(do.call(repar, list(expr = expr, trafo = .currentTrafo)))
        
        dots_eval_eval <- lapply(dots_eval, function(i) eval.parent(i, 3)) # evaluate the language in the "mytable" frame. parent1: lapply, parent2: .fun, parent3: with
        which_logical <- vapply(dots_eval_eval, function(i) {is.logical(i)}, FUN.VALUE = vector("logical", 1)) # which of the dots are logical
        logical_dots <- dots_eval[which_logical]                           # subset to matching conditions
        matching <- do.call(c, logical_dots)
        if(!is.null(matching)) { # null means no logical dots were supplied
          if (any(!matching)) {  
            return(.currentTrafo) }}
        
        args <- c(list(expr = expr, trafo = .currentTrafo), dots_eval_eval[!which_logical]) # feed the rest of the eval'd dots into repar
        return(do.call(repar, args))
      }
      .fun()
    })
    
    
  })
  names(out) <- names(mytrafo)
  if (!is.list(trafo)) out <- out[[1]]
  attr(out, "tree") <- lookuptable
  
  return(out)
  
}


#' @export
#' @rdname define
branch <- function(
    trafo,
    table = NULL,
    conditions = rownames(table),
    apply = c("nothing", "insert", "define")) {
  
  
  apply <- match.arg(apply)
  
  # --- trivial case ----------------------------------------------------------
  if (is.null(table) && is.null(conditions))
    return(trafo)
  
  # --- normalize inputs ------------------------------------------------------
  if (is.null(conditions))
    conditions <- paste0("C", seq_len(nrow(table)))
  
  if (is.null(table))
    table <- data.frame(condition = conditions, row.names = conditions)
  
  rownames(table) <- conditions
  
  # --- branch trafo -----------------------------------------------------------
  out <- setNames(lapply(conditions, function(x) trafo), conditions)
  attr(out, "tree") <- table
  
  # --- optional application of table -----------------------------------------
  if (apply == "nothing")
    return(out)
  
  for (cn in conditions) {
    row <- table[cn, , drop = FALSE]
    
    for (par in colnames(row)) {
      val <- row[[par]]
      
      if (is.na(val))
        next
      
      expr <- paste0(par, " ~ ", val)
      
      if (apply == "insert") {
        out <- insert(out, expr, conditionMatch = cn)
      }
      
      if (apply == "define") {
        out <- define(out, expr, conditionMatch = cn)
      }
    }
  }
  
  out
}



#' Reparameterization
#' 
#' @param expr character of the form `"lhs ~ rhs"` where `rhs`
#' reparameterizes `lhs`. Both `lhs` and `rhs`
#' can contain a number of symbols whose values need to be passed by the `...` argument.
#' @param trafo character or equation vector or list thereof. The object where the replacement takes place in
#' @param ... pass symbols as named arguments
#' @param reset logical. If true, the trafo element corresponding to lhs is reset according to rhs. 
#' If false, lhs wherever it occurs in the rhs of trafo is replaced by rhs of the formula.
#' @return an equation vector with the reparameterization.
#' @details Left and right-hand side of `expr` are searched for symbols. If separated by
#' "_", symbols are recognized as such, e.g. in `Delta_x` where the symbols are 
#' "Delta" and "x". Each symbol for which values (character or numbers) are passed by the
#' `...` argument is replaced.
#' @export
#' @importFrom stats as.formula
#' @examples
#' innerpars <- letters[1:3]
#' constraints <- c(a = "b + c")
#' mycondition <- "cond1"
#' 
#' trafo <- repar("x ~ x", x = innerpars)
#' trafo <- repar("x ~ y", trafo, x = names(constraints), y = constraints)
#' trafo <- repar("x ~ exp(x)", trafo, x = innerpars)
#' trafo <- repar("x ~ x + Delta_x_condition", trafo, x = innerpars, condition = mycondition)
repar <- function(expr, trafo = NULL, ..., reset = FALSE) {
  
  if (inherits(expr, "formula")) expr <- deparse(expr)
  
  parsed.expr <- as.character(stats::as.formula(gsub("_", ":", expr, fixed = TRUE)))
  lhs <- parsed.expr[2]
  lhs.symbols <- getSymbols(lhs)
  rhs <- parsed.expr[3]
  rhs.symbols <- getSymbols(rhs)
  
  # Make sure that arguments are characters
  args <- lapply(list(...), as.character)
  
  replacements <- as.data.frame(args, stringsAsFactors = FALSE)
  
  lhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], lhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  rhs <- sapply(1:nrow(replacements), function(i) {
    out <- replaceSymbols(colnames(replacements), replacements[i, ], rhs)
    gsub(":", "_", out, fixed = TRUE)
  })
  
  if (is.null(trafo)) {
    trafo <- as.eqnvec(structure(lhs, names = lhs))
  } else if (is.list(trafo) & !reset) {
    trafo <- lapply(trafo, function(t) replaceSymbols(lhs, rhs, t))
  } else if (is.character(trafo) & !reset) {
    trafo <- replaceSymbols(lhs, rhs, trafo)
  } else if (is.list(trafo) & reset) {
    trafo <- lapply(trafo, function(t) {t[lhs] <- rhs; return(t)})
  } else if (is.character(trafo) & reset) {
    trafo[lhs] <- rhs
  }
  
  return(trafo)
  
  
}

paste_ <- function(...) paste(..., sep = "_")
