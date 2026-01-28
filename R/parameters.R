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
#' @param givenCQ Named character vector of conservation quantity replacements.
#'   Applies to `method = "implicit"` only. Names are states whose ODEs will be
#'   replaced by algebraic constraints. Values are expressions for those states.
#'   Example: `c(pAKT = "totAKT - AKT")` replaces pAKT's ODE with the constraint
#'   `pAKT - (totAKT - AKT) = 0`, effectively enforcing `pAKT = totAKT - AKT`.
#'   Default NULL means no conservation constraints are applied.
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
#' @export
P <- function(x = NULL,
              method = c("explicit", "implicit"), 
              parameters=NULL, 
              deriv = TRUE,
              deriv2 = FALSE,
              fixed = NULL,
              keep.root = TRUE, 
              positive = TRUE,
              givenCQ = NULL,
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
                   givenCQ = givenCQ,
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
                              givenCQ = givenCQ,
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
  if (is.null(modelname)) modelname <- "expl_parfn"
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
      outdir = getwd(),
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
      J_outer <- outPEval$jacobian[1, , , drop = FALSE]
      dim(J_outer) <- dim(J_outer)[2:3]
      dimnames(J_outer) <- dimnames(outPEval$jacobian)[2:3]
      
      if (has_upstream) {
        # Chain rule: J_inner = J_outer %*% dP
        # einsum("ij,jk->ik", J_outer, dP_sub) = simple matrix multiply
        nonfixed <- colnames(J_outer)
        dP_sub <- dP[nonfixed, , drop = FALSE]
        myderiv <- J_outer %*% dP_sub
      } else {
        myderiv <- J_outer
      }
    }
    
    if (deriv2 && !is.null(outPEval$hessian)) {
      H_outer <- outPEval$hessian[1, , , , drop = FALSE]
      dim(H_outer) <- dim(H_outer)[2:4]
      dimnames(H_outer) <- dimnames(outPEval$hessian)[2:4]
      
      if (has_upstream && !is.null(dP2)) {
        # Chain rule for Hessian
        nonfixed <- colnames(J_outer)
        dP_sub  <- dP[nonfixed, , drop = FALSE]
        dP2_sub <- dP2[nonfixed, , , drop = FALSE]
        
        n_out <- dim(H_outer)[1]  # number of output parameters
        n_m <- dim(dP_sub)[1]     # number of intermediate parameters  
        n_k <- dim(dP_sub)[2]     # number of outer parameters
        
        # term1: einsum("imn,mj,nk->ijk", H_outer, dP_sub, dP_sub)
        # = t(dP_sub) %*% H_outer[i,,] %*% dP_sub for each i
        # Use Kronecker: vec(term1[i,,]) = (dP_sub %x% dP_sub)^T %*% vec(H_outer[i,,])
        dPkron <- dP_sub %x% dP_sub  # (n_m * n_m) x (n_k * n_k)
        H_flat <- matrix(H_outer, nrow = n_out, ncol = n_m * n_m)
        term1_flat <- H_flat %*% dPkron  # n_out x (n_k * n_k)
        term1 <- array(term1_flat, dim = c(n_out, n_k, n_k))
        
        # term2: einsum("im,mjk->ijk", J_outer, dP2_sub)
        # = J_outer %*% dP2_sub_flat, then reshape
        # dP2_sub is [n_m, n_k, n_k], reshape to n_m x (n_k * n_k)
        dP2_flat <- matrix(dP2_sub, nrow = n_m, ncol = n_k * n_k)
        term2_flat <- J_outer %*% dP2_flat  # n_out x (n_k * n_k)
        term2 <- array(term2_flat, dim = c(n_out, n_k, n_k))
        
        myderiv2 <- term1 + term2
        dimnames(myderiv2) <- list(dimnames(H_outer)[[1]], colnames(dP_sub), colnames(dP_sub))
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
#' Constructs an implicit parameter transformation for steady-state problems.
#' Solves f(x, p) = 0 for dependent states x given parameters p using the
#' implicit function theorem for derivatives.
#'
#' @param x An `eqnvec` or `eqnlist` defining the residual equations f(x,p)=0.
#' @param parameters Character vector of outer parameter names.
#' @param deriv Logical. Compute Jacobian via implicit function theorem.
#' @param deriv2 Logical. Compute Hessian (implies deriv=TRUE).
#' @param fixed Character vector of compile-time fixed parameters.
#' @param keep.root Logical. Warm-start from previous solution.
#' @param positive Logical or named logical. Enforce positive solutions via log-transform.
#' @param givenCQ Named character vector of conservation quantity replacements.
#'   Names are states to replace, values are expressions for those states.
#'   Example: `c(pAKT = "totAKT - AKT")` replaces pAKT's ODE with the constraint
#'   `pAKT - (totAKT - AKT) = 0`. Default NULL (no replacements).
#' @param parlower,parupper Named numeric vectors of bounds for multistart.
#' @param nstart Integer. Number of multistart attempts.
#' @param optionsSolver List of options for nleqslv.
#' @param attach.input Logical. Include input parameters in output.
#' @param condition Character. Condition label.
#' @param compile Logical. Compile via funCpp.
#' @param modelname Character. Base name for compiled code.
#' @param verbose Logical. Print diagnostics.
#'
#' @return A [parfn] object.
#' @seealso [Pexpl], [P]
#' @importFrom nleqslv nleqslv
#' @export
Pimpl <- function(x, parameters = NULL, deriv = TRUE, deriv2 = FALSE, fixed = NULL,
                  keep.root = TRUE, positive = TRUE, givenCQ = NULL,
                  parlower = NULL, parupper = NULL, nstart = 10L,
                  optionsSolver = list(), attach.input = FALSE,
                  condition = NULL, compile = FALSE, modelname = NULL, verbose = FALSE) {
  
  forcing <- character(0)
  
  # === Process input ===
  if (inherits(x, "eqnlist")) {
    trafo <- as.eqnvec(x)
    
    # Detect forcing states (RHS effectively 0)
    # Handles patterns like "0", "1*(0)", "-1*(0)+1*(0)", etc.
    for (st in names(trafo)) {
      rhs <- as.character(trafo[st])
      if (.is_zero_rhs(rhs)) forcing <- c(forcing, st)
    }
    if (length(forcing)) {
      trafo <- trafo[setdiff(names(trafo), forcing)]
      message("Forcing states (treated as parameters): ", paste(forcing, collapse = ", "))
    }
  } else {
    trafo <- as.eqnvec(x)
  }
  
  # === Apply given conservation quantities ===
  if (!is.null(givenCQ)) {
    if (is.null(names(givenCQ)) || any(names(givenCQ) == ""))
      stop("givenCQ must be a named character vector")
    
    unknown <- setdiff(names(givenCQ), names(trafo))
    if (length(unknown))
      stop("givenCQ references unknown states: ", paste(unknown, collapse = ", "))
    
    trafo_vec <- setNames(as.character(trafo), names(trafo))
    for (st in names(givenCQ)) {
      # Replace ODE with constraint: state - expr = 0
      trafo_vec[st] <- paste0(st, " - (", givenCQ[st], ")")
    }
    trafo <- as.eqnvec(trafo_vec)
    
    if (verbose) {
      message("Applied conservation constraints:\n",
              paste(sprintf("  %s = %s", names(givenCQ), givenCQ), collapse = "\n"))
    }
  }
  
  # === Partition states ===
  states <- names(trafo)
  nonstates <- setdiff(getSymbols(trafo), states)
  # Forcing states are treated as parameters - add them to nonstates if not already there
  # (they should be there if they appear in RHS of other states, but might not be)
  nonstates <- unique(c(nonstates, forcing))
  if (is.null(parameters)) parameters <- nonstates
  
  indep_st <- intersect(states, parameters)
  dep_st <- setdiff(states, indep_st)
  if (!length(dep_st)) stop("No dependent states to solve.")
  if (!is.null(fixed) && any(fixed %in% dep_st))
    stop("Dependent states cannot be fixed: ", paste(intersect(fixed, dep_st), collapse = ", "))
  
  pos_mask <- .normalize_positive_mask(positive, dep_st)
  tfm <- .pos_transforms(pos_mask)
  
  # === Build evaluator ===
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition)) modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  all_params <- c(states, nonstates)
  diff_params <- setdiff(all_params, fixed)
  
  FEval <- suppressWarnings(CppODE::funCpp(
    as.eqnvec(trafo[dep_st]), variables = NULL, parameters = all_params,
    fixed = fixed, compile = compile, modelname = modelname, outdir = getwd(),
    verbose = verbose, convenient = FALSE, deriv = TRUE, deriv2 = deriv2
  ))
  
  # Solver options
  opts <- list(method = "Newton", global = "dbldog",
               control = list(ftol = 1e-6, xtol = 1e-6, maxit = 1000))
  for (k in c("method", "global")) if (!is.null(optionsSolver[[k]])) opts[[k]] <- optionsSolver[[k]]
  for (k in c("ftol", "xtol", "maxit")) {
    if (!is.null(optionsSolver[[k]])) opts$control[[k]] <- optionsSolver[[k]]
    if (!is.null(optionsSolver$control[[k]])) opts$control[[k]] <- optionsSolver$control[[k]]
  }
  
  guess_env <- new.env(parent = emptyenv())
  compile_fixed <- fixed
  dep_idx <- match(dep_st, diff_params)
  n_dep <- length(dep_st)
  
  # === Transformation closure ===
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, verbose = FALSE) {
    
    # Reset warm-start on new fit context
    token <- getOption(".dMod.fit_token")
    if (!is.null(token) && !identical(guess_env$token, token)) {
      guess_env$guess <- NULL; guess_env$token <- token
    }
    if (deriv2 && !deriv) { deriv <- TRUE; warning("deriv2 requires deriv=TRUE") }
    
    p <- c(pars, fixed)
    dP <- attr(pars, "deriv"); dP2 <- attr(pars, "deriv2")
    rt_fixed <- names(fixed)
    
    # Fill missing
    miss <- setdiff(c(dep_st, indep_st, nonstates), names(p))
    if (length(miss)) p[miss] <- 1
    
    x0 <- p[dep_st]
    if (keep.root && !is.null(guess_env$guess))
      x0[names(guess_env$guess)] <- guess_env$guess[names(x0)]
    
    # Pack helper
    pack <- function(xd) {
      xf <- setNames(numeric(length(states)), states)
      xf[dep_st] <- xd; xf[indep_st] <- p[indep_st]
      c(xf, p[nonstates])
    }
    
    # Residual and Jacobian
    fn <- function(z) {
      tryCatch(FEval(NULL, pack(tfm$to_x(z)), attach.input = FALSE, deriv = FALSE)$out[1,],
               error = function(e) rep(NaN, length(z)))
    }
    jac <- function(z) {
      E <- tryCatch(FEval(NULL, pack(tfm$to_x(z)), attach.input = FALSE, deriv = TRUE),
                    error = function(e) NULL)
      if (is.null(E)) return(NULL)
      Jx <- .slice_J(E$jacobian)[, dep_idx, drop = FALSE]
      if (any(pos_mask)) sweep(Jx, 2, ifelse(pos_mask, tfm$to_x(z), 1), `*`) else Jx
    }
    
    # Solve
    if (is.null(guess_env$guess) && !is.null(parlower)) {
      res <- .multistart_solve(fn, jac, x0, dep_st, pos_mask, parlower, parupper,
                               nstart, opts, tfm, verbose)
      x_star <- res$x; root <- res$root
    } else {
      root <- .quiet_nleqslv(fn, jac, tfm$to_z(x0), opts)
      x_star <- tfm$to_x(root$x)
    }
    
    if (!root$converged || any(!is.finite(x_star)))
      stop(sprintf("[Pimpl] No convergence (code %d: %s, |f|=%.2e)",
                   root$termcd, .nleqslv_msg(root$termcd), max(abs(root$fvec))))
    if (keep.root) guess_env$guess <- setNames(x_star, dep_st)
    if (verbose) cat(sprintf("[Pimpl] Converged: |f|=%.2e, iter=%d\n",
                             max(abs(root$fvec)), root$iter))
    
    # === Implicit Function Theorem ===
    E <- FEval(NULL, pack(x_star), attach.input = FALSE, deriv = TRUE, deriv2 = deriv2)
    J <- .slice_J(E$jacobian)
    
    # Partition indices
    dep_match <- match(dep_st, diff_params)
    ns_cols <- intersect(nonstates, diff_params)
    ind_cols <- intersect(indep_st, diff_params)
    ns_match <- if (length(ns_cols)) match(ns_cols, diff_params) else integer(0)
    ind_match <- if (length(ind_cols)) match(ind_cols, diff_params) else integer(0)
    
    dfdx <- J[, dep_match, drop = FALSE]
    dimnames(dfdx) <- list(dep_st, dep_st)
    
    # Singularity check
    rc <- tryCatch(rcond(dfdx), error = function(e) 0)
    if (rc < 1e-12) {
      sv <- svd(dfdx)$d
      stop(sprintf("Jacobian singular (rcond=%.2e). SVs: %s",
                   rc, paste(format(sv, digits = 3), collapse = ", ")))
    }
    inv_dfdx <- solve(dfdx)
    
    # First derivatives: dx/dp = -inv(dfdx) * dfdp
    Dx_ns <- if (length(ns_match)) -inv_dfdx %*% J[, ns_match, drop = FALSE]
    else matrix(0, n_dep, 0)
    Dx_ind <- if (length(ind_match)) -inv_dfdx %*% J[, ind_match, drop = FALSE]
    else matrix(0, n_dep, 0)
    if (ncol(Dx_ns)) { rownames(Dx_ns) <- dep_st; colnames(Dx_ns) <- ns_cols }
    if (ncol(Dx_ind)) { rownames(Dx_ind) <- dep_st; colnames(Dx_ind) <- ind_cols }
    
    # Build output
    out_st <- setNames(numeric(length(states)), states)
    out_st[dep_st] <- x_star; out_st[indep_st] <- p[indep_st]
    out_ns <- p[intersect(nonstates, parameters)]
    out_vec <- c(out_st, out_ns)
    
    cols <- setdiff(names(pars), c(compile_fixed, rt_fixed))
    n_cols <- length(cols)
    J_out <- matrix(0, length(out_vec), n_cols, dimnames = list(names(out_vec), cols))
    
    # Identity for pass-through
    pt <- intersect(setdiff(names(out_vec), dep_st), cols)
    if (length(pt)) J_out[cbind(pt, pt)] <- 1
    
    # IFT derivatives
    for (nm in intersect(ns_cols, cols)) J_out[dep_st, nm] <- Dx_ns[, nm]
    for (nm in intersect(ind_cols, cols)) J_out[dep_st, nm] <- Dx_ind[, nm]
    
    # === Second derivatives (BLAS-optimized, no einsum) ===
    # IFT: d²x/dpj dpk = -inv(dfdx) * [H_xx:Dx:Dx + H_xp:Dx + H_px:Dx + H_pp]
    H_out <- NULL
    if (deriv2) {
      H <- .slice_H(E$hessian)
      P_cols <- c(ns_cols, ind_cols)
      n_P <- length(P_cols)
      
      if (n_P > 0) {
        D_all <- cbind(Dx_ns, Dx_ind)
        colnames(D_all) <- P_cols
        
        # H_xx term: sum_mn H[i,m,n] * D[m,j] * D[n,k]
        # = (D ⊗ D)' vec(H[i]) for each i
        H_xx <- H[, dep_match, dep_match, drop = FALSE]
        Dkron <- D_all %x% D_all
        H_xx_flat <- matrix(H_xx, nrow = n_dep, ncol = n_dep * n_dep)
        rhs_flat <- H_xx_flat %*% Dkron
        
        # Add H_xp, H_px, H_pp terms for nonstates
        if (length(ns_match)) {
          n_ns <- length(ns_cols)
          H_xp <- H[, dep_match, ns_match, drop = FALSE]
          H_pp <- H[, ns_match, ns_match, drop = FALSE]
          
          # H_xp[i,m,j]*D[m,k] + H_xp[i,m,k]*D[m,j] + H_pp[i,j,k]
          # Process per output equation (loop is unavoidable for mixed indices)
          rhs_arr <- array(rhs_flat, dim = c(n_dep, n_P, n_P))
          
          for (i in seq_len(n_dep)) {
            H_xp_i <- H_xp[i, , , drop = TRUE]
            if (!is.matrix(H_xp_i)) H_xp_i <- matrix(H_xp_i, nrow = n_dep, ncol = n_ns)
            H_pp_i <- H_pp[i, , , drop = TRUE]
            if (!is.matrix(H_pp_i)) H_pp_i <- matrix(H_pp_i, nrow = n_ns, ncol = n_ns)
            
            # Cross terms: t(D_ns) %*% H_xp + t(H_xp) %*% D_ns
            cross <- crossprod(Dx_ns, H_xp_i) + crossprod(H_xp_i, Dx_ns)
            rhs_arr[i, 1:n_ns, 1:n_ns] <- rhs_arr[i, 1:n_ns, 1:n_ns] + cross + H_pp_i
          }
          rhs_flat <- matrix(rhs_arr, nrow = n_dep, ncol = n_P * n_P)
        }
        
        d2x_flat <- -inv_dfdx %*% rhs_flat
        d2x <- array(d2x_flat, dim = c(n_dep, n_P, n_P),
                     dimnames = list(dep_st, P_cols, P_cols))
        
        H_out <- array(0, c(length(out_vec), n_cols, n_cols),
                       dimnames = list(names(out_vec), cols, cols))
        act_P <- intersect(P_cols, cols)
        if (length(act_P)) H_out[dep_st, act_P, act_P] <- d2x[, act_P, act_P]
      }
    }
    
    # === Chain rule with upstream ===
    myderiv <- if (!is.null(dP)) J_out %*% dP[cols, , drop = FALSE] else J_out
    
    myderiv2 <- NULL
    if (deriv2 && !is.null(H_out)) {
      if (!is.null(dP) && !is.null(dP2)) {
        dP_sub <- dP[cols, , drop = FALSE]
        dP2_sub <- dP2[cols, , , drop = FALSE]
        n_out <- nrow(H_out)
        n_k <- ncol(dP_sub)
        
        # term1: H_out contracted with dP twice (Kronecker method)
        dPkron <- dP_sub %x% dP_sub
        H_flat <- matrix(H_out, nrow = n_out, ncol = n_cols * n_cols)
        term1_flat <- H_flat %*% dPkron
        term1 <- array(term1_flat, dim = c(n_out, n_k, n_k))
        
        # term2: J_out %*% dP2 (flatten and matmul)
        dP2_flat <- matrix(dP2_sub, nrow = n_cols, ncol = n_k * n_k)
        term2_flat <- myderiv %*% dP2_flat
        term2 <- array(term2_flat, dim = c(n_out, n_k, n_k))
        
        myderiv2 <- term1 + term2
        dimnames(myderiv2) <- list(names(out_vec), colnames(dP_sub), colnames(dP_sub))
      } else {
        myderiv2 <- H_out
      }
    }
    
    # Attach input
    if (attach.input) {
      extras <- setdiff(names(pars), names(out_vec))
      if (length(extras)) {
        out_vec <- c(out_vec, pars[extras])
        addJ <- matrix(0, length(extras), ncol(myderiv), dimnames = list(extras, colnames(myderiv)))
        hit <- intersect(extras, colnames(myderiv))
        if (length(hit)) addJ[cbind(hit, hit)] <- 1
        myderiv <- rbind(myderiv, addJ)
      }
    }
    
    as.parvec(out_vec, deriv = if (deriv) myderiv else FALSE,
              deriv2 = if (deriv2) myderiv2 else FALSE)
  }
  
  attr(p2p, "equations") <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- unique(parameters)
  attr(p2p, "modelname") <- modelname
  parfn(p2p, parameters, condition)
}


# ============================================================================
# Internal helpers
# ============================================================================

.slice_J <- function(J4d) {
  out <- J4d[1, , , drop = FALSE]
  dim(out) <- dim(J4d)[2:3]
  dimnames(out) <- dimnames(J4d)[2:3]
  out
}

.slice_H <- function(H5d) {
  out <- H5d[1, , , , drop = FALSE]
  dim(out) <- dim(H5d)[2:4]
  dimnames(out) <- dimnames(H5d)[2:4]
  out
}

.pos_transforms <- function(pos_mask) {
  if (!any(pos_mask)) return(list(to_z = identity, to_x = identity))
  list(
    to_z = function(x) { x[pos_mask] <- log(pmax(x[pos_mask], 1e-12)); x },
    to_x = function(z) { z[pos_mask] <- exp(pmin(pmax(z[pos_mask], -700), 700)); z }
  )
}

.nleqslv_msg <- function(code) {
  c("1" = "converged", "2" = "x-tol", "3" = "no step",
    "4" = "iter limit", "5" = "singular", "6" = "ill-cond")[as.character(code)] %||% "unknown"
}

.quiet_nleqslv <- function(fn, jac, start, opts) {
  r <- suppressWarnings(nleqslv::nleqslv(start, fn, jac,
                                         method = opts$method, global = opts$global, control = opts$control))
  r$converged <- (r$termcd == 1)
  r
}

.multistart_solve <- function(fn, jac, x0, dep_st, pos_mask, lower, upper,
                              nstart, opts, tfm, verbose) {
  lo <- setNames(ifelse(pos_mask, 1e-6, -1e6), dep_st)
  hi <- setNames(rep(1e6, length(dep_st)), dep_st)
  if (!is.null(lower)) lo[names(lower)[names(lower) %in% dep_st]] <- lower[names(lower) %in% dep_st]
  if (!is.null(upper)) hi[names(upper)[names(upper) %in% dep_st]] <- upper[names(upper) %in% dep_st]
  
  best <- list(fval = Inf)
  for (i in seq_len(nstart)) {
    z_try <- tfm$to_z(setNames(runif(length(dep_st), lo, hi), dep_st))
    root <- tryCatch(.quiet_nleqslv(fn, jac, z_try, opts), error = function(e) NULL)
    if (!is.null(root) && !any(is.na(root$x))) {
      fval <- max(abs(root$fvec))
      if (fval < best$fval) {
        best <- list(root = root, fval = fval, x = tfm$to_x(root$x))
        if (verbose) cat(sprintf("  [%d] |f|=%.2e\n", i, fval))
        if (fval < 1e-10) break
      }
    }
  }
  if (is.null(best$root)) {
    best$root <- .quiet_nleqslv(fn, jac, tfm$to_z(x0), opts)
    best$x <- tfm$to_x(best$root$x)
  }
  list(x = best$x, root = best$root)
}

.normalize_positive_mask <- function(positive, dep_st) {
  if (is.logical(positive) && length(positive) == 1L && is.null(names(positive)))
    return(setNames(rep(positive, length(dep_st)), dep_st))
  if (is.logical(positive) && !is.null(names(positive))) {
    mask <- setNames(rep(FALSE, length(dep_st)), dep_st)
    mask[intersect(names(positive), dep_st)] <- positive[intersect(names(positive), dep_st)]
    return(mask)
  }
  stop("'positive' must be TRUE, FALSE, or named logical vector")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.is_zero_rhs <- function(rhs) {
  # Check if RHS is effectively zero
  # Handles: "0", "1*(0)", "-1*(0)+1*(0)", etc.
  rhs <- gsub("\\s+", "", rhs)
  if (rhs %in% c("", "0", "0.0")) return(TRUE)
  
  # Remove all ±number*(0) terms
  rhs_clean <- gsub("[+-]?[0-9.]*\\*?\\(0\\)", "", rhs)
  rhs_clean <- gsub("^[+-]+", "", rhs_clean)  # Remove leading +/-
  
  # If nothing left, it was all zeros
  if (rhs_clean == "") return(TRUE)
  
  # Also check if no symbols remain (pure numeric that equals 0)
  if (length(getSymbols(rhs)) == 0) {
    val <- tryCatch(eval(parse(text = rhs)), error = function(e) NA)
    if (!is.na(val) && abs(val) < 1e-12) return(TRUE)
  }
  
  FALSE
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
    row <- row[, !colnames(row) %in% c("condition", "conditions"), drop = FALSE]
    
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