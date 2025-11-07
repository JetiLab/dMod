## Functions to generate parameter transformation ----


#' Generate a parameter transformation function
#'
#' Constructs a parameter transformation from symbolic expressions, either
#' **explicitly** (\eqn{x = g(p)}) or **implicitly** (\eqn{f(x, p) = 0}).
#' This function serves as a high-level wrapper for \link{Pexpl} and \link{Pimpl},
#' allowing for easy generation of condition-specific parameter mappings.
#'
#' @description
#' The function \code{P()} generates a parameter transformation function that maps
#' **outer parameters** to **inner parameters** used in a model.  
#' It can operate in two modes:
#'
#' - **Explicit mode** (`method = "explicit"`, see \link{Pexpl}):  
#'   Inner parameters are directly computed from symbolic expressions.
#'
#' - **Implicit mode** (`method = "implicit"`, see \link{Pimpl}):  
#'   Inner parameters are defined as steady states satisfying \eqn{f(x, p) = 0}.
#'
#' Both types can be combined with other transformations via arithmetic operators
#' (`+`, `*`, etc.) thanks to the \link{parfn} interface.
#'
#' @param trafo Object of class \code{eqnvec} or a named character vector,
#' or a list thereof. If a list is provided, \code{P()} is called on each
#' element and conditions are taken from the list names.
#' @param method Character, either \code{"explicit"} or \code{"implicit"}.
#' Determines whether to use \link{Pexpl} or \link{Pimpl}.
#' @param parameters Character vector of outer parameters.
#' @param deriv Logical. If \code{TRUE}, compute and attach the Jacobian of the
#' transformation as attribute \code{"deriv"}.
#' @param deriv2 Logical. If \code{TRUE}, compute and attach the Hessian as
#' attribute \code{"deriv2"}. Implies \code{deriv = TRUE}.
#' @param fixed Character vector of parameter names treated as fixed (no
#' derivatives returned with respect to them).
#' @param keep.root Logical. Applies to \code{method = "implicit"} only.
#' If \code{TRUE}, reuse the root from the previous call as a warm-start guess
#' for faster convergence.
#' @param positive Logical. Applies to \code{method = "implicit"} only.
#' If \code{TRUE}, the steady-state solver is performed in log-space
#' (\eqn{x = exp(z)}) to enforce positive solutions.
#' @param attach.input Logical. Attach input parameters to the output if they
#' are not overwritten by the transformation (identity mapping).
#' @param condition Character. Condition label for which the transformation is
#' generated. If \code{trafo} is a list, this is inferred from list names.
#' @param compile Logical. If \code{TRUE}, compile the transformation via
#' \link{funCpp} for improved performance.
#' @param modelname Character. Used when \code{compile = TRUE} to define the
#' base name of the generated C code file.
#' @param verbose Logical. Print compiler output and diagnostic messages to the
#' R console.
#'
#' @return
#' An object of class \link{parfn}, representing the parameter transformation.
#' The returned function \code{p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE)}
#' computes inner parameters and attaches derivatives as attributes
#' \code{"deriv"} and \code{"deriv2"} when requested.
#'
#' @seealso
#' \link{Pexpl} for explicit parameter transformations,
#' \link{Pimpl} for implicit (steady-state) transformations,
#' and \link{parfn} for details on combining transformations.
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @export
P <- function(trafo = NULL,
              method = c("explicit", "implicit"), 
              parameters=NULL, 
              deriv = TRUE,
              deriv2 = FALSE,
              fixed = NULL,
              keep.root = TRUE, 
              positive = TRUE, 
              attach.input = FALSE,
              condition = NULL, 
              compile = FALSE, 
              modelname = NULL, 
              verbose = FALSE) {
  
  if (is.null(trafo)) return()
  if (!is.list(trafo)) {
    trafo <- list(trafo)
    names(trafo) <- condition
  }
  
  method <- match.arg(method)
  
  Reduce("+", lapply(1:length(trafo), function(i) {
  
      switch(method, 
           explicit = Pexpl(trafo = as.eqnvec(trafo[[i]]), 
                            parameters = parameters,
                            deriv = deriv,
                            deriv2 = deriv2,
                            fixed = fixed,
                            attach.input = attach.input,
                            condition = names(trafo[i]),
                            compile = compile,
                            modelname = modelname, 
                            verbose = verbose),
           implicit = Pimpl(trafo = as.eqnvec(trafo[[i]]), 
                            parameters = parameters,
                            deriv = deriv,
                            deriv2 = deriv2,
                            fixed = fixed,
                            keep.root = keep.root,
                            positive = positive,
                            attach.input = attach.input,
                            condition = names(trafo[i]), 
                            compile = compile, 
                            modelname = modelname, 
                            verbose = verbose))
    
    
  }))
  
  
}

#' Parameter transformation (explicit)
#'
#' Constructs a parameter transformation function that maps outer parameters
#' to inner parameters according to symbolic expressions. 
#' 
#' @param trafo eqnvec or named character vector. Names correspond to the parameters being fed into
#' the model (the inner parameters). The elements of \code{trafo} are equations that express 
#' the inner parameters in terms of other parameters (the outer parameters).
#' @param parameters Character vector. Optional. If given, the generated parameter
#' transformation returns values for each element in \code{parameters}. If elements of
#' \code{parameters} are not in \code{names(trafo)} the identity transformation is assumed.
#' @param deriv Logical, if \code{TRUE} the Jacobian of the transformation is precomputed 
#' symbolically and returned as attribute \code{"deriv"}.
#' @param deriv2 Logical, if \code{TRUE} the Hessian of the transformation is also precomputed 
#' symbolically and returned as attribute \code{"deriv2"}. Implies \code{deriv = TRUE}.
#' @param fixed Character vector of parameter names treated as fixed (no derivatives returned 
#' with respect to them).
#' @param attach.input Attach those incoming parameters to the output which are not overwritten by
#' the parameter transformation. 
#' @param compile Logical, compile the function (see \link{funCpp}).
#' @param condition Character, the condition for which the transformation is generated.
#' @param modelname Character, used if \code{compile = TRUE}, sets a fixed filename for the
#' C file.
#' @param verbose Print compiler output to the R console.
#' 
#' @return A function \code{p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE)} representing the 
#' parameter transformation. Here, \code{p} is a named numeric vector with the values of the 
#' outer parameters, \code{fixed} is a named numeric vector with values of the outer parameters 
#' being considered as fixed (no derivatives returned), and \code{deriv}/\code{deriv2} determine 
#' whether the Jacobian and Hessian of the transformation are attached as attributes 
#' \code{"deriv"} and \code{"deriv2"}.
#' 
#' @seealso \link{Pimpl} for implicit parameter transformations.
#' @export
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
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
  symbols <- getSymbols(trafo)
  
  if (is.null(parameters)) {
    parameters <- symbols
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
  PEval <- CppODE::funCpp(
    x          = trafo,
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
  
  # ---------------------------------------------------------------------------
  # Define returned parameter transformation function
  # ---------------------------------------------------------------------------
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, env = parent.frame()) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE; enabling deriv = TRUE.")
      deriv <- TRUE
    }
    
    # Evaluate inner parameters
    outPEval <- PEval(NULL, pars, attach.input = attach.input, deriv = deriv, deriv2 = deriv2)
    
    pinner <- setNames(as.numeric(outPEval$out), colnames(outPEval$out))
    # Sanity check
    if (any(is.nan(pinner))) {
      stop(
        paste0(
          "The following inner parameter(s) evaluate to NaN:\n\t",
          paste0(names(pinner)[is.nan(pinner)], collapse = "\n\t"),
          ".\nLikely cause: division by zero or missing inputs."
        )
      )
    }
    
    # Identify fixed vs. active parameters
    fixed_now <- if (is.null(fixed)) character(0) else names(fixed)
    nonfixed  <- setdiff(attr(PEval, "parameter"), fixed_now)
    
    # ----------------- Apply chain rules -----------------
    myderiv  <- NULL
    myderiv2 <- NULL
    
    if (deriv2) {
      
      # ----- First-order derivative -----
      J_outer <- outPEval$jacobian[,,1]
      if (!is.null(J_outer)) {
        J_outer <- J_outer[, nonfixed, drop = FALSE]
        
        dP  <- attr(pars, "deriv",  exact = TRUE)
        dP2 <- attr(pars, "deriv2", exact = TRUE)
        
        if (!is.null(dP) && !is.null(dP2)) {
          # Parameter transformation active
          dP <- dP[nonfixed, , drop = FALSE]
          myderiv <- einsum::einsum("ij,jk->ik", J_outer, dP)
        } else {
          # No transformation
          myderiv <- J_outer
        }
      }
      
      # ----- Second-order derivative -----
      H_outer <- outPEval$hessian[,,,1]
      if (!is.null(H_outer)) {
        H_outer <- H_outer[, nonfixed, nonfixed, drop = FALSE]
        
        dP  <- attr(pars, "deriv",  exact = TRUE)
        dP2 <- attr(pars, "deriv2", exact = TRUE)
        
        if (!is.null(dP) && !is.null(dP2)) {
          dP  <- dP [nonfixed, , drop = FALSE]
          dP2 <- dP2[nonfixed, , , drop = FALSE]
          
          # term1: contraction of outer Hessian with two Jacobians
          term1 <- einsum::einsum("imn,mj,nk->ijk", H_outer, dP, dP)
          # term2: contraction of outer Jacobian with second parameter derivatives
          J_outer <- attr(pinner, "jacobian")[, nonfixed, 1, drop = FALSE]
          term2 <- einsum::einsum("im,mjk->ijk", J_outer, dP2)
          # combine both
          myderiv2 <- term1 + term2
          
        } else {
          # no parameter transformation
          myderiv2 <- H_outer
        }
      }
      
    } else if (deriv) {
      
      # ----- First-order derivative only -----
      J_outer <- outPEval$jacobian[,,1]
      if (!is.null(J_outer)) {
        J_outer <- J_outer[, nonfixed, drop = FALSE]
        
        dP  <- attr(pars, "deriv",  exact = TRUE)
        dP2 <- attr(pars, "deriv2", exact = TRUE)
        
        if (!is.null(dP) && !is.null(dP2)) {
          dP <- dP[nonfixed, , drop = FALSE]
          myderiv <- einsum::einsum("ij,jk->ik", J_outer, dP)
        } else {
          myderiv <- J_outer
        }
      }
    }
    
    
    # -------------------------------------------------------------------------
    # Assemble result and return
    # -------------------------------------------------------------------------
    pinner <- as.parvec(pinner, 
                        deriv = if(deriv) myderiv else FALSE, 
                        deriv2 = if(deriv2) myderiv2 else FALSE)
    
    if (attach.input && !all(names(pars) %in% names(pinner))) {
      pinner <- c(pinner, 
                  as.parvec(pars[setdiff(names(pars), names(pinner))],
                            deriv = if(deriv) NULL else FALSE,
                            deriv2 = if(deriv2) NULL else FALSE))
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
#' Constructs an implicit parameter transformation that solves a system of
#' steady-state equations \eqn{f(x, p) = 0}. The solution \eqn{x(p)} defines
#' the mapping from outer parameters \eqn{p} to inner steady-state parameters \eqn{x}.
#'
#' @description
#' This function computes steady-state mappings using a root solver and attaches
#' first- and optionally second-order derivatives via the implicit function theorem.
#' If \code{positive = TRUE}, dependent states are solved in log-space (\eqn{x=\exp(z)})
#' to enforce positivity, but derivatives are returned in \eqn{x}-space (no extra
#' \eqn{\mathrm{diag}(x)} factors are applied).
#'
#' The input \code{trafo} may be an \code{eqnvec} (named character vector) specifying
#' residuals \eqn{f(x,p)=0} directly, or an \code{eqnlist} (stoichiometric representation).
#' For \code{eqnlist}, residuals are built via \code{as.eqnvec()}, conserved quantities
#' are detected via \code{conservedQuantities()}, and each conservation law is turned
#' into an algebraic constraint using an automatically introduced "total" parameter.
#' This removes rank deficiency without requiring users to choose independent states.
#'
#' @param trafo Residual system as \code{eqnvec}/named character vector
#'   (names are states), or an \code{eqnlist} from which the residuals are derived.
#' @param parameters Character vector of outer parameter names (optional).
#'   For \code{eqnlist} input, if omitted, this is inferred as
#'   \{model/rate parameters\} ∪ \{totals\}. If provided, it is validated.
#'   If \code{parameters} contains some of the state names, those states are treated
#'   as independent (not solved).
#' @param deriv Logical; if \code{TRUE}, attach the Jacobian as attribute \code{"deriv"}.
#' @param deriv2 Logical; if \code{TRUE}, also attach the Hessian as attribute \code{"deriv2"}.
#'   Implies \code{deriv = TRUE}.
#' @param fixed Character vector of parameter names to be treated as fixed
#'   (no derivatives returned w.r.t. them).
#' @param keep.root Logical; reuse the last root as starting guess.
#' @param positive Logical; enforce positivity of solved dependent states by solving
#'   in log-space internally. Derivatives are still taken in \eqn{x}-space.
#' @param attach.input Logical; attach incoming parameters to the output if they are
#'   not overwritten by the implicit transformation (identity mapping).
#' @param condition Character; condition label (optional).
#' @param compile Logical; compile generated functions via \link[CppODE]{funCpp}.
#' @param modelname Character; base name for compiled code files when \code{compile=TRUE}.
#' @param verbose Logical; print compiler output and diagnostics.
#'
#' @return
#' A function \code{p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE)} returning the inner
#' parameter vector with attached attributes \code{"deriv"} and optionally \code{"deriv2"}.
#' The set of expected outer parameters is stored in \code{attr(p2p, "parameters")}.
#'
#' @seealso \link{Pexpl}, \link[rootSolve]{multiroot}, \link{as.eqnvec}, \link{conservedQuantities}
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @importFrom rootSolve multiroot
#' @importFrom abind abind
#' @export
Pimpl <- function(trafo,
                  parameters   = NULL,
                  deriv        = TRUE,
                  deriv2       = FALSE,
                  fixed        = NULL,
                  keep.root    = TRUE,
                  positive     = TRUE,
                  attach.input = FALSE,
                  condition    = NULL,
                  compile      = FALSE,
                  modelname    = NULL,
                  verbose      = FALSE) {
  
  # Normalize input: accept eqnlist; insert conservation totals if needed
  totals <- character(0)
  if (inherits(trafo, "eqnlist")) {
    eq    <- trafo
    trafo <- as.eqnvec(eq)
    
    cq_df <- conservedQuantities(eq$smatrix)
    n_cq  <- if (!is.null(attr(cq_df, "n"))) attr(cq_df, "n") else 0L
    
    states_detected <- names(trafo)
    symbols_all     <- getSymbols(trafo)
    nonstates_auto  <- setdiff(symbols_all, states_detected)
    
    if (n_cq > 0L) {
      if (!is.null(parameters)) {
        cand_totals <- setdiff(parameters, nonstates_auto)
        totals <- if (length(cand_totals) == n_cq) cand_totals else paste0("total", seq_len(n_cq))
      } else {
        totals <- paste0("total", seq_len(n_cq))
      }
      # Replace some state residuals by conservation constraints "(sum) - total_i = 0"
      repl_states <- tail(states_detected, n_cq)
      for (i in seq_len(n_cq)) {
        cons_expr <- cq_df[i, 1]
        trafo[repl_states[i]] <- paste0("(", cons_expr, ") - ", totals[i])
      }
      # Recompute symbol sets and ensure totals are not counted as nonstates
      states_detected <- names(trafo)
      symbols_all     <- getSymbols(trafo)
      nonstates_auto  <- setdiff(symbols_all, c(states_detected, totals))
    }
    
    required_outer <- unique(c(nonstates_auto, totals))
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
  }
  
  # Partition symbols: states vs outer parameters
  states    <- names(trafo)
  nonstates <- setdiff(getSymbols(trafo), states)
  
  # Allow some states to be independent by listing them in `parameters`
  indep_st <- intersect(states, if (is.null(parameters)) character(0) else parameters)
  dep_st   <- setdiff(states, indep_st)
  if (length(dep_st) == 0L) {
    stop("Pimpl: no dependent states to solve; 'parameters' covers all states.")
  }
  
  # Build compiled evaluator for the dependent residuals only
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition))
    modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  FEval <- CppODE::funCpp(
    x          = as.eqnvec(trafo[dep_st]),
    variables  = NULL,
    parameters = c(states, nonstates),                # fixed order: states, then nonstates
    fixed      = NULL,
    compile    = compile,
    modelname  = modelname,
    verbose    = verbose,
    convenient = FALSE,
    deriv      = isTRUE(deriv) || isTRUE(deriv2),     # need Jacobian if either requested
    deriv2     = isTRUE(deriv2)                       # Hessian only if requested
  )
  
  # Warm-start state
  guess_env <- new.env(parent = emptyenv())
  guess_env$guess <- NULL
  
  # Returned transformation closure
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, env = parent.frame()) {
    
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE; enabling deriv = TRUE.")
      deriv <- TRUE
    }
    
    # Upstream derivative attributes (for chained transformations)
    dP  <- attr(pars, "deriv",  exact = TRUE)
    dP2 <- attr(pars, "deriv2", exact = TRUE)
    
    # Merge fixed values
    p <- pars
    if (!is.null(fixed)) {
      if (any(names(fixed) %in% names(p))) p[names(fixed)] <- fixed else p <- c(p, fixed)
    }
    
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
    
    # Solve f_dep(x_dep,·)=0; positivity via z=log(x_dep) if requested
    if (positive) {
      f_z <- function(z, .) {
        x_dep <- exp(z)
        FEval(NULL, p = pack_full(x_dep),
              attach.input = FALSE, deriv = FALSE, deriv2 = FALSE)$out[1, ]
      }
      z0         <- log(pmax(x_dep0, 1e-8))
      root       <- rootSolve::multiroot(f_z, start = z0, parms = NULL, positive = FALSE)
      x_dep_star <- exp(root$root)
    } else {
      f_x <- function(x_dep, .) {
        FEval(NULL, p = pack_full(x_dep),
              attach.input = FALSE, deriv = FALSE, deriv2 = FALSE)$out[1, ]
      }
      root       <- rootSolve::multiroot(f_x, start = x_dep0, parms = NULL, positive = FALSE)
      x_dep_star <- root$root
    }
    
    if (keep.root) guess_env$guess <- setNames(x_dep_star, dep_st)
    
    # Evaluate Jacobian/Hessian of the dependent residuals at the solution
    need_hess <- isTRUE(deriv2)
    E <- FEval(NULL, p = pack_full(x_dep_star),
               attach.input = FALSE, deriv = TRUE, deriv2 = need_hess)
    
    # FEval Jacobian is w.r.t. c(states, nonstates)
    J <- E$jacobian[,,1, drop = TRUE]
    all_cols <- c(states, nonstates)
    
    # Block partition df/d(·)
    dfdx_dep   <- J[, match(dep_st,   all_cols), drop = FALSE]
    dfdx_indep <- if (length(indep_st)) J[, match(indep_st, all_cols), drop = FALSE]
    else matrix(0, nrow = nrow(dfdx_dep), ncol = 0, dimnames = list(dep_st, character(0)))
    dfdp_ns    <- J[, match(nonstates, all_cols), drop = FALSE]
    
    rownames(dfdx_dep) <- dep_st
    colnames(dfdx_dep) <- dep_st
    if (ncol(dfdx_indep)) {
      rownames(dfdx_indep) <- dep_st
      colnames(dfdx_indep) <- indep_st
    }
    rownames(dfdp_ns) <- dep_st
    colnames(dfdp_ns) <- nonstates
    
    # First-order IFT in x-space
    inv_fxd <- solve(dfdx_dep)
    Dx_ns   <- - inv_fxd %*% dfdp_ns
    Dx_ind  <- if (ncol(dfdx_indep)) - inv_fxd %*% dfdx_indep
    else matrix(0, nrow = nrow(dfdx_dep), ncol = 0, dimnames = list(dep_st, character(0)))
    
    # Prepare the base output vector (solved states + pass-through)
    x_full_star <- setNames(numeric(length(states)), states)
    x_full_star[dep_st]   <- x_dep_star
    x_full_star[indep_st] <- x_ind
    out_vec <- c(x_full_star, r_ns)
    
    # Build Jacobian over incoming parameter basis (excluding fixed)
    cols_nonfixed <- setdiff(names(pars), names(fixed))
    J_outer <- matrix(0, nrow = length(out_vec), ncol = length(cols_nonfixed),
                      dimnames = list(names(out_vec), cols_nonfixed))
    
    # Identity for pass-through outputs already in out_vec
    passthrough <- setdiff(names(out_vec), dep_st)
    if (length(passthrough)) {
      idx <- intersect(passthrough, cols_nonfixed)
      if (length(idx)) J_outer[idx, idx] <- diag(1, length(idx))
    }
    
    # Fill derivatives for dependent states
    active_ns  <- intersect(nonstates, cols_nonfixed)
    active_ind <- intersect(indep_st,  cols_nonfixed)
    if (length(active_ns))  J_outer[dep_st, active_ns]  <- Dx_ns[,  active_ns,  drop = FALSE]
    if (length(active_ind)) J_outer[dep_st, active_ind] <- Dx_ind[, active_ind, drop = FALSE]
    
    # Second-order IFT (x-space)
    H_outer <- NULL
    if (deriv2) {
      H <- E$hessian[,,,1, drop = TRUE]
      
      idx_xdep <- match(dep_st,   all_cols)
      idx_xind <- match(indep_st, all_cols)
      idx_ns   <- match(nonstates, all_cols)
      
      H_xx   <- H[, idx_xdep, idx_xdep, drop = FALSE]
      H_xp   <- H[, idx_xdep, idx_ns,   drop = FALSE]
      H_pp   <- H[, idx_ns,   idx_ns,   drop = FALSE]
      H_xind <- if (length(indep_st)) H[, idx_xdep, idx_xind, drop = FALSE] else NULL
      H_indp <- if (length(indep_st)) H[, idx_xind, idx_ns,   drop = FALSE] else NULL
      
      P_cols <- c(nonstates, indep_st)
      D_depP <- cbind(Dx_ns, Dx_ind)
      colnames(D_depP) <- P_cols
      
      rhs <- array(0, dim = c(length(dep_st), length(P_cols), length(P_cols)),
                   dimnames = list(dep_st, P_cols, P_cols))
      
      rhs <- rhs + einsum::einsum("iab,aj,bk->ijk", H_xx, D_depP, D_depP)
      
      if (length(nonstates)) {
        D_ns <- Dx_ns
        rhs[, nonstates, nonstates] <- rhs[, nonstates, nonstates] +
          H_pp +
          einsum::einsum("iak,aj->ijk", H_xp, D_ns) +
          einsum::einsum("iaj,ak->ijk", H_xp, D_ns)
      }
      
      if (length(indep_st)) {
        D_ind <- Dx_ind
        rhs[, nonstates, indep_st] <- rhs[, nonstates, indep_st] +
          einsum::einsum("iak,aj->ijk", H_xp, D_ind) +
          einsum::einsum("iab,aj,bk->ijk", H_xx, Dx_ns, D_ind)
        
        rhs[, indep_st, nonstates] <- rhs[, indep_st, nonstates] +
          einsum::einsum("iak,aj->ijk", H_indp, D_ns) +
          einsum::einsum("iab,aj,bk->ijk", H_xx, D_ind, Dx_ns)
        
        rhs[, indep_st, indep_st] <- rhs[, indep_st, indep_st] +
          einsum::einsum("iab,aj,bk->ijk", H_xx, D_ind, D_ind) +
          einsum::einsum("iak,aj->ijk", H_xind, D_ind) +
          einsum::einsum("iaj,ak->ijk", H_xind, D_ind)
      }
      
      rhs_mat <- matrix(rhs, nrow = nrow(dfdx_dep), ncol = length(P_cols) * length(P_cols))
      sol_mat <- - solve(dfdx_dep, rhs_mat)
      d2x_dep <- array(sol_mat, dim = c(length(dep_st), length(P_cols), length(P_cols)),
                       dimnames = list(dep_st, P_cols, P_cols))
      
      # Embed into full Hessian over current outputs and incoming cols
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
        hit  <- intersect(extras, colnames(J_outer))            # only non-fixed appear as columns
        if (length(hit)) addJ[cbind(hit, hit)] <- 1
        J_outer <- rbind(J_outer, addJ)
        
        # Expand Hessian with zero rows for the new outputs
        if (deriv2) {
          addH <- array(0, dim = c(length(extras), ncol(J_outer), ncol(J_outer)),
                        dimnames = list(extras, colnames(J_outer), colnames(J_outer)))
          H_outer <- abind::abind(H_outer, addH, along = 1)     # keep [out x in x in]
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
        termB   <- einsum::einsum("im,mjk->ijk",   myderiv, dP2_use)  # use composed J
        myderiv2 <- termA + termB
      } else {
        myderiv2 <- H_outer
      }
    }
    
    # Build final parvec
    res <- as.parvec(out_vec,
                     deriv  = if (deriv)  myderiv  else FALSE,
                     deriv2 = if (deriv2) myderiv2 else FALSE)
    
    # No extra append here; we already expanded out_vec/J/H above when attach.input=TRUE
    res
  }
  
  # Metadata
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- unique(parameters)
  attr(p2p, "modelname")  <- modelname
  
  parfn(p2p, parameters, condition)
}






## Functions to simplify the creation of parameter transformations ----

#' Define parameter transformations by \code{define()}, \code{branch()} and \code{insert()}
#' 
#' @param trafo named character vector of parametric expressions or object 
#' of class \code{eqnvec}
#' @param expr character of the form \code{"lhs ~ rhs"} where both \code{lhs}
#' and \code{rhs} can contain a number of symbols for which vaues are passed
#' by the \code{...} argument
#' @param  conditionMatch optional character, Use as regular expression to apply the reparameterization only to conditions containing conditionMatch
#' @param ... used to pass values for symbols as named arguments
#' @return object of the same class as trafo or list thereof, if \code{branch()} has been
#' used.
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
#' @param table table of covariates as data frame. Rownames are used as unique identifier,
#' usually called "conditions", and columns represent covariates associated with these conditions.
#' @param conditions character vector with condition names. Overwrites the rownames of table.
branch <- function(trafo, table = NULL, conditions = rownames(table)) {
  
  if (is.null(table) & is.null(conditions)) 
    return(trafo)
  
  if (is.null(conditions)) conditions <- paste0("C", 1:nrow(table))
  if (is.null(table)) table <- data.frame(condition = conditions, row.names = conditions)
  
  out <- lapply(conditions, function(x) trafo)
  names(out) <- conditions
  
  rownames(table) <- conditions
  attr(out, "tree") <- table
  
  return(out)  
  
}


#' Reparameterization
#' 
#' @param expr character of the form \code{"lhs ~ rhs"} where \code{rhs}
#' reparameterizes \code{lhs}. Both \code{lhs} and \code{rhs}
#' can contain a number of symbols whose values need to be passed by the \code{...} argument.
#' @param trafo character or equation vector or list thereof. The object where the replacement takes place in
#' @param ... pass symbols as named arguments
#' @param reset logical. If true, the trafo element corresponding to lhs is reset according to rhs. 
#' If false, lhs wherever it occurs in the rhs of trafo is replaced by rhs of the formula.
#' @return an equation vector with the reparameterization.
#' @details Left and right-hand side of \code{expr} are searched for symbols. If separated by
#' "_", symbols are recognized as such, e.g. in \code{Delta_x} where the symbols are 
#' "Delta" and "x". Each symbol for which values (character or numbers) are passed by the
#' \code{...} argument is replaced.
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