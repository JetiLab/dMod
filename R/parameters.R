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
#' Constructs an implicit parameter transformation that solves a set of
#' steady-state equations of the form \eqn{f(x, p) = 0}. The solution
#' \eqn{x(p)} defines the implicit mapping from outer parameters \eqn{p}
#' to inner steady-state parameters \eqn{x}.
#'
#' @description
#' This function computes steady-state mappings (implicit transformations)
#' using a numerical root solver. It supports both first- and second-order
#' derivatives obtained symbolically via \link{funCpp} and combined with the
#' implicit function theorem. If \code{positive = TRUE}, the dependent variables
#' are internally reparameterized as \eqn{x = exp(z)} to guarantee positive
#' steady-state solutions.
#'
#' @param trafo Named character vector defining the equations to be set to zero
#' (\eqn{f(x, p) = 0}). Names correspond to dependent variables \eqn{x}.
#' @param parameters Character vector. Optional. If provided, determines which
#' outer parameters are expected as inputs.
#' @param deriv Logical. If \code{TRUE}, the Jacobian of the implicit transformation
#' is computed and attached as attribute \code{"deriv"}.
#' @param deriv2 Logical. If \code{TRUE}, the Hessian of the transformation is
#' also computed and attached as attribute \code{"deriv2"}. Implies
#' \code{deriv = TRUE}.
#' @param fixed Character vector of parameter names treated as fixed (no
#' derivatives are returned with respect to them).
#' @param keep.root Logical. Whether to reuse the root from the previous call as
#' initial guess for the next evaluation (warm start). Can improve convergence
#' for related parameter sets.
#' @param positive Logical. If \code{TRUE}, the dependent variables are reparameterized
#' as \eqn{x = exp(z)} internally to enforce positivity of steady-state solutions.
#' @param attach.input Logical. Attach those incoming parameters to the output
#' that are not overwritten by the implicit transformation (identity mapping).
#' @param condition Character. Condition label for the transformation (optional).
#' @param compile Logical. If \code{TRUE}, compile the function via \link{funCpp}
#' for improved performance.
#' @param modelname Character. Base name for the generated C file(s) if compiled.
#' @param verbose Logical. Print compiler output and additional diagnostic
#' messages to the R console.
#'
#' @return
#' A function \code{p2p(p, fixed = NULL, deriv = TRUE, deriv2 = FALSE)} representing
#' the implicit parameter transformation. It returns a named numeric vector of inner
#' parameters with attached Jacobian (\code{"deriv"}) and optionally Hessian
#' (\code{"deriv2"}) attributes.
#'
#' @details
#' The transformation is defined implicitly by the steady-state condition
#' \eqn{f(x, p) = 0}. Given the outer parameters \eqn{p}, the steady state
#' \eqn{x(p)} is computed by \link[rootSolve]{multiroot}.  
#' 
#' First- and second-order derivatives are obtained using the implicit function
#' theorem:
#' \deqn{dx/dp = - (df/dx)^{-1} (df/dp)}
#' \deqn{d^2x/dp^2 = - (df/dx)^{-1}(f_{xx}[dx/dp,dx/dp] + f_{xp}dx/dp + f_{xp}^Tdx/dp + f_{pp})}
#' 
#' If \code{positive = TRUE}, the root is solved in log-space (\eqn{z = log(x)}),
#' ensuring \eqn{x > 0}. The corresponding derivatives are automatically
#' transformed back to \eqn{x}-space using the chain rule:
#' \eqn{dx/dp = diag(x) \, dz/dp},  
#' \eqn{d^2x/dp^2 = diag(x)(dz/dp \, dz/dp^T + d^2z/dp^2)}.
#'
#' @seealso
#' \link{Pexpl} for explicit parameter transformations,
#' \link[rootSolve]{multiroot} for steady-state computation,
#' and \link{funCpp} for the underlying symbolic compilation.
#'
#' @examples
#' ########################################################################
#' ## Example 1: Steady-state trafo
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#'
#' p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
#' P.steadyState(p.outerValues)
#'
#' ########################################################################
#' ## Example 2: Steady-state trafo combined with log-transform
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B",
#'        B = "k1*A - k2*B")
#' P.steadyState <- Pimpl(f, "A")
#'
#' logtrafo <- c(k1 = "exp(logk1)", k2 = "exp(logk2)",
#'               A = "exp(logA)", B = "exp(logB)")
#' P.log <- P(logtrafo)
#'
#' p.outerValue <- c(logk1 = 1, logk2 = -1, logA = 0, logB = 0)
#' (P.log)(p.outerValue)
#' (P.steadyState * P.log)(p.outerValue)
#'
#' ########################################################################
#' ## Example 3: Steady-states with conserved quantities
#' ########################################################################
#' f <- c(A = "-k1*A + k2*B", B = "k1*A - k2*B")
#' replacement <- c(B = "A + B - total")
#' f[names(replacement)] <- replacement
#'
#' pSS <- Pimpl(f, "total")
#' pSS(c(k1 = 1, k2 = 2, A = 5, B = 5, total = 3))
#'
#' @importFrom CppODE funCpp
#' @importFrom einsum einsum
#' @importFrom rootSolve multiroot
#' @export
Pimpl <- function(trafo,
                  parameters = NULL,
                  deriv = TRUE,
                  deriv2 = FALSE,
                  fixed = NULL,
                  keep.root  = TRUE,
                  positive   = TRUE,
                  attach.input = FALSE,
                  condition  = NULL,
                  compile    = FALSE,
                  modelname  = NULL,
                  verbose    = FALSE) {
  
  # ---------------------------------------------------------------------------
  # Identify dependent ("states") and independent ("nonstates") variables
  # ---------------------------------------------------------------------------
  states    <- names(trafo)
  nonstates <- getSymbols(trafo, exclude = states)
  
  parameters <- c(states, nonstates)
  
  # Set model name for compiled code
  if (is.null(modelname)) modelname <- "impl_parfn"
  if (!is.null(condition))
    modelname <- paste(modelname, sanitizeConditions(condition), sep = "_")
  
  # ---------------------------------------------------------------------------
  # Build compiled evaluator for f(x,p)
  # This object can compute f, its Jacobian (df/dx, df/dp) and Hessian.
  # ---------------------------------------------------------------------------
  FEval <- CppODE::funCpp(
    x          = trafo,
    variables  = NULL,
    parameters = parameters,  # parameter order = (x, p)
    fixed      = fixed,
    compile    = compile,
    modelname  = modelname,
    verbose    = verbose,
    convenient = FALSE,
    deriv      = deriv,
    deriv2     = deriv2
  )
  
  # Create an environment to store the last successful root (for warm-start)
  guess_env <- new.env(parent = emptyenv())
  guess_env$guess <- NULL
  
  # ---------------------------------------------------------------------------
  # Returned closure: evaluates the implicit transformation for given parameters
  # ---------------------------------------------------------------------------
  p2p <- function(pars, fixed = NULL, deriv = TRUE, deriv2 = FALSE, env = parent.frame()) {
    
    # Ensure consistency between deriv and deriv2 flags
    if (deriv2 && !deriv) {
      warning("deriv2 = TRUE requires deriv = TRUE; enabling deriv = TRUE.")
      deriv <- TRUE
    }
    
    # Extract possible upstream derivatives if this transformation is nested
    dP  <- attr(pars, "deriv",  exact = TRUE)
    dP2 <- attr(pars, "deriv2", exact = TRUE)
    
    # -------------------------------------------------------------------------
    # Merge in fixed parameters (overwrite if present)
    # -------------------------------------------------------------------------
    p <- pars
    if (!is.null(fixed)) {
      if (any(names(fixed) %in% names(p)))
        p[names(fixed)] <- fixed
      else
        p <- c(p, fixed)
    }
    
    # Ensure all dependent variables exist in the input (use 1 as fallback)
    if (!all(states %in% names(p))) p[setdiff(states, names(p))] <- 1
    x0 <- p[states]                         # current guess for dependent vars
    r  <- p[setdiff(names(p), states)]      # remaining (independent) params
    
    # Use previous root as starting value if keep.root = TRUE
    if (keep.root && !is.null(guess_env$guess))
      x0[names(x0) %in% names(guess_env$guess)] <-
      guess_env$guess[names(x0) %in% names(guess_env$guess)]
    
    # -------------------------------------------------------------------------
    # Solve f(x,p)=0  -- optionally in log-space for positivity
    # -------------------------------------------------------------------------
    if (positive) {
      # Define the transformed system f(exp(z), p) = 0  → solve for z
      f_z <- function(z, parms) {
        x <- exp(z)
        FEval(NULL, p = c(x, parms),
              attach.input = FALSE, deriv = FALSE, deriv2 = FALSE)$out[1, ]
      }
      
      # Log-transform initial guess (protect against zeros)
      z0 <- log(pmax(x0, 1e-6))
      
      # Call multiroot in unconstrained z-space
      root <- rootSolve::multiroot(f_z, start = z0, parms = r, positive = FALSE)
      
      # Extract results: steady-state z and corresponding x = exp(z)
      z_star <- root$root
      x_star <- exp(z_star)
      
    } else {
      # Direct solve in x-space (no positivity enforced)
      f_x <- function(x, parms) {
        FEval(NULL, p = c(x, parms),
              attach.input = FALSE, deriv = FALSE, deriv2 = FALSE)$out[1, ]
      }
      
      root <- rootSolve::multiroot(f_x, start = x0, parms = r, positive = FALSE)
      x_star <- root$root
      z_star <- log(pmax(x_star, 1e-6))  # keep z for later chain rule
    }
    
    # Store root as warm-start guess for next call (if enabled)
    if (keep.root) guess_env$guess <- x_star
    
    # -------------------------------------------------------------------------
    # Evaluate Jacobian/Hessian of f(x,p) at the solution (symbolically)
    # -------------------------------------------------------------------------
    need_hess <- isTRUE(deriv2)
    E <- FEval(NULL, p = c(x_star, r),
               attach.input = FALSE, deriv = TRUE, deriv2 = need_hess)
    
    # Extract the full Jacobian of f wrt (x,p)
    J <- E$jacobian[,,1, drop = TRUE]
    n <- length(states)
    m <- length(nonstates)
    
    # Partition Jacobian into df/dx and df/dp blocks
    dfdx <- J[, seq_len(n), drop = FALSE]
    dfdp <- J[, n + seq_len(m), drop = FALSE]
    rownames(dfdx) <- states
    colnames(dfdx) <- states
    rownames(dfdp) <- states
    colnames(dfdp) <- nonstates
    
    # -------------------------------------------------------------------------
    # First-order implicit function theorem:  dx/dp = - (df/dx)^{-1} * df/dp
    # -------------------------------------------------------------------------
    D_full <- solve(dfdx, -dfdp)
    
    # Apply chain rule for x = exp(z):  dx/dp = diag(x) * dz/dp
    if (positive) D_full <- diag(x_star) %*% D_full
    
    # -------------------------------------------------------------------------
    # Build outer Jacobian for all parameters (x + untouched ones)
    # -------------------------------------------------------------------------
    cols_nonfixed <- setdiff(names(pars), names(fixed))
    out_vec <- c(x_star, p[setdiff(names(p), states)])  # full parameter set
    
    # Initialize Jacobian of correct size (rows = outputs, cols = active inputs)
    J_outer <- matrix(0, nrow = length(out_vec), ncol = length(cols_nonfixed),
                      dimnames = list(names(out_vec), cols_nonfixed))
    
    # Identity mapping for pass-through parameters
    passthrough <- setdiff(names(out_vec), states)
    if (length(passthrough))
      J_outer[passthrough, intersect(passthrough, cols_nonfixed)] <-
      diag(1, nrow = length(intersect(passthrough, cols_nonfixed)))
    
    # Fill in Jacobian block for implicit states
    active_cols <- intersect(nonstates, cols_nonfixed)
    if (length(active_cols))
      J_outer[states, active_cols] <- D_full[, active_cols, drop = FALSE]
    
    # Compose with upstream derivative if this trafo is nested
    if (!is.null(dP)) {
      dP_use <- dP[cols_nonfixed, , drop = FALSE]
      myderiv <- einsum::einsum("ij,jk->ik", J_outer, dP_use)
    } else {
      myderiv <- J_outer
    }
    
    # -------------------------------------------------------------------------
    # Second-order derivatives: implicit function theorem (2nd order)
    # -------------------------------------------------------------------------
    myderiv2 <- NULL
    if (deriv2) {
      
      # Extract 2nd derivatives of f wrt (x,p)
      H <- E$hessian[,,,1, drop = TRUE]
      H_xx <- H[, seq_len(n), seq_len(n), drop = FALSE]
      H_xp <- H[, seq_len(n), n + seq_len(m), drop = FALSE]
      H_pp <- H[, n + seq_len(m), n + seq_len(m), drop = FALSE]
      
      # Restrict to active columns only (for performance)
      idx_active <- match(active_cols, nonstates)
      D <- D_full[, idx_active, drop = FALSE]
      m_act <- ncol(D)
      
      # Compute the right-hand side tensor for d2x/dp^2
      # RHS = f_xx[dx/dp,dx/dp] + f_xp*dx/dp + f_xp^T*dx/dp + f_pp
      term1 <- einsum::einsum("iab,aj,bk->ijk", H_xx, D, D)
      term2 <- einsum::einsum("iak,aj->ijk", H_xp, D)
      term3 <- einsum::einsum("iaj,ak->ijk", H_xp, D)
      term4 <- H_pp[, idx_active, idx_active, drop = FALSE]
      RHS <- term1 + term2 + term3 + term4
      
      # Solve -f_x^{-1} * RHS for each (j,k)
      sol <- -solve(dfdx, matrix(RHS, nrow = n, ncol = m_act * m_act))
      d2x <- array(sol, dim = c(n, m_act, m_act),
                   dimnames = list(states, active_cols, active_cols))
      
      # Apply second-order chain rule for x = exp(z)
      # d2x = diag(x)*(dzdp dzdp^T + d2zdp2)
      if (positive) {
        termA <- einsum::einsum("ij,kj->ijk", D, D)
        d2x <- einsum::einsum("i,ijk->ijk", x_star, termA + d2x)
      }
      
      # Embed into full Hessian (states + passthrough)
      H_outer <- array(0, dim = c(length(out_vec),
                                  length(cols_nonfixed),
                                  length(cols_nonfixed)),
                       dimnames = list(names(out_vec),
                                       cols_nonfixed,
                                       cols_nonfixed))
      H_outer[states, active_cols, active_cols] <- d2x
      
      # Compose with upstream Hessians if nested
      if (!is.null(dP) && !is.null(dP2)) {
        dP_use  <- dP [cols_nonfixed, , drop = FALSE]
        dP2_use <- dP2[cols_nonfixed, , , drop = FALSE]
        termA   <- einsum::einsum("imn,mj,nk->ijk", H_outer, dP_use, dP_use)
        termB   <- einsum::einsum("im,mjk->ijk",   J_outer, dP2_use)
        myderiv2 <- termA + termB
      } else {
        myderiv2 <- H_outer
      }
    }
    
    # -------------------------------------------------------------------------
    # Assemble the final parameter vector with derivative attributes
    # -------------------------------------------------------------------------
    res <- as.parvec(out_vec,
                     deriv  = if (deriv)  myderiv  else FALSE,
                     deriv2 = if (deriv2) myderiv2 else FALSE)
    
    # Optionally attach untouched input parameters that are not in the output
    if (attach.input && !all(names(pars) %in% names(res))) {
      res <- c(res,
               as.parvec(pars[setdiff(names(pars), names(res))],
                         deriv  = if (deriv) NULL else FALSE,
                         deriv2 = if (deriv2) NULL else FALSE))
    }
    
    res
  }
  
  # ---------------------------------------------------------------------------
  # Attach metadata for consistency with Pexpl and return as parfn
  # ---------------------------------------------------------------------------
  attr(p2p, "equations")  <- as.eqnvec(trafo)
  attr(p2p, "parameters") <- parameters
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