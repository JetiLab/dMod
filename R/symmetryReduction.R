# reduceSymmetry(): constructive removal of the non-identifiable directions a
# symmetryDetection() result reports. Scaling directions are gauged exactly through
# the integer weight lattice (transversal pins + invariant monomials); curved
# directions go through module reduction and an escalating exact invariant search
# (monomial -> polynomial <= dPoly -> rational via Darboux <= dDarboux), every
# failed stage leaving a negative certificate. All arithmetic is exact: integer
# lattice kernels through the Python module, sampling over GF(p) with CRT
# reconstruction, and an independent verify prime for the final X(I) = 0 check.

# ---- coordinates and weight extraction -------------------------------------------

# the full coordinate list of the analysis; results predating info$coordinates fall
# back to the union of the direction supports (the trafo then covers only those)
.symRedCoordinates <- function(object) {
  co <- object$info$coordinates
  if (!is.null(co) && length(co)) return(as.character(co))
  warning("reduceSymmetry(): this result carries no coordinate list (older dMod); ",
          "the emitted trafo covers only coordinates appearing in some direction -- ",
          "extend it with identity entries before use in P().", call. = FALSE)
  sort(unique(unlist(lapply(object$symmetries, .symCoords))))
}

# The integer weight rows of the scaling directions, parsed PER DIRECTION (unlike
# .symWeights, which drops everything when one weight is symbolic). A scaling with
# symbolic weights (e.g. a Hill recast "-nhill") is routed to the curved machinery.
# Returns W (rows = clean scalings, columns = union of weight names), `rows` (their
# indices into syms) and `symbolic` (indices of scaling directions W cannot hold).
.symRedWeightRows <- function(syms) {
  isScal <- vapply(syms, function(d)
    isTRUE(d$type == "scaling") && !is.null(d$weights), logical(1))
  rows <- integer(0); symbolic <- integer(0); wlist <- list()
  for (i in which(isScal)) {
    w <- suppressWarnings(as.numeric(unlist(syms[[i]]$weights)))
    if (anyNA(w) || any(w != round(w))) symbolic <- c(symbolic, i)
    else { rows <- c(rows, i); wlist[[length(wlist) + 1L]] <-
             setNames(as.integer(w), names(syms[[i]]$weights)) }
  }
  cols <- unique(unlist(lapply(wlist, names)))
  W <- matrix(0L, length(wlist), length(cols), dimnames = list(NULL, cols))
  for (r in seq_along(wlist)) W[r, names(wlist[[r]])] <- wlist[[r]]
  list(W = W, rows = rows, symbolic = symbolic)
}

# ---- small exact-arithmetic utilities --------------------------------------------

# rank over GF(p) (screening; exact verdicts go through the integer kernel)
.symRedRankModP <- function(M, p = .symPrimes[1]) {
  if (!length(M) || nrow(M) == 0L || ncol(M) == 0L) return(0L)
  .symRrefModp(M %% p, p)$rank
}

# primitive integer kernel basis of an integer matrix (columns of the returned
# matrix), exact through the Python module (GF(p) + CRT + integer validation)
.symRedIntKernel <- function(M, sd) {
  if (ncol(M) == 0L) return(matrix(0L, 0L, 0L))
  # as.list per row so a length-1 row still reaches Python as a list, not a scalar
  ker <- sd$exactIntKernel(lapply(seq_len(nrow(M)), function(r)
    as.list(as.integer(M[r, ]))), ncol(M))
  if (!length(ker)) return(matrix(0L, ncol(M), 0L))
  do.call(cbind, lapply(ker, as.integer))
}

# divide each row by its gcd (primitive integer rows)
.symRedPrimitiveRows <- function(M) {
  if (!nrow(M)) return(M)
  g <- apply(abs(M), 1L, function(r) Reduce(function(a, b) {
    while (b) { t <- b; b <- a %% b; a <- t }; a }, r[r != 0], accumulate = FALSE))
  g[!is.finite(g) | g == 0] <- 1
  round(M / g)
}

# an exponent vector as a monomial string "a*b^-1" ("1" if empty); `^` display
.symRedMonoString <- function(a, vars) {
  nz <- which(a != 0L)
  if (!length(nz)) return("1")
  paste(ifelse(a[nz] == 1L, vars[nz], paste0(vars[nz], "^", a[nz])), collapse = "*")
}

# ---- scaling stage ---------------------------------------------------------------

# Apply the user's `fixed` set to one scaling block: the surviving scaling space is
# the left kernel of W[, S] (combinations c with c^T W vanishing on every fixed
# column). Exact: kernel of t(W[, S]) through the integer-kernel route. Returns the
# residual weight rows (primitive), the removal count and the redundant fixed
# columns (raise no rank, the .symCatFixing peel).
.symRedScalingFixed <- function(W, fixed, sd) {
  S <- intersect(fixed, colnames(W))
  k <- nrow(W)
  if (!length(S)) return(list(Wres = W, removed = 0L, redundant = character(0),
                              combo = diag(1L, k)))
  C <- .symRedIntKernel(t(W[, S, drop = FALSE]), sd)      # k x k' residual combos
  Wres <- if (ncol(C)) .symRedPrimitiveRows(t(C) %*% W) else W[0L, , drop = FALSE]
  red <- character(0); acc <- W[, 0L, drop = FALSE]
  for (cn in S) {
    aug <- cbind(acc, W[, cn])
    if (.symRedRankModP(aug) == .symRedRankModP(acc)) red <- c(red, cn)
    else acc <- aug
  }
  list(Wres = Wres, removed = k - ncol(C), redundant = red, combo = C)
}

# Exact certificate that T is an admissible transversal of W (k rows): A = W[, T]
# invertible with M = A^{-1} W INTEGER (so every survivor's absorbed monomial stays
# Laurent). Integer certification without sympy: round the double solve, then
# verify A %*% M == W exactly in integer arithmetic (entries far below 2^53).
.symRedCertifyTransversal <- function(W, T) {
  A <- W[, T, drop = FALSE]
  if (nrow(A) != length(T)) return(NULL)
  M <- tryCatch(solve(A, W), error = function(e) NULL)
  if (is.null(M)) return(NULL)
  Mr <- round(M)
  if (max(abs(M - Mr)) > 1e-7) return(NULL)
  if (!all(A %*% Mr == W)) return(NULL)
  storage.mode(Mr) <- "integer"
  Mr
}

# The transversal machinery of one scaling block: a certified representative T
# (deterministic: coordinates entering fewest rows first), the admissible family
# (full enumeration when small, else the matroid description "one per RREF row,
# choices independent"), and the survivor meanings z_j * prod_t z_t^{-M_tj}.
.symRedTransversal <- function(Wres, fixed) {
  k <- nrow(Wres)
  supp <- colnames(Wres)[colSums(Wres != 0L) > 0L]
  supp <- setdiff(supp, fixed)                     # fixed columns are already pinned
  out <- list(T = NULL, M = NULL, admissible = NULL, matroid = NULL,
              survivorMeaning = NULL)
  if (!k) return(out)
  Ws <- Wres[, supp, drop = FALSE]
  occ <- colSums(Ws != 0L)
  ord <- supp[order(occ, match(supp, colnames(Wres)))]

  certify <- function(T) .symRedCertifyTransversal(Wres, T)

  enumerate <- length(supp) <= 30L && choose(length(supp), k) <= 200
  if (enumerate) {
    combos <- utils::combn(supp, k, simplify = FALSE)
    adm <- Filter(function(T) .symRedRankModP(Wres[, T, drop = FALSE]) == k, combos)
    out$admissible <- adm
    for (T in adm) {                               # first certified in combn order
      M <- certify(T)
      if (!is.null(M)) { out$T <- T; out$M <- M; break }
    }
  } else {
    # greedy representative; family described through the RREF row supports
    T <- character(0)
    for (cn in ord) {
      if (length(T) == k) break
      if (.symRedRankModP(Wres[, c(T, cn), drop = FALSE]) > length(T)) T <- c(T, cn)
    }
    if (length(T) == k) { M <- certify(T); if (!is.null(M)) { out$T <- T; out$M <- M } }
    ref <- .symRrefModp(Wres[, supp, drop = FALSE] %% .symPrimes[1], .symPrimes[1])
    out$matroid <- lapply(seq_len(ref$rank), function(r)
      supp[which(ref$R[r, ] != 0)])
  }
  if (!is.null(out$T)) {
    surv <- setdiff(supp, out$T)
    out$survivorMeaning <- setNames(vapply(surv, function(j) {
      ex <- setNames(integer(length(out$T) + 1L), c(j, out$T))
      ex[j] <- 1L; ex[out$T] <- -out$M[, j]
      .symRedMonoString(ex, names(ex))
    }, character(1)), surv)
  }
  out
}

# One scaling block end to end: fixed -> residual space -> invariant lattice ->
# transversal -> pins. Status: "fixed" (nothing left), "reduced" (certified
# transversal) or "invariantOnly" (no Laurent transversal; lattice still reported).
.symRedScalingBlock <- function(W, labels, fixed, sd) {
  fx <- .symRedScalingFixed(W, fixed, sd)
  Wres <- fx$Wres
  blk <- list(labels = labels, type = "scaling", support = colnames(W)[colSums(W != 0L) > 0L],
              stage = "transversal", moduleCombos = NULL,
              removedByFixed = fx$removed, redundantFixed = fx$redundant,
              certificates = character(0), reason = NULL)
  if (!nrow(Wres)) {
    blk$status <- "fixed"; blk$stage <- "fixed"
    blk$invariants <- character(0); blk$transversal <- NULL
    blk$gaugeNote <- "removed entirely by the user's fixed coordinates"
    return(blk)
  }
  K <- .symRedIntKernel(Wres[, blk$support, drop = FALSE], sd)
  rownames(K) <- blk$support
  blk$invExps <- K
  blk$invariants <- if (ncol(K)) vapply(seq_len(ncol(K)), function(j)
    .symRedMonoString(K[, j], blk$support), character(1)) else character(0)
  tv <- .symRedTransversal(Wres, fixed)
  blk$transversal <- tv$T
  blk$admissible <- tv$admissible
  blk$matroid <- tv$matroid
  blk$survivorMeaning <- tv$survivorMeaning
  blk$Wres <- Wres
  if (!is.null(tv$T)) {
    blk$status <- "reduced"
    blk$gaugeNote <- paste0("each transversal coordinate may be pinned to any ",
                            "nonzero constant (representative uses 1); the family ",
                            "ranges over all admissible transversals")
    blk$pins <- setNames(rep("1", length(tv$T)), tv$T)
    blk$certificates <- c(blk$certificates,
      "transversal certified exactly: W[,T] invertible with integer W[,T]^-1 W")
  } else {
    blk$status <- "invariantOnly"
    blk$reason <- "no transversal with an integer (Laurent) absorption matrix"
    blk$pins <- NULL
  }
  blk
}

# connected components of the weight rows through shared nonzero columns
.symRedComponents <- function(W) {
  k <- nrow(W)
  if (!k) return(list())
  parent <- seq_len(k)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
    if (any(W[i, ] != 0L & W[j, ] != 0L)) parent[find(j)] <- find(i)
  split(seq_len(k), vapply(seq_len(k), find, integer(1)))
}

# ---- curved stage: generator prep, module reduction, block decomposition --------

# Locals table for sympify: every identifier mapped to a plain Symbol, so model
# names never collide with sympy builtins (Ci is the cosine integral, E is Euler's
# number, ...). The lookbehind keeps the exponent of "1e-5" from matching.
.symRedLocals <- function(strings, spy) {
  ids <- unique(unlist(regmatches(strings,
    gregexpr("(?<![0-9.A-Za-z_])[A-Za-z_][A-Za-z0-9_.]*", strings, perl = TRUE))))
  ids <- setdiff(ids, "")
  reticulate::dict(setNames(lapply(ids, function(x) spy$Symbol(x)), ids))
}

.symRedSympify <- function(x, spy, locals)
  spy$sympify(as.character(x), locals = locals)

# the free symbols of a sympy expression, as a character vector
.symRedFreeSyms <- function(e) {
  syms <- reticulate::iterate(e$free_symbols, function(s) as.character(s))
  sort(as.character(unlist(syms)))
}

# One curved generator canonicalised for the invariant search: components sympified,
# denominators cleared across the WHOLE generator (Q*X has the same invariants as X
# off Q = 0), re-serialised as polynomial strings. `vars` carries the support plus
# every coefficient symbol (a known input like `u` may enter an invariant).
.symRedGenPrep <- function(d, spy) {
  out <- list(ok = FALSE, comps = NULL, support = NULL, vars = NULL, degree = NA_integer_)
  locals <- .symRedLocals(unlist(lapply(d$generator, as.character)), spy)
  exprs <- tryCatch(lapply(d$generator, function(x) .symRedSympify(x, spy, locals)),
                    error = function(e) NULL)
  if (is.null(exprs)) { out$reason <- "component not parseable by sympy"; return(out) }
  dens <- lapply(exprs, function(e) spy$fraction(spy$cancel(spy$together(e)))[[2]])
  L <- Reduce(function(a, b) spy$lcm(a, b), dens, spy$Integer(1L))
  exprs <- lapply(exprs, function(e) spy$expand(spy$cancel(e * L)))
  isPoly <- vapply(exprs, function(e)
    tryCatch(isTRUE(e$is_polynomial()), error = function(err) FALSE), logical(1))
  if (!all(isPoly)) { out$reason <- "component not a rational function"; return(out) }
  comps <- vapply(exprs, function(e) as.character(e), character(1))
  names(comps) <- names(d$generator)
  nz <- comps != "0"
  out$ok <- TRUE
  out$comps <- comps[nz]
  out$support <- names(comps)[nz]
  out$vars <- sort(unique(c(out$support,
    unlist(lapply(exprs[nz], .symRedFreeSyms)))))
  degs <- vapply(exprs[nz], function(e) {
    # total_degree may hand back a sympy Integer; the string route converts both
    dg <- tryCatch(suppressWarnings(as.integer(as.character(spy$total_degree(e)))),
                   error = function(err) NA_integer_)
    if (length(dg) != 1L) NA_integer_ else dg
  }, integer(1))
  out$degree <- if (all(is.na(degs))) NA_integer_ else max(degs, na.rm = TRUE)
  out
}

# connected components of index sets through shared elements (union-find); returns
# a list of integer index groups
.symRedGroupBy <- function(sets) {
  k <- length(sets)
  if (!k) return(list())
  parent <- seq_len(k)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k)
    if (length(intersect(sets[[i]], sets[[j]]))) parent[find(j)] <- find(i)
  unname(split(seq_len(k), vapply(seq_len(k), find, integer(1))))
}

# Pivot-schedule discovery for the module reduction: RREF-style elimination on the
# value matrix over GF(p), at each step choosing the (row, column) whose elimination
# leaves the fewest nonzeros. Returns the schedule and its final fill.
.symRedReduceSchedule <- function(V, p) {
  k <- nrow(V); n <- ncol(V)
  used <- logical(k); sched <- list()
  for (step in seq_len(k)) {
    best <- NULL
    for (r in which(!used)) for (cc in which(V[r, ] != 0)) {
      Vt <- V
      piv <- Vt[r, cc]
      for (j in seq_len(k)) if (j != r && Vt[j, cc] != 0) {
        m <- .symRedMulmodP(Vt[j, cc], .symInvmod(piv, p), p)
        Vt[j, ] <- (Vt[j, ] - .symRedMulmodP(V[r, ], m, p)) %% p
      }
      fill <- sum(Vt != 0)
      if (is.null(best) || fill < best$fill) best <- list(r = r, c = cc, fill = fill, V = Vt)
    }
    if (is.null(best)) break
    used[best$r] <- TRUE
    sched[[length(sched) + 1L]] <- c(best$r, best$c)
    V <- best$V
  }
  list(sched = sched, fill = sum(V != 0))
}

# vectorised a*b mod p in doubles (wrapper around the scalar-b .symMulmod)
.symRedMulmodP <- function(a, b, p) .symMulmod(a %% p, b %% p, p)

# Module reduction of a curved block: generators may be combined with FUNCTION
# coefficients (pointwise span decides identifiability -- a module over the function
# field, not a vector space), so RREF with symbolic multipliers is legitimate and
# can collapse the support (M009: X13 - r2*X12 went from 7 to 3 coordinates). The
# pivot schedule is discovered on a value matrix over GF(p) (cheap, fill-minimising,
# cross-checked on a second prime), then replayed exactly in sympy. The heuristic
# affects only how small the support gets, never correctness.
.symRedModuleReduce <- function(preps, labels, sd, spy) {
  k <- length(preps)
  cols <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(unlist(lapply(preps, `[[`, "vars"))))
  if (k < 2L)
    return(list(preps = preps, labelsOf = as.list(seq_len(k)), combos = character(0)))

  # value matrices at one random point per prime; schedule = the lower-fill one
  scheds <- lapply(.symPrimes[1:2], function(p) {
    pool <- .symPool()
    pt <- setNames(pool(seq_along(vars) + 7L), vars)
    V <- matrix(0, k, length(cols))
    for (r in seq_len(k)) for (ci in seq_along(cols)) {
      cc <- cols[ci]
      if (!is.null(preps[[r]]$comps[cc]) && !is.na(preps[[r]]$comps[cc])) {
        v <- .symEvalModq(preps[[r]]$comps[[cc]], pt, p, sd)
        V[r, ci] <- if (is.na(v)) 0 else v
      }
    }
    .symRedReduceSchedule(V, p)
  })
  sched <- scheds[[which.min(vapply(scheds, `[[`, numeric(1), "fill"))]]$sched

  # exact symbolic replay of the schedule; C tracks the combination coefficients
  locals <- .symRedLocals(unlist(lapply(preps, `[[`, "comps")), spy)
  E <- lapply(preps, function(pr) {
    row <- setNames(vector("list", length(cols)), cols)
    for (cc in cols) row[[cc]] <- .symRedSympify(
      if (cc %in% names(pr$comps)) pr$comps[[cc]] else "0", spy, locals)
    row
  })
  C <- lapply(seq_len(k), function(i)
    lapply(seq_len(k), function(j) spy$Integer(as.integer(i == j))))
  zero <- function(e) identical(as.character(e), "0")
  for (st in sched) {
    r <- st[1]; cc <- cols[st[2]]
    piv <- E[[r]][[cc]]
    if (zero(piv)) next
    for (j in seq_len(k)) {
      if (j == r || zero(E[[j]][[cc]])) next
      m <- spy$cancel(E[[j]][[cc]] / piv)
      for (c2 in cols) E[[j]][[c2]] <- spy$cancel(E[[j]][[c2]] - m * E[[r]][[c2]])
      for (l in seq_len(k)) C[[j]][[l]] <- spy$cancel(C[[j]][[l]] - m * C[[r]][[l]])
    }
  }

  # Candidates = the original generators (identity combination) plus every reduced
  # row; a greedy pass picks a rank-complete basis with the smallest supports, so a
  # row is only ever replaced when the replacement's support is genuinely smaller.
  cand <- lapply(seq_len(k), function(j)
    list(pr = preps[[j]], srcIdx = j, combo = NULL))
  for (j in seq_len(k)) {
    comps <- vapply(cols, function(cc) as.character(E[[j]][[cc]]), character(1))
    names(comps) <- cols
    pr <- .symRedGenPrep(list(generator = as.list(comps[comps != "0"])), spy)
    if (!pr$ok || !length(pr$support)) next            # row eliminated entirely
    srcIdx <- which(vapply(seq_len(k), function(l) !zero(C[[j]][[l]]), logical(1)))
    if (length(srcIdx) == 1L && srcIdx == j &&
        identical(as.character(C[[j]][[j]]), "1")) next      # unchanged: already in
    combo <- paste0(
      paste0("(", vapply(srcIdx, function(l) as.character(C[[j]][[l]]), character(1)),
             ")*", labels[srcIdx], collapse = " + "),
      "  [support ", length(pr$support), "]")
    cand[[length(cand) + 1L]] <- list(pr = pr, srcIdx = srcIdx, combo = combo)
  }

  p <- .symPrimes[1]
  pool <- .symPool()
  pt <- setNames(pool(seq_along(vars) + 23L), vars)
  vrow <- function(pr) vapply(cols, function(cc) {
    if (!(cc %in% names(pr$comps))) return(0)
    v <- .symEvalModq(pr$comps[[cc]], pt, p, sd)
    if (is.na(v)) 0 else as.numeric(v)
  }, numeric(1))
  ord <- order(vapply(cand, function(cd) length(cd$pr$support), integer(1)))
  acc <- matrix(0, 0L, length(cols))
  outPreps <- list(); labelsOf <- list(); combos <- character(0)
  for (j in ord) {
    if (nrow(acc) == k) break
    aug <- rbind(acc, vrow(cand[[j]]$pr))
    if (.symRedRankModP(aug, p) > nrow(acc)) {
      acc <- aug
      outPreps[[length(outPreps) + 1L]] <- cand[[j]]$pr
      labelsOf[[length(labelsOf) + 1L]] <- cand[[j]]$srcIdx
      if (!is.null(cand[[j]]$combo)) combos <- c(combos, cand[[j]]$combo)
    }
  }
  list(preps = outPreps, labelsOf = labelsOf, combos = combos)
}

# ---- invariant stages ------------------------------------------------------------

# The log-derivative rows of a block at one random point mod p: one row per
# generator, entries eta_i = xi_i(z)/z_i over the moved coordinates. NULL when an
# evaluation fails (redraw the point).
.symRedEtaRows <- function(preps, vars, pt, p, sd) {
  rows <- matrix(0, 0L, length(vars))
  for (pr in preps) {
    row <- numeric(length(vars))
    for (i in seq_along(vars)) {
      v <- vars[i]
      if (!(v %in% names(pr$comps))) next
      xv <- .symEvalModq(pr$comps[[v]], pt, p, sd)
      if (is.na(xv)) return(NULL)
      row[i] <- .symMulmod(xv, .symInvmod(pt[[v]] %% p, p), p)
    }
    rows <- rbind(rows, row)
  }
  rows
}

# Common nullspace of row-sampling matrices over the 4 primes, CRT-reconstructed to
# primitive integer columns. `rowsAt(pt, p)` delivers the rows of one point (NULL =
# redraw). Sampling saturates when the GF(p) rank is flat twice. Returns the integer
# basis (columns), or NULL when pivots disagree across every prime pair or an entry
# does not reconstruct, and always the row count for the certificate text.
.symRedModularNullspace <- function(rowsAt, n, evalVars, primes = .symPrimes) {
  basisPer <- list(); pivotsPer <- list(); nRows <- 0L
  nv <- length(evalVars)
  for (pi in seq_along(primes)) {
    p <- primes[pi]
    pool <- .symPool()
    rows <- matrix(0, 0L, n); prev <- -1L; flat <- 0L
    off <- 50L * pi; draws <- 0L
    repeat {
      pt <- setNames(pool(seq_len(nv) + off), evalVars); off <- off + nv
      draws <- draws + 1L
      new <- rowsAt(pt, p)
      if (is.null(new)) { if (draws > 25L) break else next }
      rows <- rbind(rows, new)
      r <- .symRedRankModP(rows, p)
      if (r == n) break
      flat <- if (r == prev) flat + 1L else 0L
      prev <- r
      if (flat >= 2L || nrow(rows) > 4L * n + 20L) break
    }
    nRows <- max(nRows, nrow(rows))
    if (!nrow(rows)) return(list(basis = NULL, rows = 0L))   # sampling never succeeded
    ref <- .symRrefModp(rows %% p, p)
    free <- setdiff(0:(n - 1L), ref$piv)
    if (!length(free))
      return(list(basis = matrix(0L, n, 0L), rows = nRows))
    pivotsPer[[pi]] <- ref$piv
    basisPer[[pi]] <- .symNullspaceBasis(list(dim = n, pivots = ref$piv, R = ref$R),
                                         free, p)
  }
  agree <- vapply(pivotsPer, function(pv) identical(pv, pivotsPer[[1]]), logical(1))
  use <- which(agree)
  if (length(use) < 2L) return(list(basis = NULL, rows = nRows))
  nf <- ncol(basisPer[[use[1]]])
  residues <- vapply(use, function(jj) as.integer(basisPer[[jj]]),
                     integer(n * nf))
  rec <- symRatRecon(matrix(residues, n * nf, length(use)), as.integer(primes[use]))
  if (any(rec$den == "0")) return(list(basis = NULL, rows = nRows))
  num <- suppressWarnings(as.numeric(rec$num)); den <- suppressWarnings(as.numeric(rec$den))
  if (anyNA(num) || anyNA(den) || any(abs(num) > 2^40) || any(den > 2^40))
    return(list(basis = NULL, rows = nRows))
  B <- matrix(num / den, n, nf)
  for (j in seq_len(nf)) {                    # clear each column to primitive integers
    dj <- den[(j - 1L) * n + seq_len(n)]
    L <- Reduce(function(a, b) a * b / .symRedGcd(a, b), unique(dj), 1)
    B[, j] <- B[, j] * L
  }
  B <- round(B)
  storage.mode(B) <- "integer"
  list(basis = t(.symRedPrimitiveRows(t(B))), rows = nRows)
}

.symRedGcd <- function(a, b) { a <- abs(a); b <- abs(b)
  while (b) { t <- b; b <- a %% b; a <- t }; a }

# Stage 1: Laurent-monomial invariants. m = prod z_i^{a_i} is invariant under X iff
# sum_i a_i * xi_i/z_i vanishes identically -- linear in the exponents a, one row per
# (generator, point). Exponents run over the MOVED coordinates only: a coefficient
# symbol (xi identically 0) would only multiply in a trivial constant factor.
.symRedMonomialInvariants <- function(preps, sd) {
  vars <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  evalVars <- sort(unique(c(vars, unlist(lapply(preps, `[[`, "vars")))))
  n <- length(vars)
  ns <- .symRedModularNullspace(function(pt, p)
    .symRedEtaRows(preps, vars, pt, p, sd), n, evalVars)
  cert <- sprintf(paste0("monomial stage: exact modular nullspace over %d primes, ",
                         "%d sampling rows"), length(.symPrimes), ns$rows)
  if (is.null(ns$basis))
    return(list(ok = FALSE, invariants = character(0), exps = NULL,
                cert = paste0(cert, " -- reconstruction failed (inconclusive)")))
  if (!ncol(ns$basis))
    return(list(ok = FALSE, invariants = character(0), exps = NULL,
                cert = paste0("no Laurent-monomial invariant (", cert, ")")))
  inv <- vapply(seq_len(ncol(ns$basis)), function(j)
    .symRedMonoString(ns$basis[, j], vars), character(1))
  list(ok = TRUE, invariants = inv, exps = ns$basis, vars = vars, cert = cert)
}

# an integer coefficient vector over a monomial table as a polynomial string
.symRedPolyString <- function(cf, expts, vars) {
  nz <- which(cf != 0L)
  terms <- vapply(nz, function(k) {
    mono <- .symRedMonoString(expts[k, ], vars)
    a <- cf[k]
    if (mono == "1") as.character(a)
    else if (a == 1L) mono
    else if (a == -1L) paste0("-", mono)
    else paste0(a, "*", mono)
  }, character(1))
  out <- terms[1]
  for (t in terms[-1]) out <- paste0(out, if (startsWith(t, "-")) " - " else " + ",
                                     sub("^-", "", t))
  out
}

# Stage 2: polynomial invariants of total degree <= d. For I = sum_k c_k m_k the
# condition X(I) = 0 is LINEAR in c: X(m_k) = m_k * sum_i a^(k)_i xi_i/z_i, so a row
# per (generator, point) has entries m_k(z) * <a^(k), eta(z)> mod p. The ansatz runs
# over the block support plus the coefficient symbols, but monomials with no moved
# coordinate are dropped up front: X kills them identically, so they contribute only
# the trivial invariants (the observed 2/5/9 pattern) and excluding them is lossless
# up to additive constants.
.symRedPolyInvariants <- function(preps, dPoly, sd) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  expts <- .symMonoTable(length(vars), as.integer(dPoly))
  keep <- rowSums(expts[, match(moved, vars), drop = FALSE]) > 0L
  expts <- expts[keep, , drop = FALSE]
  N <- nrow(expts)
  if (!N) return(list(ok = FALSE, invariants = character(0),
                      cert = "polynomial stage skipped: empty ansatz"))
  rowsAt <- function(pt, p) {
    eta <- .symRedEtaRows(preps, vars, pt, p, sd)          # nGen x nvars
    if (is.null(eta)) return(NULL)
    mv <- symMonoResidues(expts, as.integer(pt %% p), p)   # monomial values
    t(vapply(seq_len(nrow(eta)), function(g)
      vapply(seq_len(N), function(k)
        .symRedMulmodP(mv[k], sum(.symRedMulmodP(expts[k, ], eta[g, ], p)) %% p, p),
        numeric(1)),
      numeric(N)))
  }
  ns <- .symRedModularNullspace(rowsAt, N, vars)
  cert <- sprintf(paste0("polynomial stage (degree <= %d, %d monomials): exact ",
                         "modular nullspace over %d primes, %d sampling rows"),
                  dPoly, N, length(.symPrimes), ns$rows)
  if (is.null(ns$basis))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0(cert, " -- reconstruction failed (inconclusive)")))
  if (!ncol(ns$basis))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no polynomial invariant of total degree <= ", dPoly,
                              " beyond functions of the unmoved coordinates (",
                              cert, ")")))
  sel <- .symRedIndependentPoly(ns$basis, expts, vars, moved, preps)
  inv <- vapply(sel, function(j)
    .symRedPolyString(ns$basis[, j], expts, vars), character(1))
  list(ok = TRUE, invariants = inv,
       coefs = ns$basis[, sel, drop = FALSE], expts = expts, vars = vars,
       cert = paste0(cert, sprintf("; %d functionally independent of %d",
                                   length(sel), ncol(ns$basis))))
}

# The linear nullspace over-generates: powers and products of low-degree invariants
# reappear as independent basis vectors. Select a functionally independent subset by
# a greedy Jacobian rank test at a generic point, lowest degree / fewest terms
# first; the target count is #moved - rank(xi), the dimension of the invariant
# foliation. Selection only -- correctness rests on the verify layer.
.symRedIndependentPoly <- function(B, expts, vars, moved, preps) {
  pt <- setNames(as.numeric(.symPool()(seq_along(vars) + 5L)), vars)
  Xi <- t(vapply(preps, function(pr) vapply(moved, function(v)
    if (v %in% names(pr$comps))
      eval(parse(text = pr$comps[[v]]), as.list(pt)) else 0,
    numeric(1)), numeric(length(moved))))
  target <- length(moved) - qr(Xi, tol = 1e-9)$rank
  mv <- exp(expts %*% log(pt))                       # monomial values at pt
  grad <- function(cf) vapply(seq_along(vars), function(i)
    sum(cf * mv * expts[, i] / pt[i]), numeric(1))
  deg <- apply(expts, 1L, sum)
  key <- vapply(seq_len(ncol(B)), function(j)
    max(deg[B[, j] != 0L]) * 1e4 + sum(B[, j] != 0L), numeric(1))
  J <- matrix(0, 0L, length(vars)); sel <- integer(0)
  for (j in order(key)) {
    if (length(sel) >= target && target > 0L) break
    Jc <- rbind(J, grad(B[, j]))
    if (qr(Jc, tol = 1e-9)$rank > nrow(J)) { J <- Jc; sel <- c(sel, j) }
  }
  if (!length(sel)) sel <- seq_len(min(1L, ncol(B)))
  sel
}

# symbolic determinant of the extactic matrix stays tractable only for a small
# monomial basis; above the cap the stage falls back to coordinate factors
.symRedExtacticCap <- 12L

# apply the generator X = sum_i xi_i d/dz_i of one prep to a sympy expression
.symRedApplyX <- function(pr, xiExprs, f, spy, locals) {
  out <- spy$Integer(0L)
  for (v in names(xiExprs))
    out <- out + xiExprs[[v]] * spy$diff(f, .symRedSympify(v, spy, locals))
  spy$expand(out)
}

# Exact rational value of a sympy expression at an integer point, as "num/den".
# Substitution pairs use explicit Symbol objects: string keys would re-sympify and
# collide with sympy builtins (Ci, E, ...), silently skipping the substitution.
.symRedSubsRational <- function(e, pt, spy) {
  pairs <- lapply(names(pt), function(n)
    reticulate::tuple(spy$Symbol(n), spy$Integer(as.integer(pt[[n]]))))
  as.character(spy$cancel(e$subs(pairs)))
}

# Stage 3: rational invariants via Darboux polynomials (extactic route). Every
# Darboux polynomial P (X(P) = lambda*P, lambda the polynomial cofactor) of degree
# <= d divides the extactic determinant E_d of the monomial basis B_d (Pereira), so
# spy$factor(E_d) delivers ALL candidates of that degree; candidates must divide for
# every generator of the block. Rational invariants are prod P_j^{n_j} for integer
# vectors n with sum n_j lambda_j = 0 identically -- the integer kernel of the
# cofactor matrix, sampled exactly at integer points and lcm-cleared. This subsumes
# equal-cofactor pairs and the coordinate factors z_i | xi_i. Above the basis cap
# only the coordinate factors are tried (stated in the certificate).
.symRedDarboux <- function(preps, dDarboux, sd, spy, verbose = FALSE) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xiOf <- lapply(preps, function(pr)
    lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))

  # Candidate Darboux polynomials from three sources: the coordinate factors, the
  # irreducible factors of the xi components themselves (an invariant hypersurface
  # often sits inside the vanishing of xi -- and when a first integral makes the
  # extactic degenerate, these are the candidates that remain reachable), and the
  # extactic factors when the basis stays below the cap. The exact division test
  # filters, so extra candidates never cost correctness.
  cands <- lapply(moved, function(v) .symRedSympify(v, spy, locals))
  for (xi in unlist(xiOf, recursive = FALSE)) {
    fl <- tryCatch(spy$factor_list(xi), error = function(e) NULL)
    if (!is.null(fl)) cands <- c(cands, lapply(fl[[2]], function(pair) pair[[1]]))
  }
  Nbasis <- choose(length(vars) + dDarboux, dDarboux)
  # the symbolic determinant's cost is driven by the entry degrees, which grow to
  # ~ dDarboux + (N-1)*(D-1) in the last extactic row -- cap that, not just N
  Dmax <- max(1L, vapply(preps, function(pr)
    if (is.na(pr$degree)) 3L else pr$degree, integer(1)))
  entryDeg <- dDarboux + (Nbasis - 1L) * (Dmax - 1L)
  extactic <- dDarboux >= 1L && Nbasis <= .symRedExtacticCap &&
              entryDeg <= .symRedExtacticCap
  degenerate <- FALSE
  if (extactic) {
    expts <- .symMonoTable(length(vars), as.integer(dDarboux))
    basis <- lapply(seq_len(nrow(expts)), function(k)
      .symRedSympify(.symRedMonoString(expts[k, ], vars), spy, locals))
    N <- length(basis)
    rows <- list(basis)
    for (j in seq_len(N - 1L))
      rows[[j + 1L]] <- lapply(rows[[j]], function(f)
        .symRedApplyX(preps[[1]], xiOf[[1]], f, spy, locals))
    Emat <- spy$Matrix(lapply(rows, function(r) r))
    D <- Emat$det(method = "berkowitz")               # the one expensive step
    nonzero <- any(vapply(1:2, function(t) {          # E identically 0 = degenerate
      pt <- setNames(as.list(.symPool()(seq_along(vars) + 11L * t)), vars)
      v <- tryCatch(.symRedSubsRational(D, pt, spy), error = function(e) NA_character_)
      !is.na(v) && v != "0"
    }, logical(1)))
    if (!nonzero) degenerate <- TRUE
    else {
      fl <- spy$factor_list(D)                      # convert=TRUE: already an R list
      cands <- c(cands, lapply(fl[[2]], function(pair) pair[[1]]))
    }
  }

  # keep candidates that are Darboux for EVERY generator; collect their cofactors
  Ps <- list(); lams <- list()
  for (P in cands) {
    if (isTRUE(P$is_number)) next
    lamRow <- vector("list", length(preps)); ok <- TRUE
    for (g in seq_along(preps)) {
      XP <- .symRedApplyX(preps[[g]], xiOf[[g]], P, spy, locals)
      q <- spy$cancel(XP / P)
      den <- spy$fraction(q)[[2]]
      if (!identical(as.character(den), "1")) { ok <- FALSE; break }
      lamRow[[g]] <- q
    }
    if (!ok) next
    dup <- any(vapply(Ps, function(Q)
      isTRUE(spy$cancel(P / Q)$is_number), logical(1)))
    if (!dup) { Ps[[length(Ps) + 1L]] <- P; lams[[length(lams) + 1L]] <- lamRow }
  }

  coverage <- if (degenerate)
      "extactic degenerate (E identically 0): coordinate and xi factors only"
    else if (extactic) sprintf("extactic complete for factor degree <= %d", dDarboux)
    else sprintf(paste0("extactic skipped (basis %d, projected entry degree %d, ",
                        "cap %d): coordinate and xi factors only"),
                 Nbasis, entryDeg, .symRedExtacticCap)
  certBase <- sprintf("Darboux stage (%d candidate factors; %s)", length(Ps), coverage)
  if (length(Ps) < 2L)
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no rational invariant from Darboux factors -- ",
                              certBase)))

  # integer kernel of the cofactor matrix, sampled exactly at integer points
  nMon <- choose(length(vars) + max(1L, dDarboux), length(vars)) + 3L
  rows <- list()
  for (g in seq_along(preps)) for (t in seq_len(nMon)) {
    pt <- setNames(as.list(.symPool()(seq_along(vars) + 101L * t + 7L * g)), vars)
    vals <- vapply(seq_along(Ps), function(j)
      .symRedSubsRational(lams[[j]][[g]], pt, spy), character(1))
    num <- suppressWarnings(as.numeric(sub("/.*$", "", vals)))
    den <- suppressWarnings(as.numeric(ifelse(grepl("/", vals),
                                              sub("^.*/", "", vals), "1")))
    if (anyNA(num) || anyNA(den)) next
    L <- Reduce(function(a, b) a * b / .symRedGcd(a, b), unique(den), 1)
    row <- round(num * (L / den))
    if (max(abs(row)) < 2^31) rows[[length(rows) + 1L]] <- as.integer(row)
  }
  if (!length(rows))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  K <- .symRedIntKernel(do.call(rbind, rows), sd)
  if (!ncol(K))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no rational invariant with Darboux factors of degree <= ",
                              dDarboux, " -- ", certBase)))
  Pstr <- vapply(Ps, function(P) gsub("\\*\\*", "^", as.character(P)), character(1))
  inv <- vapply(seq_len(ncol(K)), function(j) {
    nz <- which(K[, j] != 0L)
    paste(vapply(nz, function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      if (K[l, j] == 1L) base else paste0(base, "^", K[l, j])
    }, character(1)), collapse = "*")
  }, character(1))
  list(ok = TRUE, invariants = inv, cert = certBase)
}

# The invariant count a block needs: #moved coordinates minus the generic rank of
# the xi-matrix -- the dimension of the invariant foliation. Numeric, at a generic
# point; selection-only (correctness rests on the verify layer).
.symRedBlockCorank <- function(preps) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  pt <- setNames(as.numeric(.symPool()(seq_along(vars) + 5L)), vars)
  Xi <- t(vapply(preps, function(pr) vapply(moved, function(v)
    if (v %in% names(pr$comps))
      eval(parse(text = pr$comps[[v]]), as.list(pt)) else 0,
    numeric(1)), numeric(length(moved))))
  length(moved) - qr(Xi, tol = 1e-9)$rank
}

# a functionally independent subset of invariant strings (any stage's output, all
# R-parseable), greedy by string length under a numeric Jacobian rank test
.symRedIndependentSet <- function(invs, vars, target) {
  if (!length(invs) || target <= 0L) return(character(0))
  pt <- setNames(as.numeric(.symPool()(seq_along(vars) + 13L)), vars)
  grad <- function(iv) {
    ex <- parse(text = iv)[[1]]
    vapply(vars, function(v)
      tryCatch(eval(stats::D(ex, v), as.list(pt)), error = function(e) NA_real_),
      numeric(1))
  }
  J <- matrix(0, 0L, length(vars)); sel <- character(0)
  for (iv in unique(invs)[order(nchar(unique(invs)))]) {
    if (length(sel) >= target) break
    gr <- grad(iv)
    if (anyNA(gr)) next
    Jc <- rbind(J, gr)
    if (qr(Jc, tol = 1e-9)$rank > nrow(J)) { J <- Jc; sel <- c(sel, iv) }
  }
  sel
}

# run the invariant-search stages on one (already module-reduced) sub-block of
# generators -- monomial, then polynomial <= dPoly, then Darboux <= dDarboux.
.symRedSolveBlock <- function(preps, dPoly, dDarboux, sd, spy, verbose = FALSE) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  target <- .symRedBlockCorank(preps)
  certs <- character(0); found <- character(0); stage <- "monomial"
  pick <- function(cand) .symRedIndependentSet(cand, vars, target)

  # the stages ACCUMULATE: a block needs `target` functionally independent
  # invariants before a gauge pin is exact, so one early find must not stop the
  # escalation (the joint A<->B block: the monomial stage sees only A*alpha, the
  # other two live in the polynomial stage)
  mono <- .symRedMonomialInvariants(preps, sd)
  certs <- c(certs, mono$cert)
  if (mono$ok) found <- pick(mono$invariants)
  if (length(found) < target) {
    if (isTRUE(verbose)) message("  monomial: ", length(found), "/", target,
                                 "; escalating to degree <= ", dPoly)
    stage <- "polynomial"
    poly <- .symRedPolyInvariants(preps, dPoly, sd)
    certs <- c(certs, poly$cert)
    if (poly$ok) found <- pick(c(found, poly$invariants))
  }
  if (length(found) < target) {
    if (isTRUE(verbose)) message("  polynomial: ", length(found), "/", target,
                                 "; escalating to Darboux")
    stage <- "darboux"
    dar <- .symRedDarboux(preps, dDarboux, sd, spy, verbose)
    certs <- c(certs, dar$cert)
    if (dar$ok) found <- pick(c(found, dar$invariants))
  }
  list(stage = stage, invariants = found, target = target, certificates = certs,
       reason = if (length(found) >= target) NULL
         else sprintf(paste0("found %d of %d independent invariants (poly degree ",
                             "<= %d, Darboux degree <= %d); gauge pinning needs ",
                             "the full set"),
                      length(found), target, dPoly, dDarboux))
}

# The curved pipeline: prep every direction with a closed-form generator, group by
# shared support, module-reduce each group, re-split (a reduced row may decouple),
# then run the invariant stages per sub-block. Support-only directions and failed
# preps become unresolved singleton blocks.
.symRedCurved <- function(syms, idx, labels, fixed, dPoly, dDarboux, sd, spy,
                          verbose = FALSE) {
  blocks <- list()
  emit <- function(b) blocks[[length(blocks) + 1L]] <<- b
  unresolved <- function(labs, support, reason)
    list(labels = labs, type = "curved", support = support, stage = "none",
         status = "unresolved", invariants = character(0),
         certificates = character(0), reason = reason,
         moduleCombos = NULL, transversal = NULL, survivorMeaning = NULL,
         pins = NULL)

  hasGen <- vapply(idx, function(i) !is.null(syms[[i]]$generator), logical(1))
  for (i in idx[!hasGen])
    emit(unresolved(labels[i], .symCoords(syms[[i]]),
      "no closed-form generator; re-run symmetryDetection(reconstruct = TRUE)"))

  gi <- idx[hasGen]
  preps <- lapply(gi, function(i) .symRedGenPrep(syms[[i]], spy))
  bad <- !vapply(preps, `[[`, logical(1), "ok")
  for (w in which(bad))
    emit(unresolved(labels[gi[w]], .symCoords(syms[[gi[w]]]), preps[[w]]$reason))
  gi <- gi[!bad]; preps <- preps[!bad]
  if (!length(gi)) return(blocks)

  for (grp in .symRedGroupBy(lapply(preps, `[[`, "support"))) {
    glab <- labels[gi[grp]]
    if (isTRUE(verbose))
      message("curved block {", paste(glab, collapse = ", "), "}: module reduction")
    mr <- .symRedModuleReduce(preps[grp], glab, sd, spy)
    for (sub in .symRedGroupBy(lapply(mr$preps, `[[`, "support"))) {
      subLabs <- sort(unique(glab[unlist(mr$labelsOf[sub])]))
      sol <- .symRedSolveBlock(mr$preps[sub], dPoly, dDarboux, sd, spy, verbose)
      emit(list(labels = subLabs, type = "curved",
                support = sort(unique(unlist(lapply(mr$preps[sub], `[[`, "support")))),
                stage = sol$stage, target = sol$target,
                status = if (length(sol$invariants)) "invariantOnly" else "unresolved",
                invariants = sol$invariants, certificates = sol$certificates,
                reason = sol$reason, moduleCombos = mr$combos,
                transversal = NULL, survivorMeaning = NULL, pins = NULL,
                preps = mr$preps[sub]))
    }
  }
  blocks
}

# ---- verification and trafo assembly ---------------------------------------------

# Exact internal verification. Scaling blocks: the invariant exponents annihilate
# the residual weights over the integers (by construction; asserted). Curved blocks:
# X(I) = 0 proven symbolically -- the Lie derivative of every invariant along every
# generator must cancel() to literally "0", a full proof, not a sample. A failing
# invariant is dropped and the failure recorded (this guards our own algebra; it
# should never fire).
.symRedVerify <- function(blocks, spy) {
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    if (identical(b$type, "scaling")) {
      if (!is.null(b$Wres) && !is.null(b$invExps) && ncol(b$invExps) &&
          nrow(b$Wres)) {
        prod <- b$Wres[, rownames(b$invExps), drop = FALSE] %*% b$invExps
        if (any(prod != 0))
          stop("reduceSymmetry(): internal error -- scaling invariant not in the ",
               "weight kernel.", call. = FALSE)
        blocks[[bi]]$certificates <- c(b$certificates,
          "verified: integer weight annihilation, exact")
      }
      next
    }
    if (!length(b$invariants) || is.null(b$preps)) next
    locals <- .symRedLocals(c(unlist(lapply(b$preps, `[[`, "comps")),
                              b$invariants), spy)
    xiOf <- lapply(b$preps, function(pr)
      lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))
    keep <- logical(length(b$invariants))
    for (l in seq_along(b$invariants)) {
      Ie <- tryCatch(.symRedSympify(gsub("\\^", "**", b$invariants[l]), spy, locals),
                     error = function(e) NULL)
      if (is.null(Ie)) next
      keep[l] <- all(vapply(seq_along(b$preps), function(g) {
        XI <- spy$cancel(spy$together(
          .symRedApplyX(b$preps[[g]], xiOf[[g]], Ie, spy, locals)))
        identical(as.character(XI), "0")
      }, logical(1)))
    }
    if (any(!keep)) {
      blocks[[bi]]$certificates <- c(b$certificates, sprintf(
        "verification DROPPED %d invariant(s): X(I) != 0", sum(!keep)))
      blocks[[bi]]$invariants <- b$invariants[keep]
      if (!any(keep)) blocks[[bi]]$status <- "unresolved"
    } else {
      blocks[[bi]]$certificates <- c(b$certificates,
        "verified: X(I) = 0 exactly (sympy cancel) for every generator")
    }
  }
  blocks
}

# Solve a curved block's invariants into trafo entries. Each invariant is set equal
# to a new parameter carried under the name of the coordinate solved for (the
# survivor-keeps-name idiom: an entry may reference its own name, meaning the OUTER
# parameter, exactly like a log-trafo x = "exp(x)"). The system {I_l = c_l} is
# solved jointly in sympy for distinct linear-degree coordinates, the pins of every
# scaling block substituted in, and only rational solutions with integer powers are
# emitted -- a block whose invariants cannot be solved rationally stays
# invariantOnly, its invariants still reported.
.symRedSolveInvariants <- function(b, pins, spy, coords = character(0)) {
  out <- list(pins = NULL, meaning = NULL, solved = FALSE)
  if (!length(b$invariants)) return(out)
  locals <- .symRedLocals(c(b$invariants, b$support, names(pins), pins), spy)
  Ies <- lapply(b$invariants, function(iv)
    tryCatch(.symRedSympify(gsub("\\^", "**", iv), spy, locals),
             error = function(e) NULL))
  if (any(vapply(Ies, is.null, logical(1)))) return(out)
  # One carrier per invariant: a support coordinate of polynomial degree 1 in the
  # numerator that no other invariant has claimed. Candidates are tried from the
  # END of the coordinate list first: znames orders states before parameters, and a
  # STATE carrier would collide with itself in the trafo (state name = initial
  # condition, the same symbol on the right = new parameter) -- states are left to
  # the gauge pins, which are constants.
  used <- character(0); carriers <- character(0)
  for (Ie in Ies) {
    frac <- spy$fraction(spy$cancel(spy$together(Ie)))
    cand <- setdiff(b$support, used)
    if (length(coords)) cand <- cand[order(-match(cand, coords, nomatch = 0L))]
    pick <- NA_character_
    for (v in cand) {
      dnum <- tryCatch(as.integer(as.character(
        spy$degree(frac[[1]], gen = spy$Symbol(v)))), error = function(e) NA_integer_)
      dden <- tryCatch(as.integer(as.character(
        spy$degree(frac[[2]], gen = spy$Symbol(v)))), error = function(e) NA_integer_)
      # I = c is solvable rationally for v whenever N - c*D is linear in v,
      # i.e. numerator AND denominator are at most linear (either one exactly)
      if (!is.na(dnum) && !is.na(dden) && dnum <= 1L && dden <= 1L &&
          (dnum == 1L || dden == 1L)) {
        pick <- v; break
      }
    }
    if (is.na(pick)) return(out)
    used <- c(used, pick); carriers <- c(carriers, pick)
  }
  # joint solve of I_l = tmp_l for the carriers; tmp_l renamed to the carrier name
  # afterwards (the equation cannot hold the same symbol in both meanings)
  tmpN <- paste0("dModRedC", seq_along(Ies))
  eqs <- lapply(seq_along(Ies), function(l)
    Ies[[l]] - spy$Symbol(tmpN[l]))
  sol <- tryCatch(spy$solve(eqs, lapply(carriers, function(v) spy$Symbol(v)),
                            dict = TRUE), error = function(e) NULL)
  if (is.null(sol) || !length(sol)) return(out)
  sol <- sol[[1]]
  # After the rewrite the output is constant in the non-carrier moved coordinates
  # (varying them IS the orbit motion at fixed invariants), so pinning them is an
  # exact choice of orbit parametrisation -- the curved-orbit reachability caveat
  # applies only BEFORE the rewrite. Pins substitute into the solved entries.
  gauge <- setdiff(b$support, carriers)
  allPins <- c(pins, setNames(rep("1", length(gauge)), gauge))
  entries <- character(0)
  pinPairs <- lapply(names(allPins), function(nm)
    reticulate::tuple(spy$Symbol(nm), .symRedSympify(allPins[[nm]], spy, locals)))
  for (l in seq_along(carriers)) {
    # convert=TRUE turns the solution dict into an R list keyed by str(Symbol)
    e <- if (is.list(sol)) sol[[carriers[l]]]
         else reticulate::py_get_item(sol, spy$Symbol(carriers[l]), silent = TRUE)
    if (is.null(e)) return(out)
    e <- spy$cancel(e$subs(pinPairs))
    ratOK <- tryCatch(isTRUE(e$is_rational_function()), error = function(err) FALSE)
    if (!ratOK) return(out)
    entries[carriers[l]] <- gsub("\\*\\*", "^", as.character(e))
  }
  # tmp_l -> carrier name: the new parameter keeps the solved coordinate's name
  for (l in seq_along(carriers))
    entries <- setNames(gsub(paste0("\\b", tmpN[l], "\\b"), carriers[l], entries),
                        names(entries))
  out$pins <- c(entries, setNames(rep("1", length(gauge)), gauge))
  out$gauge <- gauge
  out$meaning <- setNames(b$invariants, carriers)
  out$solved <- TRUE
  out
}

# The complete P()-ready trafo: identity over every coordinate, the scaling pins,
# and the solved curved entries. Scaling pins are constants and curved entries come
# from a joint solve referencing only outer parameters, so no resolveRecurrence
# rewriting is needed; asserted: no entry references the INNER meaning of a
# different redefined coordinate through a chain a pin would break.
.symRedAssembleTrafo <- function(blocks, coords) {
  vals <- setNames(as.character(coords), coords)
  for (b in blocks) if (!is.null(b$pins)) {
    nm <- intersect(names(b$pins), names(vals))
    vals[nm] <- b$pins[nm]
  }
  as.eqnvec(structure(as.character(vals), names = names(vals)))
}

# ---- result object and print -----------------------------------------------------

.symRedResult <- function(object, blocks, trafo, coords, fixed, settings, call) {
  status <- vapply(blocks, `[[`, character(1), "status")
  labs <- lapply(blocks, `[[`, "labels")
  structure(list(
    method      = object$method,
    coordinates = coords,
    fixed       = fixed,
    blocks      = blocks,
    trafo       = trafo,
    family      = Filter(Negate(is.null), lapply(blocks, function(b)
      if (identical(b$status, "reduced") &&
          (!is.null(b$admissible) || !is.null(b$matroid)))
        list(labels = b$labels, admissible = b$admissible, matroid = b$matroid,
             gaugeNote = b$gaugeNote))),
    removed     = unlist(labs[status %in% c("fixed", "reduced")]),
    remaining   = unlist(labs[!status %in% c("fixed", "reduced")]),
    settings    = settings,
    call        = call), class = "symmetryreduction")
}

# wrap-and-indent for the print method (display columns, not bytes)
.symRedWrap <- function(items, indent, width, sep = ", ") {
  txt <- paste(items, collapse = sep)
  paste(strwrap(txt, width = max(30L, width - nchar(indent)),
                initial = "", prefix = ""), collapse = paste0("\n", indent))
}

#' @export
print.symmetryreduction <- function(x, width = getOption("width"), ...) {
  bar <- strrep("-", 60)
  nDir <- length(x$removed) + length(x$remaining)
  cat(bar, "\n", sep = "")
  cat(sprintf("reduceSymmetry  |  from: %s   directions: %d\n", x$method, nDir))
  cat(bar, "\n", sep = "")
  if (!nDir) { cat("Nothing to reduce.\n"); return(invisible(x)) }
  cat(sprintf("Reduced %d of %d direction%s%s.\n\n",
              length(x$removed), nDir, if (nDir == 1L) "" else "s",
              if (length(x$remaining))
                sprintf(" (%s remaining)", paste(x$remaining, collapse = ", "))
              else ""))
  for (b in x$blocks) {
    hd <- sprintf("%s block {%s} -- %s",
                  if (b$type == "scaling") "Scaling" else "Curved",
                  paste(b$labels, collapse = ", "), b$status)
    cat(hd, "\n", sep = "")
    ind <- "    "
    if (!is.null(b$removedByFixed) && b$removedByFixed > 0L)
      cat(ind, sprintf("fixed removed %d direction%s", b$removedByFixed,
                       if (b$removedByFixed == 1L) "" else "s"),
          if (length(b$redundantFixed))
            paste0("  (redundant: ", paste(b$redundantFixed, collapse = ", "), ")")
          else "", "\n", sep = "")
    if (length(b$moduleCombos))
      cat(ind, "module reduction: ", paste(b$moduleCombos, collapse = ";  "),
          "\n", sep = "")
    if (length(b$invariants))
      cat(ind, "invariants:  ",
          .symRedWrap(b$invariants, paste0(ind, "             "), width), "\n",
          sep = "")
    if (!is.null(b$transversal))
      cat(ind, "transversal: ", paste(paste0(b$transversal, " = 1"), collapse = ", "),
          "\n", sep = "")
    if (length(b$survivorMeaning))
      for (nm in names(b$survivorMeaning))
        cat(ind, "  ", nm, " now carries ", b$survivorMeaning[[nm]], "\n", sep = "")
    for (ct in b$certificates) cat(ind, "[", ct, "]\n", sep = "")
    if (!is.null(b$reason)) cat(ind, "reason: ", b$reason, "\n", sep = "")
    cat("\n")
  }
  nonid <- x$trafo[x$trafo != names(x$trafo)]
  if (length(nonid)) {
    cat("Trafo (non-identity entries; identities cover the other coordinates):\n")
    for (nm in names(nonid)) cat("  ", nm, " = ", nonid[[nm]], "\n", sep = "")
    cat("\n")
  }
  for (f in x$family) {
    cat("Family {", paste(f$labels, collapse = ", "), "}: ", sep = "")
    if (!is.null(f$admissible))
      cat(length(f$admissible), " admissible transversal",
          if (length(f$admissible) == 1L) "" else "s", "\n",
          "  ", .symRedWrap(vapply(f$admissible, function(T)
            paste0("{", paste(T, collapse = ","), "}"), character(1)), "  ", width),
          "\n", sep = "")
    else if (!is.null(f$matroid)) {
      cat("pick one coordinate per row; the picks must be independent\n")
      for (r in f$matroid) cat("  ", .symRedWrap(r, "  ", width), "\n", sep = "")
    }
    cat("  ", .symRedWrap(f$gaugeNote, "  ", width, sep = " "), "\n", sep = "")
  }
  if (length(x$remaining)) {
    cat("Remaining: ", paste(x$remaining, collapse = ", "),
        " -- options: a structural assumption, a separating experiment, or ",
        "prediction profiles.\n", sep = "")
  }
  if (length(x$fixed))
    cat("Note: `fixed` coordinates stay identity entries here -- keep passing ",
        "fixed = c(...) to downstream calls.\n", sep = "")
  invisible(x)
}

# ---- the user-facing function ----------------------------------------------------

#' Constructive removal of detected symmetries
#'
#' Turns the non-identifiable directions reported by [symmetryDetection()] into a
#' parameter reparametrisation. Scaling directions are gauged exactly: the integer
#' weight lattice yields the invariant monomials, and a certified transversal (one
#' coordinate per independent weight row, pinned to 1) removes the directions while
#' every surviving coordinate keeps its name and absorbs an invariant product. When
#' a whole space of transformations exists, the admissible family is reported
#' instead of hiding an arbitrary choice. Curved (polynomial/general) directions go
#' through module reduction and an escalating exact invariant search -- monomial,
#' polynomial up to `dPoly`, rational via Darboux polynomials up to `dDarboux` --
#' and every failed stage leaves a precise negative certificate.
#'
#' @param object A `symmetrydetection` result, from [symmetryDetection()].
#' @param fixed Character vector of coordinates the user pins at known values
#'   beforehand (same semantics as `summary(object, fixed = )`): scaling directions
#'   removed by the fixing drop out of the reparametrisation, and the fixed
#'   coordinates never enter a transversal. Unknown names warn and are ignored.
#' @param dPoly Total-degree cap of the polynomial invariant search (stage 2).
#' @param dDarboux Degree cap of the Darboux factors in the rational invariant
#'   search (stage 3). The cofactor degree is bounded structurally (at most the
#'   generator degree minus one), not by this cap.
#' @param verbose Logical: print progress per block and stage.
#' @param ... Reserved.
#'
#' @return An object of class `symmetryreduction`:
#'   \describe{
#'     \item{`blocks`}{one entry per coupled block of directions: `labels` (the
#'       `X` labels as printed by `print(object)`), `type`, `support`, `stage`
#'       reached, `invariants` (exact strings), `certificates` (positive and
#'       negative), `transversal`/`pins`, `survivorMeaning` (which invariant each
#'       surviving coordinate now carries), `moduleCombos`, `status`
#'       (`"fixed"`, `"reduced"`, `"invariantOnly"` or `"unresolved"`), `reason`.}
#'     \item{`trafo`}{a complete [eqnvec] over all coordinates (identity plus the
#'       transversal pins and solved entries), ready for [P()] or
#'       `symmetryDetection(trafo = )`; `NULL` when nothing was reducible.}
#'     \item{`family`}{the general form: admissible transversals (or their matroid
#'       description) per reduced block, plus the gauge note that pins may be any
#'       nonzero constant.}
#'     \item{`removed`, `remaining`}{direction labels by outcome.}
#'     \item{`coordinates`, `fixed`, `settings`, `call`}{provenance.}
#'   }
#'   `print()` shows the verdict, the per-block reductions, the non-identity trafo
#'   entries, the family and the remaining directions.
#'
#' @seealso [symmetryDetection()]
#' @example inst/examples/reduceSymmetry.R
#' @export
reduceSymmetry <- function(object, fixed = NULL, dPoly = 3L, dDarboux = 2L,
                           verbose = FALSE, ...) {
  if (!inherits(object, "symmetrydetection"))
    stop("reduceSymmetry(): `object` must be a symmetrydetection result.",
         call. = FALSE)
  if (length(list(...)))
    warning("reduceSymmetry(): unused argument(s) ignored.", call. = FALSE)
  .require_ns("reticulate", "reduceSymmetry()")
  .symCall <- match.call()
  settings <- list(dPoly = as.integer(dPoly), dDarboux = as.integer(dDarboux),
                   primes = .symPrimes, verifyPrime = .symVerifyPrime)

  coords <- .symRedCoordinates(object)
  fixed <- unique(as.character(fixed))
  unknown <- setdiff(fixed, coords)
  if (length(unknown))
    warning("reduceSymmetry(): no effect: ", paste(unknown, collapse = ", "),
            " -- not a coordinate of the analysis.", call. = FALSE)

  if (isTRUE(object$identifiable) || !length(object$symmetries))
    return(.symRedResult(object, list(), NULL, coords, fixed, settings, .symCall))

  code_dir <- system.file("code", package = "dMod")
  sysmod <- reticulate::import("sys", convert = TRUE)
  if (!(code_dir %in% sysmod$path)) sysmod$path <- c(code_dir, sysmod$path)
  sd <- reticulate::import("symmetryDetection", convert = TRUE)
  spy <- reticulate::import("sympy", convert = TRUE)

  o <- .symOrdered(object)
  wr <- .symRedWeightRows(o$syms)

  # A scaling whose support overlaps a curved direction cannot be gauged on its
  # own: the curved invariants need not be scale-invariant, so the two reductions
  # would contradict each other (each assumes the other's coordinates fixed). Such
  # scalings join the curved block as generators (xi_i = w_i z_i, polynomial); the
  # scaling stage keeps only the components disjoint from every curved support.
  curvedIdx0 <- setdiff(seq_along(o$syms), wr$rows)
  curvedSupp <- unique(unlist(lapply(o$syms[curvedIdx0], .symCoords)))
  scalRows <- seq_len(nrow(wr$W))
  demoted <- integer(0)
  if (nrow(wr$W) && length(curvedSupp)) {
    repeat {                       # transitive closure: demotion can extend the overlap
      overlaps <- vapply(scalRows, function(r)
        any(colnames(wr$W)[wr$W[r, ] != 0L] %in% curvedSupp), logical(1))
      if (!any(overlaps)) break
      hit <- scalRows[overlaps]
      demoted <- c(demoted, hit)
      curvedSupp <- unique(c(curvedSupp,
        unlist(lapply(hit, function(r) colnames(wr$W)[wr$W[r, ] != 0L]))))
      scalRows <- setdiff(scalRows, hit)
    }
  }

  blocks <- list()
  if (length(scalRows)) {
    Wk <- wr$W[scalRows, , drop = FALSE]
    for (cp in .symRedComponents(Wk)) {
      rows <- scalRows[cp]
      Wb <- wr$W[rows, colSums(wr$W[rows, , drop = FALSE] != 0L) > 0L, drop = FALSE]
      if (isTRUE(verbose))
        message("scaling block {", paste(o$labels[wr$rows[rows]], collapse = ", "),
                "}: ", nrow(Wb), " direction(s)")
      blocks[[length(blocks) + 1L]] <-
        .symRedScalingBlock(Wb, o$labels[wr$rows[rows]], fixed, sd)
    }
  }
  curvedIdx <- sort(c(curvedIdx0, wr$rows[demoted]))
  if (length(curvedIdx))
    blocks <- c(blocks, .symRedCurved(o$syms, curvedIdx, o$labels, fixed,
                                      dPoly, dDarboux, sd, spy, verbose))

  blocks <- .symRedVerify(blocks, spy)

  # solve the curved blocks' invariants into trafo entries, with every scaling pin
  # substituted in so a pinned coordinate cannot re-enter through a solution
  scalPins <- unlist(lapply(blocks, function(b)
    if (identical(b$type, "scaling")) b$pins))
  if (is.null(scalPins)) scalPins <- character(0)
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    if (!identical(b$type, "curved") || !length(b$invariants)) next
    if (!is.null(b$target) && length(b$invariants) < b$target) next  # partial set:
    sol <- .symRedSolveInvariants(b, scalPins, spy, coords)     # gauge pin would be lossy
    if (sol$solved) {
      blocks[[bi]]$pins <- sol$pins
      blocks[[bi]]$transversal <- sol$gauge
      blocks[[bi]]$gaugeNote <- if (length(sol$gauge)) paste0(
        "after the rewrite the output is constant in ",
        paste(sol$gauge, collapse = ", "),
        "; pinning them to 1 is an exact parametrisation choice") else NULL
      blocks[[bi]]$survivorMeaning <- sol$meaning
      blocks[[bi]]$status <- "reduced"
      blocks[[bi]]$certificates <- c(b$certificates,
        "solved rationally: each invariant carried under its solved coordinate's name")
    } else if (isTRUE(verbose))
      message("curved block {", paste(b$labels, collapse = ", "),
              "}: invariants not rationally solvable -- reported as invariantOnly")
  }

  trafo <- .symRedAssembleTrafo(blocks, coords)
  blocks <- lapply(blocks, function(b) { b$preps <- NULL; b$Wres <- NULL
    b$invExps <- NULL; b })
  .symRedResult(object, blocks, trafo, coords, fixed, settings, .symCall)
}
