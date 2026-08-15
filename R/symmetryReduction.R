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
# number, ...). The lookbehind keeps the exponent of "1e-5" from matching; an
# identifier followed by "(" is a function call (exp, log, ...) and must stay
# resolvable, so it never enters the table. The first lookahead pins the match to
# the full identifier so backtracking cannot smuggle a truncated name past the
# call check.
.symRedLocals <- function(strings, spy) {
  ids <- unique(unlist(regmatches(strings,
    gregexpr("(?<![0-9.A-Za-z_])[A-Za-z_][A-Za-z0-9_.]*(?![A-Za-z0-9_.])(?!\\s*\\()",
             strings, perl = TRUE))))
  ids <- setdiff(ids, "")
  reticulate::dict(setNames(lapply(ids, function(x) spy$Symbol(x)), ids))
}

# Sympify memo: every locals table in this file maps every non-call identifier to
# a plain Symbol, so the parse of a given string is the same expression whichever
# table covers it -- a global string-keyed cache is exact. The same strings recur
# across prep, module reduction, Darboux, exp, verify and solve.
.symRedSympifyCache <- new.env(parent = emptyenv())
.symRedSympify <- function(x, spy, locals) {
  key <- as.character(x)
  hit <- .symRedSympifyCache[[key]]
  if (!is.null(hit)) return(hit)
  val <- spy$sympify(key, locals = locals)
  if (length(ls(.symRedSympifyCache, all.names = TRUE)) > 4096L)
    rm(list = ls(.symRedSympifyCache, all.names = TRUE),
       envir = .symRedSympifyCache)
  .symRedSympifyCache[[key]] <- val
  val
}

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
# leaves the fewest nonzeros. A Markowitz-style count cannot serve here: the point of
# the reduction is to discover symbolic cancellations, which appear only as numeric
# zeros AFTER a trial elimination. The fill delta is therefore computed exactly, but
# only over the affected rows and without matrix copies; a pivot in a singleton
# column eliminates nothing and is skipped. Returns the schedule and its final fill.
.symRedReduceSchedule <- function(V, p) {
  k <- nrow(V); n <- ncol(V)
  used <- logical(k); sched <- list()
  for (step in seq_len(k)) {
    nzR <- rowSums(V != 0); nzC <- colSums(V != 0)
    best <- NULL
    for (r in which(!used)) for (cc in which(V[r, ] != 0)) {
      if (nzC[cc] < 2) next
      pinv <- .symInvmod(V[r, cc], p)
      delta <- 0
      for (j in which(V[, cc] != 0)) {
        if (j == r) next
        m <- .symRedMulmodP(V[j, cc], pinv, p)
        delta <- delta + sum((V[j, ] - .symRedMulmodP(V[r, ], m, p)) %% p != 0) - nzR[j]
      }
      if (is.null(best) || delta < best$delta) best <- list(r = r, c = cc, delta = delta)
    }
    if (is.null(best)) break
    r <- best$r; cc <- best$c
    piv <- V[r, cc]
    for (j in seq_len(k)) if (j != r && V[j, cc] != 0) {
      m <- .symRedMulmodP(V[j, cc], .symInvmod(piv, p), p)
      V[j, ] <- (V[j, ] - .symRedMulmodP(V[r, ], m, p)) %% p
    }
    used[r] <- TRUE
    sched[[length(sched) + 1L]] <- c(r, cc)
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
  has <- lapply(seq_len(k), function(r) which(cols %in% names(preps[[r]]$comps)))
  scheds <- lapply(.symPrimes[1:2], function(p) {
    pool <- .symPool()
    pt <- setNames(pool(seq_along(vars) + 7L), vars)
    exprs <- unlist(lapply(seq_len(k), function(r)
      vapply(cols[has[[r]]], function(cc) preps[[r]]$comps[[cc]], character(1))))
    vals <- .symRedEvalBatch(exprs, pt, p, sd)
    vals[is.na(vals)] <- 0
    V <- matrix(0, k, length(cols))
    off <- 0L
    for (r in seq_len(k)) {
      nv <- length(has[[r]])
      if (nv) { V[r, has[[r]]] <- vals[off + seq_len(nv)]; off <- off + nv }
    }
    .symRedReduceSchedule(V, p)
  })
  sched <- scheds[[which.min(vapply(scheds, `[[`, numeric(1), "fill"))]]$sched

  # exact symbolic replay of the schedule, batched into one Python call (the
  # per-entry cancel() loop through reticulate dominated M009-scale blocks);
  # 'combo' tracks the combination coefficients over the original rows
  rr <- sd$moduleReduceReplay(lapply(preps, function(pr) as.list(pr$comps)),
                              as.list(cols), lapply(sched, as.list))

  # Candidates = the original generators (identity combination) plus every reduced
  # row; a greedy pass picks a rank-complete basis with the smallest supports, so a
  # row is only ever replaced when the replacement's support is genuinely smaller.
  cand <- lapply(seq_len(k), function(j)
    list(pr = preps[[j]], srcIdx = j, combo = NULL))
  for (j in seq_len(k)) {
    comps <- unlist(rr$rows[[j]])                      # nonzero entries only
    if (!length(comps)) next                           # row eliminated entirely
    pr <- .symRedGenPrep(list(generator = as.list(comps)), spy)
    if (!pr$ok || !length(pr$support)) next
    Crow <- unlist(rr$combo[[j]])
    srcIdx <- which(Crow != "0")
    if (length(srcIdx) == 1L && srcIdx == j && identical(Crow[[j]], "1"))
      next                                             # unchanged: already in
    combo <- paste0(
      paste0("(", Crow[srcIdx], ")*", labels[srcIdx], collapse = " + "),
      "  [support ", length(pr$support), "]")
    cand[[length(cand) + 1L]] <- list(pr = pr, srcIdx = srcIdx, combo = combo)
  }

  p <- .symPrimes[1]
  pool <- .symPool()
  pt <- setNames(pool(seq_along(vars) + 23L), vars)
  vrow <- function(pr) {
    out <- numeric(length(cols))
    hs <- which(cols %in% names(pr$comps))
    if (length(hs)) {
      v <- .symRedEvalBatch(vapply(cols[hs], function(cc) pr$comps[[cc]],
                                   character(1)), pt, p, sd)
      v[is.na(v)] <- 0
      out[hs] <- v
    }
    out
  }
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

# Batch evaluation of rational expression strings at one integer point mod p: one
# Python round trip for the whole vector (the scalar .symEvalModq re-parses per
# call, which dominates the sampling loops at M009 scale). NA where the
# denominator vanishes mod p or the evaluation fails.
.symRedEvalBatch <- function(exprs, pt, p, sd) {
  if (!length(exprs)) return(numeric(0))
  v <- tryCatch(sd$evalRationalModBatch(as.list(as.character(exprs)),
                                        as.list(names(pt)),
                                        list(as.list(as.numeric(pt))),
                                        as.integer(p)),
                error = function(e) NULL)
  if (is.null(v)) return(rep(NA_real_, length(exprs)))
  out <- as.numeric(unlist(v[[1]]))
  out[out < 0] <- NA_real_
  out
}

# The log-derivative rows of a block at one random point mod p: one row per
# generator, entries eta_i = xi_i(z)/z_i over the moved coordinates. NULL when an
# evaluation fails (redraw the point).
.symRedEtaRows <- function(preps, vars, pt, p, sd) {
  idx <- lapply(preps, function(pr) which(vars %in% names(pr$comps)))
  exprs <- unlist(lapply(seq_along(preps), function(g)
    vapply(vars[idx[[g]]], function(v) preps[[g]]$comps[[v]], character(1))))
  vals <- .symRedEvalBatch(exprs, pt, p, sd)
  if (anyNA(vals)) return(NULL)
  rows <- matrix(0, length(preps), length(vars))
  off <- 0L
  for (g in seq_along(preps)) {
    for (k in seq_along(idx[[g]])) {
      i <- idx[[g]][k]
      rows[g, i] <- .symMulmod(vals[off + k], .symInvmod(pt[[vars[i]]] %% p, p), p)
    }
    off <- off + length(idx[[g]])
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
    # exponents are tiny and eta < p < 2^31, so the inner products stay far below
    # 2^53 and the plain matrix product is exact; only the outer multiply needs
    # the split-multiply (elementwise-safe for equal-length vectors)
    inner <- (expts %*% t(eta)) %% p                       # N x nGen
    t(vapply(seq_len(nrow(eta)), function(g)
      .symRedMulmodP(mv, inner[, g], p), numeric(N)))
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

# Exact integer sampling rows of a cofactor matrix: one row per (generator, point),
# column j evaluating cols[[j]][[g]] (sympy expressions or strings, rational in
# `vars`; values through the exact-Fraction batch evaluator, one Python call per
# generator). Each row is lcm-cleared to integers; rows with a failed evaluation
# or entries >= 2^31 are dropped. Serves the Darboux kernel [lambda], the exp
# stage [lambda | mu] and the integrating-factor solve [lambda | mu | div X].
.symRedCofactorRows <- function(cols, vars, nGen, nMon, sd) {
  rows <- list()
  for (g in seq_len(nGen)) {
    exprs <- vapply(cols, function(cl) as.character(cl[[g]]), character(1))
    pts <- lapply(seq_len(nMon), function(t)
      as.list(.symPool()(seq_along(vars) + 101L * t + 7L * g)))
    v <- tryCatch(sd$evalRationalBatch(as.list(exprs), as.list(vars), pts),
                  error = function(e) NULL)
    if (is.null(v)) next
    for (t in seq_len(nMon)) {
      num <- suppressWarnings(as.numeric(unlist(v$num[[t]])))
      den <- suppressWarnings(as.numeric(unlist(v$den[[t]])))
      if (anyNA(num) || anyNA(den) || any(den == 0)) next
      L <- Reduce(function(a, b) a * b / .symRedGcd(a, b), unique(den), 1)
      row <- round(num * (L / den))
      if (max(abs(row)) < 2^31) rows[[length(rows) + 1L]] <- as.integer(row)
    }
  }
  rows
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
  # The extactic basis runs over the MOVED coordinates only, with every
  # coefficient symbol treated as a constant of the ground field: Pereira's
  # theorem holds over any characteristic-0 field, so completeness is "factor
  # degree <= d in the moved coordinates, arbitrary degree in the parameters" --
  # a stronger statement than the all-variables basis, at a fraction of the size
  # (C(nMoved + d, d) instead of C(nvars + d, d)); module-reduced blocks with a
  # handful of moved coordinates now pass the cap that used to skip them.
  Nbasis <- choose(length(moved) + dDarboux, dDarboux)
  # the symbolic determinant's cost is driven by the entry degrees, which grow to
  # ~ dDarboux + (N-1)*(D-1) in the last extactic row -- cap that, not just N
  Dmax <- max(1L, vapply(preps, function(pr)
    if (is.na(pr$degree)) 3L else pr$degree, integer(1)))
  entryDeg <- dDarboux + (Nbasis - 1L) * (Dmax - 1L)
  extactic <- dDarboux >= 1L && Nbasis <= .symRedExtacticCap &&
              entryDeg <= .symRedExtacticCap
  degenerate <- FALSE
  if (extactic) {
    expts <- .symMonoTable(length(moved), as.integer(dDarboux))
    basis <- lapply(seq_len(nrow(expts)), function(k)
      .symRedSympify(.symRedMonoString(expts[k, ], moved), spy, locals))
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

  # keep candidates that are Darboux for EVERY generator, with their cofactors and
  # proportionality dedup -- batched into one Python call
  dc <- sd$darbouxCofactors(lapply(cands, as.character),
                            lapply(preps, function(pr) as.list(pr$comps)))
  Ps <- lapply(dc$Ps, function(s) .symRedSympify(s, spy, locals))
  lams <- lapply(dc$lams, function(lr)
    lapply(lr, function(s) .symRedSympify(s, spy, locals)))

  coverage <- if (degenerate)
      "extactic degenerate (E identically 0): coordinate and xi factors only"
    else if (extactic)
      sprintf(paste0("extactic complete for factor degree <= %d in the moved ",
                     "coordinates"), dDarboux)
    else sprintf(paste0("extactic skipped (moved-coordinate basis %d, projected ",
                        "entry degree %d, cap %d): coordinate and xi factors only"),
                 Nbasis, entryDeg, .symRedExtacticCap)
  certBase <- sprintf("Darboux stage (%d candidate factors; %s)", length(Ps), coverage)
  # Ps/lams/coverage travel on every exit path: the exp stage consumes the factor
  # basis and its cofactors even when no purely rational combination exists here.
  fail <- function(cert) list(ok = FALSE, invariants = character(0), cert = cert,
                              Ps = Ps, lams = lams, coverage = coverage)
  if (length(Ps) < 2L)
    return(fail(paste0("no rational invariant from Darboux factors -- ", certBase)))

  # integer kernel of the cofactor matrix, sampled exactly at integer points
  nMon <- choose(length(vars) + max(1L, dDarboux), length(vars)) + 3L
  rows <- .symRedCofactorRows(lams, vars, length(preps), nMon, sd)
  if (!length(rows))
    return(fail(paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  K <- .symRedIntKernel(do.call(rbind, rows), sd)
  if (!ncol(K))
    return(fail(paste0("no rational invariant with Darboux factors of degree <= ",
                       dDarboux, " -- ", certBase)))
  Pstr <- vapply(Ps, function(P) gsub("\\*\\*", "^", as.character(P)), character(1))
  inv <- vapply(seq_len(ncol(K)), function(j) {
    nz <- which(K[, j] != 0L)
    paste(vapply(nz, function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      if (K[l, j] == 1L) base else paste0(base, "^", K[l, j])
    }, character(1)), collapse = "*")
  }, character(1))
  list(ok = TRUE, invariants = inv, cert = certBase,
       Ps = Ps, lams = lams, coverage = coverage)
}

# denominator-candidate cap of the exp stage (h = 1, all P_j, then all P_j^2)
.symRedExpDenCap <- 13L
# unknown-count budget per exp-stage linear system (g- plus mu-coefficients)
.symRedExpSizeCap <- 250L

# Stage 4: exponential factors E = exp(g/h) with polynomial cofactor mu = X(g/h)
# (Prelle-Singer: every Liouvillian first integral is of Darboux form once these
# are admitted). For a fixed denominator h = prod P_j^{k_j} built from the Darboux
# factors (cofactor Lambda_h = sum k_j lambda_j) the condition
#   X(g) - g*Lambda_h = mu*h
# is LINEAR in the coefficients of g (total degree <= dExp, shared across the
# block's generators) and of the per-generator cofactors mu (structural bound
# deg mu <= dExp + D - 1 - deg h), so the same exact modular-nullspace sampling as
# the polynomial stage applies. Every candidate is then certified symbolically
# (the defining identity must cancel() to 0 for every generator), and invariants
# combine Darboux and exponential factors through the integer kernel of the
# extended cofactor matrix [lambda | mu]:
#   I = prod P_j^{n_j} * exp(sum m_k g_k/h_k),  sum n lambda + sum m mu = 0.
.symRedExpFactors <- function(preps, darb, dExp, sd, spy, verbose = FALSE) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  K <- length(preps)
  D <- max(1L, vapply(preps, function(pr)
    if (is.na(pr$degree)) 1L else pr$degree, integer(1)))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xiOf <- lapply(preps, function(pr)
    lapply(pr$comps, function(x) .symRedSympify(x, spy, locals)))
  zero <- spy$Integer(0L)

  # denominator candidates: 1, then every Darboux factor, then their squares (so
  # the cap trims the squares first)
  hC <- list(list(h = spy$Integer(1L), cof = rep(list(zero), K), deg = 0L))
  hC2 <- list()
  for (j in seq_along(darb$Ps)) {
    P <- darb$Ps[[j]]
    dj <- tryCatch(as.integer(as.character(spy$total_degree(P))),
                   error = function(e) NA_integer_)
    if (is.na(dj) || dj < 1L) next
    hC[[length(hC) + 1L]] <- list(h = P, cof = darb$lams[[j]], deg = dj)
    hC2[[length(hC2) + 1L]] <- list(h = spy$expand(P * P),
      cof = lapply(darb$lams[[j]], function(l) spy$Integer(2L) * l), deg = 2L * dj)
  }
  hC <- c(hC, hC2)
  capped <- length(hC) > .symRedExpDenCap
  if (capped) hC <- hC[seq_len(.symRedExpDenCap)]

  # g-ansatz: moved monomials plus the constant (exp(c/h) is a genuine factor)
  exptsG <- .symMonoTable(length(vars), as.integer(dExp))
  keepG <- rowSums(exptsG[, match(moved, vars), drop = FALSE]) > 0L |
    rowSums(exptsG) == 0L
  exptsG <- exptsG[keepG, , drop = FALSE]
  Ng <- nrow(exptsG)

  factors <- list()
  seen <- character(0)
  sizeSkipped <- 0L
  for (hc in hC) {
    # deg mu <= D - 1 is the DEFINING bound of an exponential factor (Christopher/
    # Llibre), not just a structural one -- and it is lossless for the combination
    # step: higher-degree cofactors could only cancel among themselves, and the
    # corresponding g-sum lies in this ansatz with a small cofactor already. It
    # also keeps the h = 1 system from being solved by every g.
    degMu <- min(D - 1L, dExp + D - 1L - hc$deg)
    exptsMu <- if (degMu >= 0L) .symMonoTable(length(vars), degMu) else
      matrix(0L, 0L, length(vars))
    Nmu <- nrow(exptsMu)
    n <- Ng + K * Nmu
    if (n > .symRedExpSizeCap) { sizeSkipped <- sizeSkipped + 1L; next }
    hStr <- as.character(hc$h)
    lamStr <- vapply(hc$cof, as.character, character(1))
    rowsAt <- function(pt, p) {
      eta <- .symRedEtaRows(preps, vars, pt, p, sd)
      if (is.null(eta)) return(NULL)
      hl <- .symRedEvalBatch(c(hStr, lamStr), pt, p, sd)
      if (anyNA(hl)) return(NULL)
      mvG <- symMonoResidues(exptsG, as.integer(pt %% p), p)
      mvMu <- if (Nmu) symMonoResidues(exptsMu, as.integer(pt %% p), p) else numeric(0)
      inner <- (exptsG %*% t(eta)) %% p                 # Ng x K
      muNeg <- if (Nmu) (p - .symRedMulmodP(mvMu, hl[1], p)) %% p else numeric(0)
      t(vapply(seq_len(K), function(g) {
        gcols <- (.symRedMulmodP(mvG, inner[, g], p) -
                    .symRedMulmodP(mvG, hl[1L + g], p)) %% p
        mucols <- numeric(K * Nmu)
        if (Nmu) mucols[(g - 1L) * Nmu + seq_len(Nmu)] <- muNeg
        c(gcols, mucols)
      }, numeric(n)))
    }
    ns <- .symRedModularNullspace(rowsAt, n, vars)
    if (is.null(ns$basis) || !ncol(ns$basis)) next
    for (jc in seq_len(ncol(ns$basis))) {
      gvec <- ns$basis[seq_len(Ng), jc]
      if (all(gvec == 0L)) next
      gExpr <- .symRedSympify(.symRedPolyString(gvec, exptsG, vars), spy, locals)
      q <- spy$cancel(gExpr / hc$h)
      if (isTRUE(q$is_number)) next                     # exp(const)
      if (hc$deg > 0L &&
          isTRUE(tryCatch(q$is_polynomial(), error = function(e) FALSE)))
        next                                            # h | g: the h = 1 run owns it
      mu <- lapply(seq_len(K), function(g) {
        if (!Nmu) return(zero)
        mv <- ns$basis[Ng + (g - 1L) * Nmu + seq_len(Nmu), jc]
        if (all(mv == 0L)) zero else
          .symRedSympify(.symRedPolyString(mv, exptsMu, vars), spy, locals)
      })
      # exact certification of the defining identity, per generator
      okAll <- all(vapply(seq_len(K), function(g) {
        Xg <- .symRedApplyX(preps[[g]], xiOf[[g]], gExpr, spy, locals)
        identical(as.character(spy$cancel(spy$together(
          Xg - gExpr * hc$cof[[g]] - mu[[g]] * hc$h))), "0")
      }, logical(1)))
      if (!okAll) next
      key <- as.character(q)
      if (key %in% seen) next
      seen <- c(seen, key)
      factors[[length(factors) + 1L]] <- list(g = gExpr, h = hc$h, mu = mu)
    }
  }

  certBase <- sprintf(paste0("exp stage (numerator degree <= %d over %d ",
                             "denominator candidate(s)%s%s): %d exponential ",
                             "factor(s)"),
                      dExp, length(hC), if (capped) ", capped" else "",
                      if (sizeSkipped) sprintf(", %d skipped over the %d-unknown budget",
                                               sizeSkipped, .symRedExpSizeCap) else "",
                      length(factors))
  if (!length(factors))
    return(list(ok = FALSE, invariants = character(0), factors = list(),
                cert = paste0("no exponential factor with numerator degree <= ",
                              dExp, " -- ", certBase)))
  if (isTRUE(verbose)) message("  exp: ", length(factors), " factor(s) certified")

  # integer kernel of the extended cofactor matrix [lambda | mu]
  J <- length(darb$Ps)
  cols <- c(darb$lams, lapply(factors, `[[`, "mu"))
  nMon <- choose(length(vars) + max(1L, dExp + D - 1L), length(vars)) + 3L
  rows <- .symRedCofactorRows(cols, vars, K, nMon, sd)
  if (!length(rows))
    return(list(ok = FALSE, invariants = character(0), factors = factors,
                cert = paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  Kk <- .symRedIntKernel(do.call(rbind, rows), sd)
  useCol <- if (ncol(Kk)) which(vapply(seq_len(ncol(Kk)), function(j)
    any(Kk[J + seq_along(factors), j] != 0L), logical(1))) else integer(0)
  if (!length(useCol))
    return(list(ok = FALSE, invariants = character(0), factors = factors,
                cert = paste0("no invariant combining the exponential factors -- ",
                              certBase)))
  Pstr <- vapply(darb$Ps, function(P) gsub("\\*\\*", "^", as.character(P)),
                 character(1))
  inv <- vapply(useCol, function(j) {
    w <- Kk[, j]
    parts <- vapply(which(w[seq_len(J)] != 0L), function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      if (w[l] == 1L) base else paste0(base, "^", w[l])
    }, character(1))
    arg <- Reduce(`+`, lapply(which(w[J + seq_along(factors)] != 0L), function(e2)
      spy$Integer(w[J + e2]) * factors[[e2]]$g / factors[[e2]]$h))
    argStr <- gsub("\\*\\*", "^",
                   as.character(spy$cancel(spy$together(arg))))
    paste(c(parts, paste0("exp(", argStr, ")")), collapse = "*")
  }, character(1))
  list(ok = TRUE, invariants = inv, cert = certBase, factors = factors)
}

# Stage 5 (single generator, exactly two moved coordinates): integrating factor
# in Darboux form. Singer: a Liouvillian first integral exists iff an integrating
# factor M = prod P^n * exp(sum m g/h) does, with RATIONAL exponents solving the
# inhomogeneous cofactor equation sum n lambda + sum m mu = -div X. That is one
# extra column in the existing sampling kernel: integer vectors of
# [lambda | mu | div X] with a nonzero last entry, normalised by it. The first
# integral itself comes from sympy quadrature of dI = M*(xi2 dz1 - xi1 dz2) and
# is accepted only when X(I) = 0 checks exactly.
.symRedIntFactor <- function(preps, darb, exFactors, sd, spy, verbose = FALSE) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  locals <- .symRedLocals(c(unlist(lapply(preps, `[[`, "comps")), vars), spy)
  xi <- lapply(preps[[1]]$comps[moved], function(x) .symRedSympify(x, spy, locals))
  divX <- spy$expand(Reduce(`+`, lapply(moved, function(v)
    spy$diff(xi[[v]], spy$Symbol(v)))))
  J <- length(darb$Ps); E <- length(exFactors)
  cols <- c(darb$lams, lapply(exFactors, `[[`, "mu"), list(list(divX)))
  nMon <- choose(length(vars) + 2L, length(vars)) + 3L
  rows <- .symRedCofactorRows(cols, vars, 1L, nMon, sd)
  certBase <- sprintf("integrating-factor stage (%d Darboux + %d exponential factor(s))",
                      J, E)
  if (!length(rows))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("cofactor sampling failed (inconclusive) -- ", certBase)))
  Kk <- .symRedIntKernel(do.call(rbind, rows), sd)
  last <- J + E + 1L
  useCol <- if (ncol(Kk)) which(Kk[last, ] != 0L) else integer(0)
  if (!length(useCol))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("no Darboux-form integrating factor -- ", certBase)))
  Pstr <- vapply(darb$Ps, function(P) gsub("\\*\\*", "^", as.character(P)),
                 character(1))
  inv <- character(0)
  for (j in useCol) {
    w <- Kk[, j]; wL <- w[last]
    parts <- vapply(which(w[seq_len(J)] != 0L), function(l) {
      base <- if (grepl("[+*/ -]", Pstr[l])) paste0("(", Pstr[l], ")") else Pstr[l]
      paste0(base, "^(", w[l], "/", wL, ")")
    }, character(1))
    expPart <- if (E && any(w[J + seq_len(E)] != 0L)) {
      arg <- Reduce(`+`, lapply(which(w[J + seq_len(E)] != 0L), function(e2)
        spy$Rational(w[J + e2], wL) * exFactors[[e2]]$g / exFactors[[e2]]$h))
      paste0("exp(", as.character(spy$cancel(spy$together(arg))), ")")
    } else character(0)
    Mstr <- paste(c(parts, expPart), collapse = "*")
    if (!nchar(Mstr)) Mstr <- "1"
    Iv <- tryCatch(sd$integratingFactorIntegral(
      gsub("\\^", "**", Mstr), as.character(xi[[moved[1]]]),
      as.character(xi[[moved[2]]]), moved[1], moved[2]),
      error = function(e) NULL)
    if (!is.null(Iv)) {
      inv <- c(inv, gsub("\\*\\*", "^", Iv))
      break
    }
  }
  if (!length(inv))
    return(list(ok = FALSE, invariants = character(0),
                cert = paste0("integrating factor found, quadrature failed -- ",
                              certBase)))
  if (isTRUE(verbose)) message("  intfactor: first integral by quadrature")
  list(ok = TRUE, invariants = inv, cert = paste0(certBase, ": solved by quadrature"))
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
  pool <- as.numeric(.symPool()(seq_along(vars) + 13L))
  # exp-carrying invariants overflow at prime-sized points; small distinct values
  # keep every gradient finite without losing genericity (rank only)
  if (any(grepl("exp(", invs, fixed = TRUE))) pool <- 1 + pool / (max(pool) + 1)
  pt <- setNames(pool, vars)
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
    if (!all(is.finite(gr))) next
    Jc <- rbind(J, gr)
    if (qr(Jc, tol = 1e-9)$rank > nrow(J)) { J <- Jc; sel <- c(sel, iv) }
  }
  sel
}

# run the invariant-search stages on one (already module-reduced) sub-block of
# generators -- monomial, then polynomial <= dPoly, then Darboux <= dDarboux,
# then exponential factors with numerator degree <= dExp.
.symRedSolveBlock <- function(preps, dPoly, dDarboux, dExp, sd, spy,
                              verbose = FALSE) {
  moved <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  vars <- sort(unique(c(moved, unlist(lapply(preps, `[[`, "vars")))))
  target <- .symRedBlockCorank(preps)
  certs <- character(0); found <- character(0); stage <- "monomial"
  # an invariant free of every moved coordinate is trivially constant on the
  # orbits and must never count toward the target (a Darboux factor with
  # cofactor 0, e.g. a pure parameter sum, would otherwise fake a complete
  # invariant set and an exactness claim); monomial/polynomial stages exclude
  # such invariants by construction, the factor stages cannot
  movedRe <- paste0("(?<![0-9A-Za-z_.])(", paste(moved, collapse = "|"),
                    ")(?![0-9A-Za-z_.])")
  pick <- function(cand) .symRedIndependentSet(
    cand[grepl(movedRe, cand, perl = TRUE)], vars, target)

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
  if (length(found) < target) {
    if (dExp > 0L) {
      if (isTRUE(verbose)) message("  darboux: ", length(found), "/", target,
                                   "; escalating to exponential factors")
      stage <- "exp"
      ex <- .symRedExpFactors(preps, dar, dExp, sd, spy, verbose)
      certs <- c(certs, ex$cert)
      if (ex$ok) found <- pick(c(found, ex$invariants))
      if (length(found) < target && length(preps) == 1L &&
          length(unique(unlist(lapply(preps, `[[`, "support")))) == 2L) {
        stage <- "intfactor"
        itf <- .symRedIntFactor(preps, dar, ex$factors, sd, spy, verbose)
        certs <- c(certs, itf$cert)
        # a 2-coordinate single-generator block has target 1 and `found` is empty
        # here, so the quadrature integral (already proven X(I) = 0 in Python) is
        # appended directly -- pick()'s numeric gradient cannot evaluate the
        # atan/Abs forms this stage may produce
        if (itf$ok) found <- c(found, itf$invariants[1])
      }
    } else
      certs <- c(certs, "exp stage skipped (dExp = 0)")
  }
  list(stage = stage, invariants = found, target = target, certificates = certs,
       reason = if (length(found) >= target) NULL
         else sprintf(paste0("found %d of %d independent invariants (poly degree ",
                             "<= %d, Darboux degree <= %d, exp numerator degree ",
                             "<= %d); gauge pinning needs the full set"),
                      length(found), target, dPoly, dDarboux, dExp))
}

# Fixed-coordinate elimination for one curved group: over the function field the
# admissible generators are the combinations whose xi vanishes on every fixed
# coordinate (a fixed coordinate cannot move). Exact symbolic Gaussian elimination
# on the fixed columns; each pivot row cannot combine away from the fixed
# coordinates and is removed by the fixing -- the curved mirror of
# .symRedScalingFixed. Kept rows are exactly zero on every fixed column (the
# elimination is symbolic), so the fixed coordinates leave their support and can
# never become carriers or gauge pins. `sources` tracks which original directions
# enter each kept combination.
.symRedFixedReduce <- function(preps, fixed, spy) {
  k <- length(preps)
  cols <- sort(unique(unlist(lapply(preps, `[[`, "support"))))
  fixedCols <- intersect(fixed, cols)
  if (!length(fixedCols) || !k)
    return(list(preps = preps, sources = as.list(seq_len(k)),
                removed = 0L, fixedCols = character(0)))
  locals <- .symRedLocals(unlist(lapply(preps, `[[`, "comps")), spy)
  E <- lapply(preps, function(pr) {
    row <- setNames(vector("list", length(cols)), cols)
    for (cc in cols) row[[cc]] <- .symRedSympify(
      if (cc %in% names(pr$comps)) pr$comps[[cc]] else "0", spy, locals)
    row
  })
  sources <- as.list(seq_len(k))
  zero <- function(e) identical(as.character(e), "0")
  usedPivot <- logical(k)
  for (fc in fixedCols) {
    piv <- NA_integer_
    for (r in seq_len(k)) if (!usedPivot[r] && !zero(E[[r]][[fc]])) {
      piv <- r; break
    }
    if (is.na(piv)) next
    for (j in seq_len(k)) {
      if (j == piv || zero(E[[j]][[fc]])) next
      m <- spy$cancel(E[[j]][[fc]] / E[[piv]][[fc]])
      for (c2 in cols) E[[j]][[c2]] <- spy$cancel(E[[j]][[c2]] - m * E[[piv]][[c2]])
      sources[[j]] <- sort(unique(c(sources[[j]], sources[[piv]])))
    }
    usedPivot[piv] <- TRUE
  }
  outP <- list(); outS <- list()
  for (j in which(!usedPivot)) {
    comps <- vapply(cols, function(cc) as.character(E[[j]][[cc]]), character(1))
    names(comps) <- cols
    pr <- .symRedGenPrep(list(generator = as.list(comps[comps != "0"])), spy)
    if (!pr$ok || !length(pr$support)) next          # row eliminated entirely
    outP[[length(outP) + 1L]] <- pr
    outS[[length(outS) + 1L]] <- sources[[j]]
  }
  list(preps = outP, sources = outS, removed = sum(usedPivot),
       removedIdx = which(usedPivot), fixedCols = fixedCols)
}

# The curved pipeline: prep every direction with a closed-form generator, group by
# shared support, fixed-eliminate and module-reduce each group, re-split (a
# reduced row may decouple), then run the invariant stages per sub-block.
# Support-only directions and failed preps become unresolved singleton blocks.
.symRedCurved <- function(syms, idx, labels, fixed, dPoly, dDarboux, dExp, sd,
                          spy, verbose = FALSE) {
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
    fx <- .symRedFixedReduce(preps[grp], fixed, spy)
    fixCert <- if (fx$removed > 0L)
      sprintf(paste0("fixed-coordinate elimination: %d of %d direction(s) cannot ",
                     "combine away from {%s}"), fx$removed, length(grp),
              paste(fx$fixedCols, collapse = ", "))
    if (fx$removed > 0L && isTRUE(verbose))
      message("curved block {", paste(glab, collapse = ", "), "}: fixing removed ",
              fx$removed, " direction(s)")
    # pivot directions that survive nowhere get their own "fixed" block, so the
    # direction bookkeeping (removed/remaining) stays complete
    fixedLabs <- if (fx$removed > 0L)
      setdiff(glab[fx$removedIdx], glab[unlist(fx$sources)]) else character(0)
    if (length(fixedLabs) || !length(fx$preps)) {
      fl <- if (length(fx$preps)) fixedLabs else glab
      emit(list(labels = fl, type = "curved",
                support = sort(unique(unlist(lapply(preps[grp], `[[`, "support")))),
                stage = "fixed", status = "fixed",
                removedByFixed = fx$removed, invariants = character(0),
                certificates = fixCert, reason = NULL, moduleCombos = NULL,
                transversal = NULL, survivorMeaning = NULL, pins = NULL,
                gaugeNote = "removed entirely by the user's fixed coordinates"))
    }
    if (!length(fx$preps)) next
    glabK <- vapply(fx$sources, function(s)
      if (length(s) == 1L) glab[s]
      else paste0("(", paste(glab[s], collapse = "+"), ")"), character(1))
    if (isTRUE(verbose))
      message("curved block {", paste(glab, collapse = ", "), "}: module reduction")
    mr <- .symRedModuleReduce(fx$preps, glabK, sd, spy)
    for (sub in .symRedGroupBy(lapply(mr$preps, `[[`, "support"))) {
      subLabs <- sort(unique(glab[unlist(fx$sources[unlist(mr$labelsOf[sub])])]))
      sol <- .symRedSolveBlock(mr$preps[sub], dPoly, dDarboux, dExp, sd, spy,
                               verbose)
      emit(list(labels = subLabs, type = "curved",
                support = sort(unique(unlist(lapply(mr$preps[sub], `[[`, "support")))),
                stage = sol$stage, target = sol$target,
                status = if (length(sol$invariants)) "invariantOnly" else "unresolved",
                removedByFixed = if (fx$removed > 0L && !length(fixedLabs))
                  fx$removed,
                invariants = sol$invariants,
                certificates = c(fixCert, sol$certificates),
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
.symRedVerify <- function(blocks, sd) {
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
    keep <- as.logical(unlist(tryCatch(sd$verifyInvariants(
      as.list(gsub("\\^", "**", b$invariants)),
      lapply(b$preps, function(pr) as.list(pr$comps))),
      error = function(e) rep(FALSE, length(b$invariants)))))
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

# iterate over a python iterable, tolerating reticulate's auto-conversion to list
.symRedIter <- function(x, f)
  tryCatch(reticulate::iterate(x, f), error = function(e) lapply(x, f))

# Whether a solved entry stays inside the emitted trafo language: rational
# operations, Rational powers, exp and log -- nothing else (LambertW & friends
# mean the carrier choice was wrong, not that the block is reducible). Returns
# "no", "yes" or "root" (yes, with a non-integer power needing a branch note).
.symRedEntryClass <- function(e, spy) {
  fnames <- tryCatch(unlist(.symRedIter(e$atoms(spy$Function), function(f)
    as.character(reticulate::py_get_attr(reticulate::py_get_attr(f, "func"),
                                         "__name__")))), error = function(err) "?")
  if (length(fnames) && !all(fnames %in% c("exp", "log"))) return("no")
  root <- FALSE
  pows <- tryCatch(.symRedIter(e$atoms(spy$Pow), function(pw) pw),
                   error = function(err) list())
  for (pw in pows) {
    ex <- pw$exp
    if (!isTRUE(ex$is_Rational)) return("no")
    if (!isTRUE(ex$is_Integer)) root <- TRUE
  }
  if (root) "root" else "yes"
}

# Positivity certificate: TRUE when the expression is a ratio of polynomials
# whose coefficients share one sign -- such an entry maps ANY positive outer
# point to a positive inner value, so the chart covers the whole positive
# orthant. Sufficient, not necessary.
.symRedPosCert <- function(e, spy) {
  fr <- tryCatch(spy$fraction(spy$cancel(spy$together(e))),
                 error = function(err) NULL)
  if (is.null(fr)) return(FALSE)
  sgn <- function(pp) {
    ex <- tryCatch(spy$expand(pp), error = function(err) NULL)
    if (is.null(ex)) return(0L)
    if (isTRUE(ex$is_number)) {
      v <- suppressWarnings(as.numeric(as.character(spy$Float(ex))))
      if (is.na(v) || v == 0) return(0L)
      return(if (v > 0) 1L else -1L)
    }
    cf <- tryCatch(unlist(.symRedIter(spy$Poly(ex)$coeffs(), function(cc)
      suppressWarnings(as.numeric(as.character(spy$Float(cc)))))),
      error = function(err) NULL)
    if (is.null(cf) || !length(cf) || anyNA(cf)) return(0L)
    if (all(cf > 0)) 1L else if (all(cf < 0)) -1L else 0L
  }
  s1 <- sgn(fr[[1]])
  s1 != 0L && s1 == sgn(fr[[2]])
}

# Gauge-section candidates for a curved block: monomial balances "m1 = m2" over
# the support, simple ones first. A balance section can intersect every positive
# orbit (a constant pin cannot -- the curved orbit may not reach it); which one
# actually yields a positive chart is decided by the certificate above.
.symRedSectionCands <- function(support) {
  mons <- data.frame(m = support, d = 1L, stringsAsFactors = FALSE)
  if (length(support) > 1L) {
    pr <- utils::combn(support, 2L)
    mons <- rbind(mons, data.frame(m = paste(pr[1, ], pr[2, ], sep = "*"),
                                   d = 2L, stringsAsFactors = FALSE))
  }
  if (nrow(mons) < 2L) return(list())
  idx <- utils::combn(nrow(mons), 2L)
  out <- lapply(seq_len(ncol(idx)), function(j)
    c(mons$m[idx[1, j]], mons$m[idx[2, j]]))
  out[order(mons$d[idx[1, ]] + mons$d[idx[2, ]])]
}

# Monotone-transversality pre-filter for a balance candidate m1 = m2. Along an
# orbit the ratio r = m1/m2 moves as X(log r) = (X(m1)*m2 - m1*X(m2))/(m1*m2);
# when that numerator is sign-pure for every generator (zero allowed, nonzero
# for at least one), r is strictly monotone on every positive orbit segment:
# the section is crossed at most once and the solve has a single branch. Mixed
# signs allow tangencies and double crossings -- no global positive chart --
# so failing candidates are dropped before any solve is attempted.
.symRedSectionMonotone <- function(pr, preps, xiOf, spy, locals) {
  ms <- lapply(pr, function(m) .symRedSympify(gsub("\\^", "**", m), spy, locals))
  moved <- FALSE
  for (g in seq_along(preps)) {
    Xm <- lapply(ms, function(m)
      .symRedApplyX(preps[[g]], xiOf[[g]], m, spy, locals))
    d <- tryCatch(spy$expand(Xm[[1]] * ms[[2]] - ms[[1]] * Xm[[2]]),
                  error = function(e) NULL)
    if (is.null(d)) return(FALSE)
    if (isTRUE(d$is_zero)) next
    if (isTRUE(d$is_number)) { moved <- TRUE; next }
    cf <- tryCatch(unlist(.symRedIter(spy$Poly(d)$coeffs(), function(cc)
      suppressWarnings(as.numeric(as.character(spy$Float(cc)))))),
      error = function(e) NULL)
    if (is.null(cf) || !length(cf) || anyNA(cf)) return(FALSE)
    if (!(all(cf > 0) || all(cf < 0))) return(FALSE)
    moved <- TRUE
  }
  moved
}

# Solve a curved block's invariants into trafo entries. Each invariant is set equal
# to a new parameter carried under the name of the coordinate solved for (the
# survivor-keeps-name idiom: an entry may reference its own name, meaning the OUTER
# parameter, exactly like a log-trafo x = "exp(x)"). The system {I_l = c_l} is
# solved jointly in sympy; per invariant a carrier coordinate is chosen that makes
# the solve exact in the trafo language:
#   - "power": v enters the rational part only in degrees {0, k} for one k >= 1
#     (k = 1 is the plain linear case) -> N - c*D = alpha*v^k + beta, a k-th root;
#   - "exp": v sits only inside the single exponential factor, linearly in its
#     argument's numerator -> a log solve.
# Chained solutions (one carrier's entry referencing another) are resolved by
# substitution, the pins of every scaling block are substituted in, and only
# solutions built from rational operations, Rational powers, exp and log are
# emitted -- a block whose invariants cannot be solved in that language stays
# invariantOnly, its invariants still reported.
.symRedSolveInvariants <- function(b, pins, spy, coords = character(0)) {
  out <- list(pins = NULL, meaning = NULL, solved = FALSE)
  if (!length(b$invariants)) return(out)
  locals <- .symRedLocals(c(b$invariants, b$support, names(pins), pins), spy)
  Ies <- lapply(b$invariants, function(iv)
    tryCatch(.symRedSympify(gsub("\\^", "**", iv), spy, locals),
             error = function(e) NULL))
  if (any(vapply(Ies, is.null, logical(1)))) return(out)
  # the v-degrees present in a polynomial (NULL when Poly refuses)
  degsIn <- function(p, v) {
    pp <- tryCatch(spy$Poly(p, spy$Symbol(v)), error = function(e) NULL)
    if (is.null(pp)) return(NULL)
    dg <- tryCatch(unlist(.symRedIter(pp$monoms(), function(m) as.integer(m[[1]]))),
                   error = function(e) NULL)
    if (is.null(dg)) NULL else sort(unique(dg))
  }
  # One carrier per invariant that no other invariant has claimed. Candidates are
  # tried from the END of the coordinate list first: znames orders states before
  # parameters, and states are better left to the gauge pins, which are constants.
  used <- character(0); carriers <- character(0)
  for (Ie in Ies) {
    cand <- setdiff(b$support, used)
    if (length(coords)) cand <- cand[order(-match(cand, coords, nomatch = 0L))]
    expA <- tryCatch(.symRedIter(Ie$atoms(spy$exp), function(a) a),
                     error = function(e) list())
    if (length(expA) > 1L) return(out)
    G <- if (length(expA)) {
      a <- expA[[1]]$args                # tuple may auto-convert to an R list
      if (is.list(a)) a[[1]] else reticulate::py_get_item(a, 0L)
    } else NULL
    R <- if (length(expA)) spy$cancel(Ie / expA[[1]]) else Ie
    fracR <- spy$fraction(spy$cancel(spy$together(R)))
    gSyms <- if (is.null(G)) character(0) else .symRedFreeSyms(G)
    rSyms <- .symRedFreeSyms(R)
    # pure-power test: the v-degrees of numerator and denominator inside {0, k}
    # for a single k make N - c*D = alpha*v^k + beta, solvable by a k-th root
    # (k = 1 is the plain linear case)
    powerOK <- function(frac, v) {
      dn <- degsIn(frac[[1]], v); dd <- degsIn(frac[[2]], v)
      if (is.null(dn) || is.null(dd)) return(FALSE)
      length(setdiff(unique(c(dn, dd)), 0L)) == 1L
    }
    fracG <- if (!is.null(G)) spy$fraction(spy$cancel(spy$together(G)))
    pick <- NA_character_
    for (v in cand) {
      if (v %in% gSyms) {
        # log solve: v only inside the exponential; G = log(c) is then a plain
        # rational equation, so the same pure-power test applies to G
        if (v %in% rSyms) next
        if (!powerOK(fracG, v)) next
        pick <- v; break
      }
      if (powerOK(fracR, v)) { pick <- v; break }
    }
    if (is.na(pick)) return(out)
    used <- c(used, pick); carriers <- c(carriers, pick)
  }
  # joint solve of I_l = tmp_l for the carriers; tmp_l renamed to the carrier name
  # afterwards (the equation cannot hold the same symbol in both meanings)
  tmpN <- paste0("dModRedC", seq_along(Ies))
  eqs <- lapply(seq_along(Ies), function(l)
    Ies[[l]] - spy$Symbol(tmpN[l]))

  # ---- positive gauge sections -----------------------------------------------
  # A constant gauge pin is only a LOCAL chart on a curved orbit (the orbit may
  # not reach the pinned value, and the solved entries can leave the positive
  # orthant). Before falling back to pins, search for a monomial-balance section
  # m1 = m2 whose joint solve with the invariants gives entries certified
  # positive on the whole positive orthant -- a global, log-fittable chart.
  gauge0 <- setdiff(b$support, carriers)
  if (length(gauge0) >= 1L && length(gauge0) <= 2L &&
      !any(grepl("exp(", b$invariants, fixed = TRUE))) {
    scalPinPairs <- lapply(names(pins), function(nm)
      reticulate::tuple(spy$Symbol(nm), .symRedSympify(pins[[nm]], spy, locals)))
    allUnk <- c(carriers, gauge0)
    # balances may involve unmoved coordinates appearing in the invariants (on
    # each orbit they are constants, so the section still cuts transversally);
    # a pair touching no moved coordinate cannot constrain the orbit and is
    # dropped
    invSyms <- unique(unlist(lapply(Ies, .symRedFreeSyms)))
    secVars <- sort(unique(c(b$support,
      if (length(coords)) intersect(invSyms, coords) else invSyms)))
    cands <- .symRedSectionCands(secVars)
    touches <- function(pr) any(vapply(b$support, function(v)
      grepl(paste0("(?<![0-9A-Za-z_.])", v, "(?![0-9A-Za-z_.])"),
            paste(pr, collapse = " "), perl = TRUE), logical(1)))
    cands <- Filter(touches, cands)
    # restrict the search to certified-monotone dials: only balances whose ratio
    # moves one way along every orbit reach the solve at all
    if (!is.null(b$preps) && length(cands)) {
      secLoc <- .symRedLocals(c(unlist(lapply(b$preps, `[[`, "comps")),
                                unlist(cands)), spy)
      xiOf <- lapply(b$preps, function(g)
        lapply(g$comps, function(x) .symRedSympify(x, spy, secLoc)))
      need <- if (length(gauge0) == 1L) 40L else 10L
      keep <- list()
      for (pr in cands[seq_len(min(length(cands), 200L))]) {
        if (.symRedSectionMonotone(pr, b$preps, xiOf, spy, secLoc))
          keep[[length(keep) + 1L]] <- pr
        if (length(keep) >= need) break
      }
      cands <- keep
    }
    sets <- if (length(gauge0) == 1L) lapply(cands, list)
      else if (length(cands) > 1L) {
        cmb <- utils::combn(min(length(cands), 10L), 2L)
        lapply(seq_len(ncol(cmb)), function(j)
          list(cands[[cmb[1, j]]], cands[[cmb[2, j]]]))
      } else list()
    getE <- function(solDict, v)
      if (is.list(solDict)) solDict[[v]] else
        reticulate::py_get_item(solDict, spy$Symbol(v), silent = TRUE)
    for (st in sets[seq_len(min(length(sets), 40L))]) {
      secEqs <- lapply(st, function(pr)
        .symRedSympify(pr[1], spy, locals) - .symRedSympify(pr[2], spy, locals))
      sol2 <- tryCatch(spy$solve(c(eqs, secEqs),
                                 lapply(allUnk, function(v) spy$Symbol(v)),
                                 dict = TRUE), error = function(e) NULL)
      if (is.null(sol2) || !length(sol2)) next
      for (bi in seq_along(sol2)) {
        br <- sol2[[bi]]
        ent <- character(0); okBr <- TRUE
        for (v in allUnk) {
          e0 <- getE(br, v)
          if (is.null(e0) || any(allUnk %in% .symRedFreeSyms(e0))) {
            okBr <- FALSE; break
          }
          e <- spy$cancel(e0$subs(scalPinPairs))
          if (!identical(.symRedEntryClass(e, spy), "yes") ||
              !.symRedPosCert(e, spy)) { okBr <- FALSE; break }
          ent[v] <- gsub("\\*\\*", "^", as.character(e))
        }
        if (!okBr) next
        for (l in seq_along(carriers))
          ent <- setNames(gsub(paste0("\\b", tmpN[l], "\\b"), carriers[l], ent),
                          names(ent))
        out$pins <- ent
        out$gauge <- gauge0
        out$section <- vapply(st, function(pr) paste(pr[1], "=", pr[2]),
                              character(1))
        out$meaning <- setNames(b$invariants, carriers)
        out$solved <- TRUE
        return(out)
      }
    }
  }

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
  pinPairs <- lapply(names(allPins), function(nm)
    reticulate::tuple(spy$Symbol(nm), .symRedSympify(allPins[[nm]], spy, locals)))
  es <- vector("list", length(carriers))
  for (l in seq_along(carriers)) {
    # convert=TRUE turns the solution dict into an R list keyed by str(Symbol)
    e <- if (is.list(sol)) sol[[carriers[l]]]
         else reticulate::py_get_item(sol, spy$Symbol(carriers[l]), silent = TRUE)
    if (is.null(e)) return(out)
    es[[l]] <- e
  }
  # chained solutions: an entry referencing another carrier means that carrier's
  # INNER value -- substitute its solved expression until every entry references
  # tmps, pins and outer symbols only
  for (round in seq_along(carriers)) {
    dirty <- FALSE
    for (l in seq_along(carriers)) {
      hit <- setdiff(which(carriers %in% .symRedFreeSyms(es[[l]])), l)
      if (length(hit)) {
        dirty <- TRUE
        es[[l]] <- es[[l]]$subs(lapply(hit, function(l2)
          reticulate::tuple(spy$Symbol(carriers[l2]), es[[l2]])))
      }
    }
    if (!dirty) break
  }
  entries <- character(0); root <- FALSE
  for (l in seq_along(carriers)) {
    if (carriers[l] %in% .symRedFreeSyms(es[[l]])) return(out)   # unresolved chain
    e <- spy$cancel(es[[l]]$subs(pinPairs))
    cls <- .symRedEntryClass(e, spy)
    if (identical(cls, "no")) return(out)
    if (identical(cls, "root")) root <- TRUE
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
  out$rootNote <- if (root)
    "a fractional power in a solved entry assumes the positive branch" else NULL
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
          (!is.null(b$admissible) || !is.null(b$matroid) ||
           !is.null(b$gaugeNote)))
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
    hd <- sprintf("%s block {%s} -- %s%s",
                  if (b$type == "scaling") "Scaling" else "Curved",
                  paste(b$labels, collapse = ", "), b$status,
                  if (!is.null(b$stage) &&
                      !b$stage %in% c("transversal", "fixed", "none"))
                    sprintf("  [stage: %s]", b$stage) else "")
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
    if (!is.null(b$transversal) && length(b$transversal))
      cat(ind, "transversal: ",
          paste(paste0(b$transversal, " = ",
                       if (!is.null(b$pins)) b$pins[b$transversal] else "1"),
                collapse = ", "), "\n", sep = "")
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
    } else cat("\n")
    if (!is.null(f$gaugeNote))
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
#' polynomial up to `dPoly`, rational via Darboux polynomials up to `dDarboux`,
#' Liouvillian via exponential factors `exp(g/h)` with numerator degree up to
#' `dExp` -- and every failed stage leaves a precise negative certificate. By
#' Prelle--Singer the last stage completes the search class: every elementary or
#' Liouvillian first integral is of Darboux form once exponential factors are
#' admitted, so within the degree caps nothing expressible in closed form is
#' missed.
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
#' @param dExp Numerator degree cap of the exponential-factor search (stage 4):
#'   factors `exp(g/h)` with `deg(g) <= dExp` and `h` a product of Darboux
#'   factors. `dExp = 0` skips the stage. Invariants from this stage are
#'   transcendental; the emitted trafo entries may contain `log`/`exp` (fine for
#'   [P()], but not accepted by `symmetryDetection(trafo = )`, which requires
#'   rational entries).
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
                           dExp = 2L, verbose = FALSE, ...) {
  if (!inherits(object, "symmetrydetection"))
    stop("reduceSymmetry(): `object` must be a symmetrydetection result.",
         call. = FALSE)
  if (length(list(...)))
    warning("reduceSymmetry(): unused argument(s) ignored.", call. = FALSE)
  .require_ns("reticulate", "reduceSymmetry()")
  .symCall <- match.call()
  settings <- list(dPoly = as.integer(dPoly), dDarboux = as.integer(dDarboux),
                   dExp = as.integer(dExp),
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
                                      dPoly, dDarboux, dExp, sd, spy, verbose))

  blocks <- .symRedVerify(blocks, sd)

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
      blocks[[bi]]$section <- sol$section
      note <- c(
        if (!is.null(sol$section)) paste0(
          "gauge section ", paste(sol$section, collapse = ", "),
          " -- entries positive for all positive outer values")
        else if (length(sol$gauge)) paste0(
          "gauge pin ", paste(sol$gauge, collapse = ", "),
          " = 1 -- chart valid only where the entries stay positive"),
        sol$rootNote)
      blocks[[bi]]$gaugeNote <- if (length(note)) paste(note, collapse = "; ")
        else NULL
      blocks[[bi]]$survivorMeaning <- sol$meaning
      blocks[[bi]]$status <- "reduced"
      blocks[[bi]]$certificates <- c(b$certificates,
        "solved exactly: each invariant carried under its solved coordinate's name",
        if (!is.null(sol$section))
          "section pre-certified: balance ratio strictly monotone along every orbit")
    } else if (isTRUE(verbose))
      message("curved block {", paste(b$labels, collapse = ", "),
              "}: invariants not rationally solvable -- reported as invariantOnly")
  }

  trafo <- .symRedAssembleTrafo(blocks, coords)
  blocks <- lapply(blocks, function(b) { b$preps <- NULL; b$Wres <- NULL
    b$invExps <- NULL; b })
  .symRedResult(object, blocks, trafo, coords, fixed, settings, .symCall)
}
