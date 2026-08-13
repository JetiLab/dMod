# Behavioral tests for reduceSymmetry() (constructive symmetry reduction).
#
# Assertions are mathematical, never string comparison: invariants are checked by
# X(I) = 0 (exact via .symExprEqual / numeric tangents), reparametrisations by
# re-running symmetryDetection with the emitted trafo and reading the rank.

redquiet <- function(...) suppressWarnings(reduceSymmetry(...))
symdet2 <- function(...) symmetryDetection(..., verbose = FALSE)

# a synthetic symmetrydetection result carrying explicit generators
.mkdir <- function(gen, type = "polynomial", weights = NULL)
  structure(list(type = type, generator = gen, weights = weights,
                 support = names(gen), explicit = TRUE),
            class = "symmetrygenerator")
.mkobj <- function(dirs, coords)
  structure(list(method = "observability", identifiable = FALSE,
                 rank = 1L, dim = length(coords), symmetries = dirs,
                 info = list(engine = "modular", settings = list(),
                             coordinates = coords),
                 call = NULL), class = "symmetrydetection")

# X(I) at a random point for a generator given as named component strings
.lieAt <- function(gen, inv, pt, eps = 1e-7) {
  f <- function(z) eval(parse(text = inv), as.list(z))
  sum(vapply(names(gen), function(v) {
    zp <- pt; zp[v] <- zp[v] + eps
    eval(parse(text = gen[[v]]), as.list(pt)) * (f(zp) - f(pt)) / eps
  }, numeric(1)))
}


test_that("one scaling: transversal pin, family, end-to-end identifiable", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  g <- eqnvec(y = "p")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  b <- red$blocks[[1]]
  expect_identical(b$type, "scaling")
  expect_identical(b$status, "reduced")
  # the invariant lattice is spanned by ktx*ktl and ktx/m (weights 1,-1,1)
  expect_true(any(vapply(b$invariants, function(iv)
    .symExprEqual(gsub("\\^", "**", iv), "ktx*ktl") ||
    .symExprEqual(gsub("\\^", "**", iv), "1/(ktx*ktl)"), logical(1))))
  # every listed transversal is one coordinate of the support
  expect_true(all(unlist(red$family[[1]]$admissible) %in% b$support))
  # the emitted trafo removes the direction
  res2 <- symdet2(f, g, method = "observability", trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("two chained scalings with a shared coordinate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  W1 <- list(a = "1", b = "-1")            # a/b trade-off
  W2 <- list(c = "1", a = "-1")            # c/a trade-off (shares a)
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "-b"), "scaling", W1),
                     .mkdir(list(c = "c", a = "-a"), "scaling", W2)),
                c("a", "b", "c", "other"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_length(b$transversal, 2L)
  # {b, c} is not admissible alone? rank({b,c}) = 2 (rows are independent there),
  # but {a} plus nothing cannot cover both rows: every admissible pair has rank 2
  for (T in red$family[[1]]$admissible)
    expect_identical(dMod:::.symRedRankModP(
      rbind(c(1L, -1L, 0L), c(-1L, 0L, 1L))[, match(T, c("a", "b", "c")),
                                            drop = FALSE]), 2L)
  # fixing the shared coordinate removes only one direction
  red2 <- redquiet(obj, fixed = "a")
  b2 <- red2$blocks[[1]]
  expect_identical(b2$removedByFixed, 1L)
  expect_identical(b2$status, "reduced")
})


test_that("fixed: rank verdict, redundancy, unknown names", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  W <- list(x = "1", y = "-1", z = "-1")
  obj <- .mkobj(list(.mkdir(list(x = "x", y = "-y", z = "-z"), "scaling", W)),
                c("x", "y", "z"))
  expect_warning(reduceSymmetry(obj, fixed = "nothere"), "no effect")
  red <- redquiet(obj, fixed = c("x", "y"))
  b <- red$blocks[[1]]
  expect_identical(b$status, "fixed")
  expect_identical(b$removedByFixed, 1L)
  expect_identical(b$redundantFixed, "y")     # x already kills the single row
  # trafo is pure identity: fixed coordinates are the user's, not pinned here
  expect_true(all(red$trafo == names(red$trafo)))
})


test_that("curved polynomial invariant: solved and identifiable end-to-end", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "k1 + u*k2 - kdeg*A")
  g <- eqnvec(y = "A")
  res <- symdet2(f, g, method = "observability", fixed = "u", reconstruct = TRUE)
  red <- redquiet(res)
  b <- red$blocks[[1]]
  expect_identical(b$status, "reduced")
  expect_identical(b$stage, "polynomial")
  expect_true(any(vapply(b$invariants, function(iv)
    .symExprEqual(gsub("\\^", "**", iv), "k1 + u*k2"), logical(1))))
  # the monomial stage certified negative before escalating
  expect_true(any(grepl("no Laurent-monomial invariant", b$certificates)))
  res2 <- symdet2(f, g, method = "observability", fixed = "u", trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("overlapping scaling and curved directions merge into one block", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g <- eqnvec(y = "alpha*A")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  # the scaling (A, B, alpha) shares B with the affine direction: one joint block
  expect_length(red$blocks, 1L)
  b <- red$blocks[[1]]
  expect_identical(b$type, "curved")
  expect_identical(b$status, "reduced")
  expect_length(b$invariants, 3L)             # dim 5 - rank 2
  res2 <- symdet2(f, g, method = "observability", trafo = red$trafo)
  expect_true(res2$identifiable)
})


test_that("rational invariant via Darboux; degree cap gives a certificate", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  obj <- .mkobj(list(.mkdir(list(x = "x**2", y = "y**2"))), c("x", "y"))
  # cap 0: coordinate factors only, honest negative certificate
  red0 <- redquiet(obj, dDarboux = 0L)
  b0 <- red0$blocks[[1]]
  expect_identical(b0$status, "unresolved")
  expect_true(any(grepl("coordinate and xi factors only", b0$certificates)))
  # cap 1: the extactic factors x, y, x - y give (x - y)/(x*y)
  red1 <- redquiet(obj, dDarboux = 1L)
  b1 <- red1$blocks[[1]]
  expect_identical(b1$stage, "darboux")
  expect_length(b1$invariants, 1L)
  pt <- c(x = 1.31, y = 0.77)
  expect_lt(abs(.lieAt(list(x = "x^2", y = "y^2"), b1$invariants[1], pt)), 1e-4)
})


test_that("module reduction collapses the M009-shaped pair", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  X12 <- list(C = "Ci*rr + Ci", Ci = "-(Ci*rr - Ci)",
              k2i = "-(k2*rr - k2 + k2i*rr + k2i)",
              r1 = "r1*rr + r1", r2 = "r2*rr + r2 + rr")
  X13 <- list(C = "Ci*r2*rr + Ci*r2", Ci = "-(Ci*r2*rr - Ci*r2)",
              k2i = "-(k2*r2*rr - k2*r2 + k2i*r2*rr + k2i*r2)",
              r1 = "-r1*rr",
              rr = "-(r2*rr**3 - 2*r2*rr**2 - r2*rr - rr**3 - rr**2)")
  obj <- .mkobj(list(.mkdir(X12), .mkdir(X13)),
                c("C", "Ci", "k2", "k2i", "r1", "r2", "rr"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  # the discovered combination (-r2)*X1 + X2 lives on 3 coordinates
  expect_true(any(grepl("support 3", b$moduleCombos)))
})


test_that("multi-generator block: only common invariants survive", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  # X1 scales (a, b); X2 scales (b, c): the common monomial lattice is spanned by
  # a*b*c (exponents alpha = beta = gamma) -- a*b is X1-invariant but not
  # X2-invariant, so the joint block must NOT report it
  obj <- .mkobj(list(.mkdir(list(a = "a", b = "-b"), "general"),
                     .mkdir(list(b = "b", c = "-c"), "general")),
                c("a", "b", "c"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_length(b$invariants, 1L)
  expect_true(.symExprEqual(gsub("\\^", "**", b$invariants[1]), "a*b*c") ||
              .symExprEqual(gsub("\\^", "**", b$invariants[1]), "1/(a*b*c)"))
})


test_that("support-only directions and stripped coordinates degrade gracefully", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  d <- structure(list(type = "general", generator = NULL, weights = NULL,
                      support = c("p", "q"), explicit = FALSE),
                 class = "symmetrygenerator")
  obj <- .mkobj(list(d), c("p", "q", "r"))
  red <- redquiet(obj)
  b <- red$blocks[[1]]
  expect_identical(b$status, "unresolved")
  expect_match(b$reason, "reconstruct = TRUE")
  # no coordinate list: support-union fallback with a warning
  obj$info$coordinates <- NULL
  expect_warning(red2 <- reduceSymmetry(obj), "no coordinate list")
  expect_setequal(names(red2$trafo), c("p", "q"))
})


test_that("trafo entries parse and carry integer powers only", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g <- eqnvec(y = "alpha*A")
  res <- symdet2(f, g, method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  expect_true(inherits(red$trafo, "eqnvec"))
  expect_setequal(names(red$trafo), red$coordinates)
  for (v in red$trafo) {
    expect_error(parse(text = v), NA)
    expect_false(grepl("\\^\\s*\\(?\\s*[0-9]*\\.[0-9]", v))   # no fractional powers
  }
})


test_that("identifiable results return an empty reduction", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(x = "-k*x")
  res <- symdet2(f, eqnvec(y = "x"), method = "observability")
  red <- redquiet(res)
  expect_length(red$blocks, 0L)
  expect_null(red$trafo)
  expect_output(print(red), "Nothing to reduce")
})


test_that("print shows blocks, trafo and family", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  res <- symdet2(f, eqnvec(y = "p"), method = "observability", reconstruct = TRUE)
  red <- redquiet(res)
  out <- capture.output(print(red))
  expect_true(any(grepl("Scaling block", out)))
  expect_true(any(grepl("transversal:", out)))
  expect_true(any(grepl("admissible transversal", out)))
  expect_true(any(grepl("Trafo", out)))
})


test_that("EGF/MEK/ERK cascade: all four directions reduced end-to-end", {
  if (!.sympy_works()) skip("reticulate/sympy not available")
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR")   |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR")   |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK")      |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK")              |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK")          |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                        pERK_obs = "scale_pERK * pERK")
  egf <- symdet2(reactions, observables, method = "observability",
                 reduceCQ = FALSE, reconstruct = TRUE)
  red <- redquiet(egf)
  expect_length(red$remaining, 0L)
  # two clean unit scalings plus the receptor scaling merged with the rational
  # confounder into one curved block
  types <- vapply(red$blocks, `[[`, character(1), "type")
  expect_identical(sort(types), c("curved", "scaling", "scaling"))
  cb <- red$blocks[[which(types == "curved")]]
  expect_length(cb$labels, 2L)
  expect_length(cb$invariants, 4L)
  # the confounder invariant totalEGF*totalEGFR/EGF_EGFR^2 is among them
  expect_true(any(vapply(cb$invariants, function(iv) .symExprEqual(
    gsub("\\^", "**", iv),
    "(EGF + EGF_EGFR)*(EGFR + EGF_EGFR)/EGF_EGFR**2"), logical(1))))
  egf2 <- symdet2(reactions, observables, method = "observability",
                  reduceCQ = FALSE, trafo = red$trafo)
  expect_true(egf2$identifiable)
})
