# Does the chart symmetryReduction() emits cover the whole positive orthant?
#
# The reduction certifies its chart in one direction -- every emitted entry is
# positive for every admissible value of the outer parameters. This script tests the
# other one: sample z in the positive orthant, read off the carriers, push them back
# through the chart, and ask whether the result is a positive point on the SAME
# orbit. A point that is not is a model the reparametrisation cannot express.
#
# A carrier whose invariant takes both signs is declared `real` by the reduction and
# is scored against that domain, not against positivity.
#
# Usage:
#   Rscript bench/symmetry_chart_coverage.R
#   source("bench/symmetry_chart_coverage.R"); res <- chart_coverage_run()
#   chart_coverage_save(res); chart_coverage_compare(res)

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(dMod.sym.verbose = FALSE)

.cc_outdir <- normalizePath(".")


# ---- fixtures --------------------------------------------------------------
# Chosen for the shape of their orbit space, not for size: a scaling block, a
# curved block with a radical chart, a bounded orbit window, and the two blocks
# whose invariant is a difference (one positive by inspection, one indefinite).

chart_cov_fixtures <- function() {
  ab <- eqnlist() |> addReaction("A", "B", "k1*A") |> addReaction("B", "A", "k2*B")
  egfmm <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK / (Km_MEK + MEK)") |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK / (Km_pMEK + pMEK)")       |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK / (Km_ERK + ERK)")     |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK / (Km_pERK + pERK)")

  list(
    list(name = "ab-scaling", f = ab, g = eqnvec(Aobs = "alpha*A"),
         args = list(reduceCQ = FALSE)),
    list(name = "enzyme", f = eqnvec(S = "-kcat*Etot*S/(Km + S)"),
         g = eqnvec(y = "s*S"), args = list()),
    list(name = "gene", f = eqnvec(mR = "ktx - dm*mR", prot = "ktl*mR - dp*prot"),
         g = eqnvec(y = "prot"), args = list()),
    list(name = "hill", f = eqnlist() |>
           addReaction("0",  "FB", "k_pr_FB")                     |>
           addReaction("FB", "0",  "d_FB * FB")                   |>
           addReaction("0",  "x",  "k_pr_x * K^n / (K^n + FB^n)") |>
           addReaction("x",  "0",  "d_x * x"),
         g = eqnvec(xobs = "scale * x"), args = list()),
    # bounded orbit window: the orbit is a quarter arc of x1^2 + x2^2 = const
    list(name = "rotation", f = eqnvec(x1 = "-x2", x2 = "x1"),
         g = eqnvec(y = "x1^2 + x2^2"), args = list(reduceCQ = FALSE)),
    # invariant is a SUM, positive by inspection -> carrier stays log-fittable
    list(name = "moiety-shift", f = eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B"),
         g = eqnvec(y = "A + B"), args = list(reduceCQ = FALSE)),
    # invariant is a DIFFERENCE with no sign: needs a real carrier and the
    # sum-of-squares pin, a constant pin covers only part of the orthant
    list(name = "difference", f = eqnvec(x = "-(a - b)*x"), g = eqnvec(y = "x"),
         args = list()),
    # a curved block whose chart carries radicals
    list(name = "egf-mm-steadystate", f = egfmm,
         g = eqnvec(pMEK_obs = "pMEK", pERK_obs = "pERK"),
         events = eventlist(var = "EGF", time = 0, value = "1", method = "add"),
         trafo = c(
           k_phos_MEK = "k_dephos_MEK*pMEK*(Km_MEK + MEK)/(EGF_EGFR*MEK*(Km_pMEK + pMEK))",
           k_phos_ERK = "k_dephos_ERK*pERK*(ERK + Km_ERK)/(ERK*pMEK*(Km_pERK + pERK))",
           EGF        = "EGF_EGFR*k_unbind/(EGFR*k_bind)"),
         args = list(reduceCQ = FALSE))
  )
}


# ---- one fixture -----------------------------------------------------------

.cc_eval <- function(e, env)
  suppressWarnings(eval(parse(text = gsub("\\*\\*", "^", e)), env))

chart_coverage_fixture <- function(fx, n = 20000, seed = 5) {
  none <- function(status) data.frame(fixture = fx$name, carriers = "", real = 0L,
    coverage = NA_character_, inChart = NA_real_, onOrbit = NA_real_, status = status)

  r <- try(do.call(symmetryDetection, c(list(
    f = fx$f, g = fx$g, method = "observability", reconstruct = TRUE,
    events = fx$events,
    trafo = if (is.null(fx$trafo)) NULL else as.eqnvec(fx$trafo),
    verbose = FALSE), fx$args)), silent = TRUE)
  if (inherits(r, "try-error")) return(none("detection failed"))
  if (!length(r$symmetries)) return(none("identifiable"))
  red <- try(symmetryReduction(r, verbose = FALSE), silent = TRUE)
  if (inherits(red, "try-error")) return(none("reduction failed"))
  meaning <- unlist(lapply(red$blocks, `[[`, "survivorMeaning"))
  if (!length(meaning)) return(none("nothing reduced"))

  co <- red$coordinates
  set.seed(seed)
  Z <- as.data.frame(matrix(exp(runif(n * length(co), -4, 4)), ncol = length(co)))
  names(Z) <- co

  Q <- setNames(lapply(meaning, .cc_eval, env = Z), names(meaning))
  surv <- setdiff(getSymbols(as.character(red$trafo)), names(meaning))
  Zp <- lapply(red$trafo, .cc_eval, env = c(Q, as.list(Z)[intersect(surv, co)]))
  Zp <- lapply(Zp, function(v) if (length(v) == 1L) rep(v, n) else v)

  dom <- unlist(lapply(red$blocks, `[[`, "carrierDomain"))
  isReal <- function(k) identical(unname(dom[k]), "real")
  okQ <- Reduce(`&`, lapply(names(Q), function(k)
    if (isReal(k)) is.finite(Q[[k]]) else is.finite(Q[[k]]) & Q[[k]] > 0))
  okZ <- Reduce(`&`, lapply(Zp, function(v) is.finite(v) & v > 0))
  # the chart point must sit on the orbit its carriers name
  Zpd <- setNames(as.data.frame(Zp), names(red$trafo))
  onOrbit <- Reduce(`&`, lapply(names(meaning), function(k) {
    v <- .cc_eval(meaning[[k]], Zpd)
    is.finite(v) & abs(v - Q[[k]]) <= 1e-6 * (1 + abs(Q[[k]])) }))

  data.frame(fixture = fx$name,
    carriers = paste(names(meaning), collapse = ","),
    real = sum(vapply(names(meaning), isReal, logical(1))),
    coverage = paste(unique(unlist(lapply(red$blocks, `[[`, "coverage"))),
                     collapse = "/"),
    inChart = 100 * mean(okQ & okZ),
    onOrbit = 100 * mean(okQ & okZ & onOrbit),
    status = "")
}


# ---- driver ----------------------------------------------------------------

chart_coverage_run <- function(fixtures = chart_cov_fixtures(), only = NULL,
                               tol = 99.9, verbose = TRUE) {
  if (is.null(only)) only <- strsplit(Sys.getenv("DMOD_CHARTCOV_ONLY", ""), ",")[[1]]
  only <- only[nzchar(only)]
  if (length(only)) fixtures <- Filter(function(fx) fx$name %in% only, fixtures)
  res <- do.call(rbind, lapply(fixtures, chart_coverage_fixture))
  res$status[res$status == ""] <- ifelse(
    is.na(res$inChart[res$status == ""]) |
      res$onOrbit[res$status == ""] >= tol, "ok", "PARTIAL")
  attr(res, "tol") <- tol
  if (verbose) chart_coverage_report(res)
  invisible(res)
}

chart_coverage_report <- function(res) {
  cat("\n", strrep("-", 78), "\n", sep = "")
  cat(sprintf("chart coverage of the positive orthant   pass at >= %.1f%%\n",
              attr(res, "tol")))
  cat(strrep("-", 78), "\n", sep = "")
  pr <- res
  pr$inChart <- ifelse(is.finite(pr$inChart), sprintf("%6.2f%%", pr$inChart), "-")
  pr$onOrbit <- ifelse(is.finite(pr$onOrbit), sprintf("%6.2f%%", pr$onOrbit), "-")
  print(pr, right = FALSE)
  bad <- sum(res$status == "PARTIAL")
  cat("\n", if (bad)
    sprintf("%d chart(s) cover only part of the positive orthant.\n", bad)
    else "every chart covers the whole positive orthant.\n", sep = "")
  invisible(res)
}

chart_coverage_save <- function(res, file = file.path(.cc_outdir, "bench", "baselines",
                                                      "symmetry_chart_coverage.rds")) {
  saveRDS(res, file); cat("baseline written to ", file, "\n", sep = ""); invisible(file)
}

chart_coverage_compare <- function(res, file = file.path(.cc_outdir, "bench", "baselines",
                                                         "symmetry_chart_coverage.rds")) {
  if (!file.exists(file)) stop("no baseline at ", file)
  old <- readRDS(file)
  m <- merge(old[c("fixture", "onOrbit", "status")],
             res[c("fixture", "onOrbit", "status")], by = "fixture",
             suffixes = c("_old", ""))
  m$changed <- m$status_old != m$status
  print(m, right = FALSE)
  invisible(m)
}

if (!interactive() && sys.nframe() == 0L) {
  res <- chart_coverage_run()
  if (any(res$status == "PARTIAL")) quit(status = 1L)
}
