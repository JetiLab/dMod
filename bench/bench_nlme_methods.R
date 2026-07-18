## ============================================================================
## Per-method NLME runtime harness + prediction correctness net.
##
## Two jobs:
##   1. TIMING/PROFILE ("baseline") -- run selected nlmeFit methods, capture
##      wall-clock, OFV, argument, Omega, iters and top Rprof self-time; write a
##      baseline .rds for cross-commit comparison.
##   2. CORRECTNESS NET ("snapshot"/"compare") -- evaluate the composed
##      prediction g*x*p and the error model e at a fixed parameter vector and
##      persist the full prdframe (values + deriv + deriv2 arrays, with
##      dimnames). "compare" reloads and asserts BIT-IDENTICAL. This is the net
##      that guards the Tier-1 glue-caching / Tier-2 Rcpp rewrites: any change
##      claimed to be numerically neutral must leave this snapshot unchanged.
##
## Usage (env-configurable):
##   DMOD_BENCH_MODEL   warfarin | theoph        (default warfarin)
##   DMOD_BENCH_METHODS comma list of methods    (default focei)
##                        focei,quadrature,foceiQuadrature,laplaceEM,saem
##   DMOD_BENCH_MODE    baseline | snapshot | compare   (default baseline)
##   DMOD_BENCH_OUTDIR  scratch dir              (default tempdir()/dMod_bench_nlme)
##
##   Rscript bench/bench_nlme_methods.R
##   DMOD_BENCH_MODE=snapshot Rscript bench/bench_nlme_methods.R
##   DMOD_BENCH_MODE=compare  Rscript bench/bench_nlme_methods.R
##   DMOD_BENCH_METHODS=focei,laplaceEM,saem Rscript bench/bench_nlme_methods.R
## ============================================================================

.dmod_root <- "/home/simon/Documents/Projects/dMod"
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(.dmod_root, quiet = TRUE)
} else {
  library(dMod)
}

model_name   <- Sys.getenv("DMOD_BENCH_MODEL",   "warfarin")
methods_spec <- Sys.getenv("DMOD_BENCH_METHODS", "focei")
mode         <- Sys.getenv("DMOD_BENCH_MODE",    "baseline")
# Stable across processes (snapshot writes, a later compare process reads) --
# NOT under tempdir(), which differs per R session.
scratch_dir  <- Sys.getenv("DMOD_BENCH_OUTDIR", "/tmp/dMod_bench_nlme")
dir.create(scratch_dir, recursive = TRUE, showWarnings = FALSE)
methods <- trimws(strsplit(methods_spec, ",")[[1]])

wd <- file.path(tempdir(), "bench_nlme_wd")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
oldwd <- setwd(wd); on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)
set.seed(1)

## ----------------------------------------------------------------------------
## Model builders. Each returns list(obj, init, prd, e, pars_probe).
## ----------------------------------------------------------------------------
build_theoph <- function() {
  data(Theoph, package = "datasets")
  Theoph$Subject <- as.character(Theoph$Subject)
  Theoph_pos <- Theoph[Theoph$Time > 0, , drop = FALSE]
  subjects <- sort(unique(Theoph_pos$Subject))
  doses <- vapply(subjects, function(s) {
    rec <- Theoph[Theoph$Subject == s, ][1, ]; rec$Dose * rec$Wt }, 0.0)
  dlist <- as.datalist(data.frame(
    name = "y", time = Theoph_pos$Time, value = log(Theoph_pos$conc),
    sigma = NA_real_, condition = Theoph_pos$Subject, stringsAsFactors = FALSE))
  reactions <- eqnlist()
  reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",     "absorption")
  reactions <- addReaction(reactions, "",   "Cc", "Ka * Ag / V", "appearance")
  reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",   "elimination")
  m <- odemodel(reactions, modelname = "bt_ode", compile = FALSE,
                solver = "CppODE", deriv2 = TRUE)
  x <- Xs(m)
  g <- Y(c(y = "log(Cc + 1e-9)"), x, modelname = "bt_obs",
         compile = FALSE, deriv2 = TRUE)
  e <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE, deriv2 = TRUE,
         compile = FALSE, modelname = "bt_err")
  trafo <- eqnvec(Ka = "exp(tka + eta_Ka)", V = "exp(tv + eta_V)",
                  Cl = "exp(tcl + eta_Cl)", Ag = "Ag_init", Cc = "0",
                  sigma_add = "exp(log_sigma_add)")
  subj_table <- data.frame(
    eta_Ka = paste0("eta_Ka_", subjects), eta_V = paste0("eta_V_", subjects),
    eta_Cl = paste0("eta_Cl_", subjects), Ag_init = doses,
    row.names = subjects, stringsAsFactors = FALSE)
  p <- P(branch(trafo, table = subj_table, apply = "insert"),
         method = "explicit", compile = FALSE, modelname = "bt_p", deriv2 = TRUE)
  prd <- g * x * p
  compile(prd, e, cores = 4)
  om  <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
  obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
           constraintL2(mu = 0, Omega = om)
  init <- nlmeInit(c(tka = 0, tv = 3, tcl = 0.7, log_sigma_add = log(0.2)),
                   om, sd = 0.3)
  list(obj = obj, init = init, prd = prd, e = e, om = om,
       pars_probe = init)
}

build_warfarin <- function() {
  if (!requireNamespace("nlmixr2data", quietly = TRUE))
    stop("warfarin model needs nlmixr2data.")
  data(warfarin, package = "nlmixr2data")
  warfarin$id <- as.character(warfarin$id)
  obs <- warfarin[warfarin$evid == 0 &
                    !(warfarin$dvid == "cp" & warfarin$dv == 0), , drop = FALSE]
  subjects <- sort(unique(obs$id))
  doses <- vapply(subjects, function(s)
    warfarin[warfarin$id == s & warfarin$evid == 1, , drop = FALSE][1, ]$amt, 0.0)
  data_df <- rbind(
    data.frame(name = "y_cp", time = obs$time[obs$dvid == "cp"],
               value = log(obs$dv[obs$dvid == "cp"]), sigma = NA_real_,
               condition = obs$id[obs$dvid == "cp"], stringsAsFactors = FALSE),
    data.frame(name = "y_pca", time = obs$time[obs$dvid == "pca"],
               value = obs$dv[obs$dvid == "pca"], sigma = NA_real_,
               condition = obs$id[obs$dvid == "pca"], stringsAsFactors = FALSE))
  dlist <- as.datalist(data_df)
  reactions <- eqnlist()
  reactions <- addReaction(reactions, "Gut", "", "Ka * Gut", "absorption")
  reactions <- addReaction(reactions, "", "Center", "Ka * Gut", "appearance")
  reactions <- addReaction(reactions, "Center", "", "(Cl/V) * Center", "elimination")
  reactions <- addReaction(reactions, "", "PCA", "kin * IC50 / (IC50 + Center/V)",
                           "PCA_synthesis")
  reactions <- addReaction(reactions, "PCA", "", "kout * PCA", "PCA_degradation")
  m <- odemodel(reactions, modelname = "bw_ode", compile = FALSE,
                solver = "CppODE", deriv2 = TRUE)
  x <- Xs(m, optionsSens = list(atol = 1e-4, rtol = 1e-4))
  g <- Y(c(y_cp = "log(Center/V + 1e-3)", y_pca = "PCA"), x,
         modelname = "bw_obs", compile = FALSE, deriv2 = TRUE)
  e <- Y(eqnvec(y_cp = "sigma_cp", y_pca = "sigma_pca"), g,
         attach.input = FALSE, deriv2 = TRUE, compile = FALSE, modelname = "bw_err")
  trafo <- eqnvec(Ka = "exp(tka)", V = "exp(tv + eta_V)", Cl = "exp(tcl + eta_Cl)",
                  kout = "exp(tkout + eta_kout)", IC50 = "exp(tic50)",
                  kin = "exp(tpca0 + eta_pca0) * exp(tkout + eta_kout)",
                  Gut = "Gut_init", Center = "0", PCA = "exp(tpca0 + eta_pca0)",
                  sigma_cp = "exp(log_sigma_cp)", sigma_pca = "exp(log_sigma_pca)")
  subj_table <- data.frame(
    eta_V = paste0("eta_V_", subjects), eta_Cl = paste0("eta_Cl_", subjects),
    eta_pca0 = paste0("eta_pca0_", subjects), eta_kout = paste0("eta_kout_", subjects),
    Gut_init = doses, row.names = subjects, stringsAsFactors = FALSE)
  p <- P(branch(trafo, table = subj_table, apply = "insert"),
         method = "explicit", compile = FALSE, modelname = "bw_p", deriv2 = TRUE)
  prd <- g * x * p
  compile(prd, e, cores = 8)
  om <- omega(eta = c("eta_V", "eta_Cl", "eta_pca0", "eta_kout"), subjects = subjects)
  obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
           constraintL2(mu = 0, Omega = om)
  init <- nlmeInit(c(tka = 0, tv = 2, tcl = -2, tpca0 = 4.6, tkout = -3,
                     tic50 = 0.7, log_sigma_cp = log(0.2), log_sigma_pca = log(7.4)),
                   om, sd = 0.3)
  list(obj = obj, init = init, prd = prd, e = e, om = om, pars_probe = init)
}

cat(sprintf("[bench] model=%s methods=%s mode=%s\n",
            model_name, methods_spec, mode))
M <- switch(model_name, theoph = build_theoph(), warfarin = build_warfarin(),
            stop("unknown DMOD_BENCH_MODEL: ", model_name))

## Control per method (kept modest so runs finish; identical across commits).
ctrl_for <- function(method) switch(method,
  focei = list(focei = list(innerControl = list(iterlim = 50, fterm = 1e-7,
                                                mterm = 1e-7),
                            trustControl = list(iterlim = 200))),
  foceiQuadrature = list(focei = list(trustControl = list(iterlim = 200))),
  laplaceEM = list(quadrature = list(maxEcm = 60L)),
  saem = list(saem = list(nBurnin = 80L, nEM = 80L, nMcmc = 10L)),
  list())

## ----------------------------------------------------------------------------
## Correctness-net snapshot: full prdframe of prd() and e() at pars_probe.
## ----------------------------------------------------------------------------
snap_path <- file.path(scratch_dir, paste0("pred_snapshot_", model_name, ".rds"))

## Extract ONLY the numeric content (values + deriv + deriv2 arrays with
## dimnames) per condition, so the bit-identical compare is not tripped up by
## environments/closures hiding in prdframe attributes.
clean_prdlist <- function(pr) {
  lapply(pr, function(cond) {
    pv <- unclass(cond)
    attr(pv, "deriv") <- attr(pv, "deriv2") <- attr(pv, "parameters") <- NULL
    list(prediction = pv,
         deriv  = attr(cond, "deriv"),
         deriv2 = attr(cond, "deriv2"))
  })
}
capture_pred <- function() {
  times <- switch(model_name, warfarin = seq(0, 120, by = 1),
                  theoph = seq(0, 25, by = 0.5), seq(0, 24, by = 0.5))
  # prd (g*x*p) expects the FULL expanded parameter set: structural + cholPars
  # (in init) plus every per-subject eta. Snapshot at etas = 0.
  all_etas <- as.vector(M$om$subjectEtas)
  pars_full <- c(M$pars_probe,
                 setNames(rep(0, length(all_etas)), all_etas))
  pr <- M$prd(times = times, pars = pars_full, deriv = TRUE, deriv2 = TRUE)
  clean_prdlist(pr)
}

if (mode == "snapshot") {
  snap <- capture_pred()
  saveRDS(snap, snap_path)
  cat("Wrote prediction snapshot:", snap_path, "\n")
  cat("conditions:", length(snap), " | names:",
      paste(head(names(snap), 3), collapse = ","), "...\n")
  cat("cond[1] dims: pred", paste(dim(snap[[1]]$prediction), collapse = "x"),
      "| deriv", paste(dim(snap[[1]]$deriv), collapse = "x"),
      "| deriv2", paste(dim(snap[[1]]$deriv2), collapse = "x"), "\n")
  quit(save = "no")
}

if (mode == "compare") {
  if (!file.exists(snap_path)) stop("No snapshot at ", snap_path,
                                    " -- run mode=snapshot first.")
  ref <- readRDS(snap_path)
  cur <- capture_pred()
  ident <- isTRUE(all.equal(ref, cur, tolerance = 0))
  if (ident) {
    cat("PREDICTION SNAPSHOT: BIT-IDENTICAL ✓\n")
  } else {
    cat("PREDICTION SNAPSHOT: DIFFERS ✗\n")
    print(all.equal(ref, cur, tolerance = 0))
    quit(save = "no", status = 1)
  }
  quit(save = "no")
}

## ----------------------------------------------------------------------------
## Baseline: run + time + profile each requested method.
## ----------------------------------------------------------------------------
run_method <- function(method) {
  rprof <- file.path(scratch_dir, paste0("prof_", model_name, "_", method, ".out"))
  Rprof(rprof, interval = 0.01, line.profiling = FALSE)
  t <- system.time(
    fit <- tryCatch(
      nlmeFit(M$obj, M$init, method = method, control = ctrl_for(method),
              verbose = FALSE),
      error = function(err) err))[["elapsed"]]
  Rprof(NULL)
  if (inherits(fit, "error")) {
    cat(sprintf("  [%s] ERROR: %s\n", method, conditionMessage(fit)))
    return(list(method = method, elapsed = t, error = conditionMessage(fit)))
  }
  ps <- summaryRprof(rprof)
  top <- head(ps$by.self[, c("self.time", "self.pct"), drop = FALSE], 8)
  cat(sprintf("  [%s] wall=%.1fs OFV=%.3f iter=%s\n",
              method, t, fit$value, format(fit$iterations %||% NA)))
  list(method = method, elapsed = t, value = fit$value,
       argument = fit$argument, Omega = fit$Omega,
       iterations = fit$iterations, top_self = top)
}

cat("[bench] running methods...\n")
results <- lapply(methods, run_method)
names(results) <- methods

base_rds <- file.path(scratch_dir, paste0("baseline_", model_name, ".rds"))
record <- list(timestamp = Sys.time(),
               dmod_sha = tryCatch(system(paste("git -C", shQuote(.dmod_root),
                                                "rev-parse --short HEAD"),
                                          intern = TRUE),
                                   error = function(e) NA_character_),
               model = model_name, results = results)
saveRDS(record, base_rds)
cat("\nWrote baseline:", base_rds, " (sha ", record$dmod_sha, ")\n", sep = "")

cat("\n================ Summary ================\n")
for (r in results) {
  if (!is.null(r$error)) { cat(sprintf("%-16s ERROR\n", r$method)); next }
  cat(sprintf("%-16s %7.1fs  OFV=%.3f  iter=%s\n",
              r$method, r$elapsed, r$value, format(r$iterations %||% NA)))
}
cat("\nTop self-time (first method):\n")
if (length(results) && is.null(results[[1]]$error)) print(results[[1]]$top_self)
