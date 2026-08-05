## ============================================================================
## dMod NLME vs nlmixr2: Warfarin PK/PD cross-tool validation + timing
## ============================================================================
##
## Fits the classic Warfarin PK/PD dataset (32 subjects, single oral dose,
## `nlmixr2data::warfarin`) with dMod's FOCEI and with nlmixr2's FOCEI on
## *identical* data and an identical structural model, then answers two
## questions:
##
##   1. CONVERGENCE. Do dMod and nlmixr2 land on the same population estimates,
##      and are those estimates in the published warfarin ballpark? A single
##      reference fit from each tool is tabulated against nominal literature
##      values (Section 4).
##
##   2. THROUGHPUT / ROBUSTNESS. Run N multistart fits (perturbed structural
##      starts) in each tool, in parallel over `cores`, and compare total
##      wall-clock, per-fit time, and the fraction of starts that reach the
##      global optimum (Sections 5-6). dMod multistart is msEM(); nlmixr2
##      has no native multistart, so N nlmixr2() calls are forked via
##      parallel::mclapply with rxode2 threads pinned to 1 per worker.
##
## Structural model (both tools) -- coupled 1-cmt oral PK + indirect-response
## PD with inhibition of PCA synthesis (Imax fixed to 1, baseline pca0 =
## kin / kout):
##   dGut/dt    = -Ka * Gut
##   dCenter/dt =  Ka * Gut - (Cl/V) * Center            cp  = Center / V
##   dPCA/dt    =  kin * IC50/(IC50 + cp) - kout * PCA    PCA(0) = pca0
##   kin = pca0 * kout
##   V   = exp(tv    + eta_V),    Cl   = exp(tcl   + eta_Cl)
##   pca0 = exp(tpca0 + eta_pca0), kout = exp(tkout + eta_kout)
##   Ka  = exp(tka),  IC50 = exp(tic50)          (population only, N = 32)
##
## Residuals (matched exactly across tools):
##   cp  : log(DV) = log(cp) + eps,  eps ~ N(0, sigma_cp^2)   (dMod additive-on-
##         log <-> nlmixr2 `cp ~ lnorm(sigma_cp)`), so sigma_cp is one number.
##   pca : DV = PCA + eps,  eps ~ N(0, sigma_pca^2)           (both additive).
##
## OFV NOTE: the two tools' objective values are NOT on a common scale and are
## never differenced. dMod models log(cp) with an additive residual (its -2LL
## lives on the log scale); nlmixr2's lnorm() lives on cp and differs by the
## change-of-variables Jacobian 2*sum(log cp_obs). The pca endpoint is additive
## in both, so it contributes no Jacobian. Cross-tool validation is on the
## *estimates*; OFVs are compared only within a tool (the multistart waterfall).
##
## Requires: nlmixr2 (+ rxode2), nlmixr2data. Intended use: source()
## interactively, or `NFITS=10 CORES=10 Rscript inst/examples/compareNlmixr2Warfarin.R`.
## Config via environment variables (defaults in parentheses):
##   NFITS  (50)  multistart fits per tool
##   CORES  (10)  parallel workers
##   PERTSD (0.4) SD of the rnorm perturbation added to the structural starts
## Nothing is written outside tempdir().
## ============================================================================

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(dMod)
}
if (!requireNamespace("nlmixr2", quietly = TRUE))
  stop("This example needs the 'nlmixr2' package. install.packages('nlmixr2').")
if (!requireNamespace("nlmixr2data", quietly = TRUE))
  stop("This example needs 'nlmixr2data' for the warfarin dataset.")

N_FITS <- as.integer(Sys.getenv("NFITS",  "50"))
CORES  <- as.integer(Sys.getenv("CORES",  "10"))
PERT_SD <- as.numeric(Sys.getenv("PERTSD", "0.4"))
cat(sprintf("Config: N_FITS=%d  CORES=%d  PERT_SD=%.2f\n",
            N_FITS, CORES, PERT_SD))

oldwd <- setwd(tempdir())
on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)
set.seed(1)

## ----------------------------------------------------------------------------
## 1. Data (shared source: nlmixr2data::warfarin)
##
## Drop dosing rows for the observation frames, and the 4 cp==0 observations
## (log undefined). pca is kept unchanged. Both tools see exactly these rows.
## ----------------------------------------------------------------------------
data(warfarin, package = "nlmixr2data")
warfarin$id <- as.character(warfarin$id)
obs <- warfarin[warfarin$evid == 0 &
                  !(warfarin$dvid == "cp" & warfarin$dv == 0), , drop = FALSE]
subjects <- sort(unique(obs$id))
N <- length(subjects)
doses <- vapply(subjects, function(s) {
  warfarin[warfarin$id == s & warfarin$evid == 1, , drop = FALSE][1, ]$amt
}, 0.0)
cat(sprintf("Warfarin: N=%d subjects, %d obs (%d cp + %d pca)\n",
            N, nrow(obs), sum(obs$dvid == "cp"), sum(obs$dvid == "pca")))

## ============================================================================
## 2. dMod model (identical setup to inst/examples/warfarin.R)
## ============================================================================
data_df <- rbind(
  data.frame(name = "y_cp",  time = obs$time[obs$dvid == "cp"],
             value = log(obs$dv[obs$dvid == "cp"]), sigma = NA_real_,
             condition = obs$id[obs$dvid == "cp"], stringsAsFactors = FALSE),
  data.frame(name = "y_pca", time = obs$time[obs$dvid == "pca"],
             value = obs$dv[obs$dvid == "pca"], sigma = NA_real_,
             condition = obs$id[obs$dvid == "pca"], stringsAsFactors = FALSE))
dlist <- as.datalist(data_df)

reactions <- eqnlist()
reactions <- addReaction(reactions, "Gut",    "",       "Ka * Gut",  "absorption")
reactions <- addReaction(reactions, "",       "Center", "Ka * Gut",  "appearance")
reactions <- addReaction(reactions, "Center", "",       "(Cl/V) * Center",
                         "elimination")
reactions <- addReaction(reactions, "",       "PCA",
                         "kin * IC50 / (IC50 + Center/V)", "PCA_synthesis")
reactions <- addReaction(reactions, "PCA",    "",       "kout * PCA",
                         "PCA_degradation")
m <- odemodel(reactions, modelname = "warf_ode", compile = FALSE,
              backend = "cppDE", deriv2 = TRUE)
x <- Xs(m, optionsSens = list(atol = 1e-4, rtol = 1e-4))
g <- Y(c(y_cp  = "log(Center/V + 1e-3)", y_pca = "PCA"),
       x, modelname = "warf_obs", compile = FALSE, deriv2 = TRUE)
e <- Y(eqnvec(y_cp = "sigma_cp", y_pca = "sigma_pca"),
       g, attach.input = FALSE, deriv2 = TRUE, compile = FALSE,
       modelname = "warf_err")

trafo <- eqnvec(
  Ka        = "exp(tka)",
  V         = "exp(tv  + eta_V)",
  Cl        = "exp(tcl + eta_Cl)",
  kout      = "exp(tkout + eta_kout)",
  IC50      = "exp(tic50)",
  kin       = "exp(tpca0 + eta_pca0) * exp(tkout + eta_kout)",
  Gut       = "Gut_init",
  Center    = "0",
  PCA       = "exp(tpca0 + eta_pca0)",
  sigma_cp  = "exp(log_sigma_cp)",
  sigma_pca = "exp(log_sigma_pca)")
subj_table <- data.frame(
  eta_V    = paste0("eta_V_",    subjects),
  eta_Cl   = paste0("eta_Cl_",   subjects),
  eta_pca0 = paste0("eta_pca0_", subjects),
  eta_kout = paste0("eta_kout_", subjects),
  Gut_init = doses, row.names = subjects, stringsAsFactors = FALSE)
p <- P(branch(trafo, table = subj_table, apply = "insert"),
       method = "explicit", compile = FALSE, modelname = "warf_p",
       deriv2 = TRUE)
prd <- g * x * p
compile(prd, e, cores = min(CORES, 8L))

om  <- omega(eta = c("eta_V", "eta_Cl", "eta_pca0", "eta_kout"),
             subjects = subjects)
obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)
init <- emInit(c(tka = 0.0, tv = 2.0, tcl = -2.0, tpca0 = 4.6,
                   tkout = -3.0, tic50 = 0.7,
                   log_sigma_cp = log(0.2), log_sigma_pca = log(7.4)),
                 om, sd = 0.3)

focei_ctrl <- list(focei = list(
  innerControl = list(iterlim = 50, fterm = 1e-7, mterm = 1e-7),
  trustControl = list(iterlim = 200)))

## ============================================================================
## 3. nlmixr2 model + event data (same rows, dvid = cp/pca -> endpoints)
## ============================================================================
d_nlmixr <- warfarin[warfarin$evid == 1 |
                       (warfarin$evid == 0 &
                          !(warfarin$dvid == "cp" & warfarin$dv == 0)),
                     c("id", "time", "amt", "dv", "dvid", "evid")]

warf_model <- function() {
  ini({
    tka   <-  0.0
    tv    <-  2.0
    tcl   <- -2.0
    tpca0 <-  4.6
    tkout <- -3.0
    tic50 <-  0.7
    eta.v    ~ 0.09
    eta.cl   ~ 0.09
    eta.pca0 ~ 0.09
    eta.kout ~ 0.09
    cp.sd  <- 0.2
    pca.sd <- 7.4
  })
  model({
    ka   <- exp(tka)
    v    <- exp(tv + eta.v)
    cl   <- exp(tcl + eta.cl)
    kout <- exp(tkout + eta.kout)
    ic50 <- exp(tic50)
    pca0 <- exp(tpca0 + eta.pca0)
    kin  <- pca0 * kout
    d/dt(gut)    = -ka * gut
    d/dt(center) =  ka * gut - (cl / v) * center
    cp = center / v
    d/dt(pca) = kin * ic50 / (ic50 + cp) - kout * pca
    pca(0) = pca0
    cp  ~ lnorm(cp.sd)
    pca ~ add(pca.sd)
  })
}
theta0 <- c(tka = 0.0, tv = 2.0, tcl = -2.0, tpca0 = 4.6,
            tkout = -3.0, tic50 = 0.7)

## One perturbed nlmixr2 FOCEI fit. Runs in an mclapply worker, so pin rxode2
## to a single thread to avoid oversubscription when many fork at once.
run_one_nlmixr <- function(theta) {
  rxode2::setRxThreads(1L)
  ui <- suppressMessages(
    do.call(rxode2::ini,
            c(list(nlmixr2est::nlmixr2(warf_model)), as.list(theta))))
  t <- system.time(
    f <- try(suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(ui, d_nlmixr, est = "focei",
                          control = nlmixr2est::foceiControl(print = 0L)))),
      silent = TRUE))[["elapsed"]]
  if (inherits(f, "try-error"))
    return(list(ok = FALSE, elapsed = t, objf = NA_real_,
                fixef = NULL, omega = NULL))
  list(ok = TRUE, elapsed = t, objf = as.numeric(f$objf),
       fixef = nlmixr2est::fixef(f), omega = f$omega)
}

## ============================================================================
## 4. Single reference fit from each tool -> convergence table
## ============================================================================
cat("\n== dMod EM(method='focei') reference ==\n")
t_dmod1 <- system.time(
  fit_dmod <- EM(obj, init, method = "focei", control = focei_ctrl,
                      verbose = FALSE))[["elapsed"]]

cat("== nlmixr2 FOCEI reference ==\n")
ref_nlmixr <- run_one_nlmixr(theta0)

## dMod estimates on the natural (pharma) scale.
a  <- fit_dmod$argument
Om <- fit_dmod$Omega
dmod_est <- c(
  Ka_pop = exp(a[["tka"]]), V_pop = exp(a[["tv"]]), Cl_pop = exp(a[["tcl"]]),
  pca0_pop = exp(a[["tpca0"]]), kout_pop = exp(a[["tkout"]]),
  IC50_pop = exp(a[["tic50"]]),
  sd_eta_V    = sqrt(Om["eta_V",    "eta_V"]),
  sd_eta_Cl   = sqrt(Om["eta_Cl",   "eta_Cl"]),
  sd_eta_pca0 = sqrt(Om["eta_pca0", "eta_pca0"]),
  sd_eta_kout = sqrt(Om["eta_kout", "eta_kout"]),
  sigma_cp = exp(a[["log_sigma_cp"]]), sigma_pca = exp(a[["log_sigma_pca"]]))

nx_est <- function(r) {
  th <- r$fixef; om <- r$omega
  c(Ka_pop = exp(th[["tka"]]), V_pop = exp(th[["tv"]]), Cl_pop = exp(th[["tcl"]]),
    pca0_pop = exp(th[["tpca0"]]), kout_pop = exp(th[["tkout"]]),
    IC50_pop = exp(th[["tic50"]]),
    sd_eta_V    = sqrt(om["eta.v",    "eta.v"]),
    sd_eta_Cl   = sqrt(om["eta.cl",   "eta.cl"]),
    sd_eta_pca0 = sqrt(om["eta.pca0", "eta.pca0"]),
    sd_eta_kout = sqrt(om["eta.kout", "eta.kout"]),
    sigma_cp = th[["cp.sd"]], sigma_pca = th[["pca.sd"]])
}

## Nominal literature ballpark for this warfarin PK/PD parameterisation
## (O'Reilly warfarin; e.g. the Monolix / nlmixr2 warfarin PK-PD demos). These
## are round reference magnitudes, not a specific published fit.
ref <- c(Ka_pop = 0.5, V_pop = 8.0, Cl_pop = 0.13, pca0_pop = 96,
         kout_pop = 0.05, IC50_pop = 1.2,
         sd_eta_V = NA, sd_eta_Cl = NA, sd_eta_pca0 = NA, sd_eta_kout = NA,
         sigma_cp = NA, sigma_pca = NA)

tab <- data.frame(
  Nominal       = ref,
  dMod_FOCEI    = dmod_est[names(ref)],
  nlmixr2_FOCEI = if (isTRUE(ref_nlmixr$ok)) nx_est(ref_nlmixr)[names(ref)]
                  else NA_real_,
  check.names = FALSE)
tab$`rel.diff.%` <- 100 * (tab$dMod_FOCEI - tab$nlmixr2_FOCEI) / tab$nlmixr2_FOCEI

cat("\n================= Population estimates =================\n")
print(round(tab, 4))
cat(sprintf(
  "\ndMod  reference OFV (-2LL, log-cp scale) : %10.2f  (%.1fs)\n",
  fit_dmod$value, t_dmod1))
if (isTRUE(ref_nlmixr$ok))
  cat(sprintf(
    "nlmixr2 reference objf (lnorm-cp scale)  : %10.2f  (%.1fs)\n",
    ref_nlmixr$objf, ref_nlmixr$elapsed))
cat("  (OFVs on different scales; see header -- compare the estimates.)\n")

## ============================================================================
## 5. Multistart timing: N_FITS fits per tool, CORES workers, matched PERT_SD
## ============================================================================
cat(sprintf("\n== Multistart: %d fits x %d cores (perturb SD = %.2f) ==\n",
            N_FITS, CORES, PERT_SD))

cat("-- dMod msEM --\n")
t_dmod_ms <- system.time(
  ms_dmod <- msEM(obj, init, method = "focei", control = focei_ctrl,
                       fits = N_FITS, cores = CORES, sd = PERT_SD))[["elapsed"]]
pf_dmod <- as.parframe(ms_dmod)

cat("-- nlmixr2 x N via mclapply --\n")
inits_nlmixr <- lapply(seq_len(N_FITS), function(i)
  setNames(theta0 + rnorm(length(theta0), sd = PERT_SD), names(theta0)))
t_nlmixr_ms <- system.time(
  ms_nlmixr <- parallel::mclapply(inits_nlmixr, run_one_nlmixr,
                                  mc.cores = CORES,
                                  mc.preschedule = FALSE))[["elapsed"]]

## ============================================================================
## 6. Report: within-tool waterfall + throughput
## ============================================================================
frac_converged <- function(ofv, tol = 1) {
  ofv <- ofv[is.finite(ofv)]
  if (!length(ofv)) return(c(best = NA, n_ok = 0, n_at_best = 0, frac = NA))
  best <- min(ofv)
  c(best = best, n_ok = length(ofv),
    n_at_best = sum(ofv <= best + tol),
    frac = mean(ofv <= best + tol))
}
dmod_ofv   <- pf_dmod$value
nlmixr_ofv <- vapply(ms_nlmixr,
                     function(r) if (isTRUE(r$ok)) r$objf else NA_real_, 0.0)
per_fit_time_nlmixr <- vapply(ms_nlmixr,
                              function(r) r$elapsed, 0.0)

cat("\n================= Multistart summary =================\n")
smry <- rbind(
  dMod    = c(fits = N_FITS,
              total_s = round(t_dmod_ms, 1),
              per_fit_s = round(t_dmod_ms * CORES / N_FITS, 1),
              frac_converged(dmod_ofv)),
  nlmixr2 = c(fits = N_FITS,
              total_s = round(t_nlmixr_ms, 1),
              per_fit_s = round(median(per_fit_time_nlmixr, na.rm = TRUE), 1),
              frac_converged(nlmixr_ofv)))
print(round(smry, 3))
cat("\n  total_s   = wall-clock for all N fits on `cores` workers\n")
cat("  per_fit_s = dMod: total*cores/N; nlmixr2: median single-fit elapsed\n")
cat("  best      = lowest within-tool OFV (NOT comparable across tools)\n")
cat("  frac      = share of fits within 1 OFV unit of that tool's best\n")

cat("\n-- dMod OFV waterfall (sorted) --\n")
print(round(sort(dmod_ofv), 2))
cat("\n-- nlmixr2 OFV waterfall (sorted) --\n")
print(round(sort(nlmixr_ofv), 2))

invisible(list(fit_dmod = fit_dmod, ref_nlmixr = ref_nlmixr,
               estimates = tab, ms_dmod = ms_dmod, ms_nlmixr = ms_nlmixr,
               timing = c(dmod_total = t_dmod_ms, nlmixr_total = t_nlmixr_ms)))
