## ============================================================================
## dMod NLME vs nlmixr2: Theophylline cross-tool validation
## ============================================================================
##
## Fits the classic Theophylline dataset (12 subjects, 1-compartment oral PK)
## with dMod's FOCEI and with nlmixr2's FOCEI and SAEM, on *identical* data
## and an identical model, and tabulates the population estimates, OFV, and
## wall-clock runtime side by side. This is the runnable validation that
## backs up the "parity with NONMEM / nlmixr2" claim: dMod's niche is exact
## derivatives through the whole ODE + sensitivity + FOCEI-aggregation
## pipeline, and this script shows the estimates land in the same place.
##
## Model (both tools):
##   dAg/dt   = -Ka * Ag                        (gut depot)
##   dCc/dt   =  Ka * Ag / V - (Cl/V) * Cc       (central concentration)
##   Ka = exp(tka + eta.ka), V = exp(tv + eta.v), Cl = exp(tcl + eta.cl)
##   log(DV) = log(Cc/ ... ) + eps,  eps ~ N(0, sigma^2)   (lognormal residual)
##
## The residual is additive on the log scale in dMod (log(DV) = log(Cc) + eps)
## and the exactly-equivalent lognormal error `cp ~ lnorm(add.sd)` in nlmixr2.
##
## OFV: the cross-tool objective values are NOT on a common scale here and are
## printed with an explicit caveat rather than differenced. dMod models log(DV)
## (its -2LL lives on the log scale), whereas nlmixr2's lnorm() likelihood lives
## on DV and differs by the change-of-variables Jacobian 2*sum(log DV); SAEM adds
## its own quadrature constant. The population *estimates* are the cross-tool
## validation. (For a genuine OFV match, run both tools on the *same* likelihood
## scale; dMod's value_nonmem = value - N_obs*log(2*pi) then matches a NONMEM
## OBJ of the same model.)
##
## Reference values (Pinheiro & Bates 2000):
##   Ka_pop ~ 1.50   V_pop ~ 32   Cl_pop ~ 3.0
##   sd(eta_Ka) ~ 0.5   sd(eta_V) ~ 0.1   sd(eta_Cl) ~ 0.3
##
## Requires: nlmixr2 (+ rxode2). saemix is optional and not used here.
## Intended use: source() interactively. Nothing is written to disk.
## ============================================================================

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(quiet = TRUE)
} else {
  library(dMod)
}

have_nlmixr2 <- requireNamespace("nlmixr2", quietly = TRUE)
if (!have_nlmixr2)
  stop("This example needs the 'nlmixr2' package. install.packages('nlmixr2').")

oldwd <- setwd(tempdir())
on.exit(setwd(oldwd), add = TRUE)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE),
       force = TRUE)
set.seed(1)

## ----------------------------------------------------------------------------
## 1. Data (shared source: datasets::Theoph)
## ----------------------------------------------------------------------------
data(Theoph, package = "datasets")
Theoph$Subject <- as.character(Theoph$Subject)
Theoph_pos <- Theoph[Theoph$Time > 0, , drop = FALSE]
subjects   <- sort(unique(Theoph_pos$Subject))
N          <- length(subjects)
doses <- vapply(subjects, function(s) {
  rec <- Theoph[Theoph$Subject == s, ][1, ]
  rec$Dose * rec$Wt
}, 0.0)

## ============================================================================
## 2. dMod FOCEI fit
## ============================================================================
dlist <- as.datalist(data.frame(
  name = "y", time = Theoph_pos$Time, value = log(Theoph_pos$conc),
  sigma = NA_real_, condition = Theoph_pos$Subject, stringsAsFactors = FALSE))

reactions <- eqnlist()
reactions <- addReaction(reactions, "Ag", "",  "Ka * Ag",      "absorption")
reactions <- addReaction(reactions, "",  "Cc", "Ka * Ag / V",  "appearance")
reactions <- addReaction(reactions, "Cc", "",  "Cl/V * Cc",    "elimination")
m <- odemodel(reactions, modelname = "cmp_theoph_ode", compile = FALSE,
              backend = "cppDE", deriv2 = TRUE)
x <- Xs(m)
g <- Y(c(y = "log(Cc + 1e-9)"), x, modelname = "cmp_theoph_obs",
       compile = FALSE, deriv2 = TRUE)
e <- Y(eqnvec(y = "sigma_add"), g, attach.input = FALSE, deriv2 = TRUE,
       compile = FALSE, modelname = "cmp_theoph_err")

trafo <- eqnvec(Ka        = "exp(tka + eta_Ka)",
                V         = "exp(tv  + eta_V)",
                Cl        = "exp(tcl + eta_Cl)",
                Ag        = "Ag_init",
                Cc        = "0",
                sigma_add = "exp(log_sigma_add)")
subj_table <- data.frame(
  eta_Ka = paste0("eta_Ka_", subjects),
  eta_V  = paste0("eta_V_",  subjects),
  eta_Cl = paste0("eta_Cl_", subjects),
  Ag_init = doses, row.names = subjects, stringsAsFactors = FALSE)
p   <- P(branch(trafo, table = subj_table, apply = "insert"),
         method = "explicit", compile = FALSE, modelname = "cmp_theoph_p",
         deriv2 = TRUE)
prd <- g * x * p
compile(prd, e, cores = 4)

om  <- omega(eta = c("eta_Ka", "eta_V", "eta_Cl"), subjects = subjects)
obj <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) +
         constraintL2(mu = 0, Omega = om)
init <- emInit(c(tka = 0.0, tv = 3.0, tcl = 0.7, log_sigma_add = log(0.2)),
                 om, sd = 0.3)

cat("\n== dMod EM(method = 'focei') ==\n")
t_dmod <- system.time(
  fit_dmod <- EM(obj, init, method = "focei",
                      control = list(focei = list(
                        innerControl = list(iterlim = 30),
                        trustControl = list(iterlim = 100))),
                      verbose = FALSE))[["elapsed"]]

## ============================================================================
## 3. nlmixr2 FOCEI + SAEM on the same data
## ============================================================================
# rxode2/nlmixr2 event data built from the SAME datasets::Theoph rows.
dose_rows <- data.frame(ID = subjects, TIME = 0, DV = NA_real_,
                        AMT = doses[subjects], EVID = 1L, CMT = "depot",
                        stringsAsFactors = FALSE)
obs_rows  <- data.frame(ID = Theoph_pos$Subject, TIME = Theoph_pos$Time,
                        DV = Theoph_pos$conc, AMT = 0, EVID = 0L, CMT = "center",
                        stringsAsFactors = FALSE)
nlmixr_data <- rbind(dose_rows, obs_rows)
nlmixr_data <- nlmixr_data[order(match(nlmixr_data$ID, subjects),
                                 nlmixr_data$TIME, -nlmixr_data$EVID), ]

theo_model <- function() {
  ini({
    tka <- 0.0;  tv <- 3.0;  tcl <- 0.7
    eta.ka ~ 0.09; eta.v ~ 0.09; eta.cl ~ 0.09     # var = 0.3^2
    add.sd <- 0.2
  })
  model({
    ka <- exp(tka + eta.ka)
    v  <- exp(tv  + eta.v)
    cl <- exp(tcl + eta.cl)
    d/dt(depot)  = -ka * depot
    d/dt(center) =  ka * depot - (cl / v) * center
    cp = center / v
    cp ~ lnorm(add.sd)     # log(DV) = log(cp) + eps -- matches dMod exactly
  })
}

run_nlmixr <- function(est) {
  ctl <- switch(est,
                focei = nlmixr2est::foceiControl(print = 0L),
                saem  = nlmixr2est::saemControl(print = 0L),
                list())
  t <- system.time(
    f <- suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(theo_model, nlmixr_data, est = est, control = ctl)))
  )[["elapsed"]]
  list(fit = f, elapsed = t)
}

cat("== nlmixr2 FOCEI ==\n"); nx_focei <- run_nlmixr("focei")
cat("== nlmixr2 SAEM ==\n");  nx_saem  <- try(run_nlmixr("saem"), silent = TRUE)

## ============================================================================
## 4. Comparison table
## ============================================================================
# dMod estimates (log-scale structural pars -> natural scale; Omega SD).
a  <- fit_dmod$argument
Om <- fit_dmod$Omega
dmod_est <- c(
  Ka_pop = exp(a[["tka"]]), V_pop = exp(a[["tv"]]), Cl_pop = exp(a[["tcl"]]),
  sd_eta_Ka = sqrt(Om["eta_Ka", "eta_Ka"]),
  sd_eta_V  = sqrt(Om["eta_V",  "eta_V"]),
  sd_eta_Cl = sqrt(Om["eta_Cl", "eta_Cl"]),
  sigma = exp(a[["log_sigma_add"]]))

# nlmixr2 estimates: fixed effects on the same log scale, omega SDs from $omega.
nx_estimates <- function(f) {
  th <- nlmixr2est::fixef(f)       # named: tka, tv, tcl, add.sd (declared scale)
  om <- f$omega                    # K x K, dimnames eta.ka / eta.v / eta.cl
  c(Ka_pop = exp(th[["tka"]]), V_pop = exp(th[["tv"]]), Cl_pop = exp(th[["tcl"]]),
    sd_eta_Ka = sqrt(om["eta.ka", "eta.ka"]),
    sd_eta_V  = sqrt(om["eta.v",  "eta.v"]),
    sd_eta_Cl = sqrt(om["eta.cl", "eta.cl"]),
    sigma = th[["add.sd"]])
}

ref <- c(Ka_pop = 1.50, V_pop = 32, Cl_pop = 3.0,
         sd_eta_Ka = 0.5, sd_eta_V = 0.1, sd_eta_Cl = 0.3, sigma = NA)

tab <- data.frame(
  `PinheiroBates` = ref,
  `dMod_FOCEI`    = dmod_est[names(ref)],
  `nlmixr2_FOCEI` = nx_estimates(nx_focei$fit)[names(ref)],
  check.names = FALSE)
if (!inherits(nx_saem, "try-error"))
  tab$`nlmixr2_SAEM` <- nx_estimates(nx_saem$fit)[names(ref)]

cat("\n================= Population estimates =================\n")
print(round(tab, 3))

cat("\n================= Runtime (seconds) ==================\n")
rt <- c(dMod_FOCEI    = t_dmod,
        nlmixr2_FOCEI = nx_focei$elapsed,
        nlmixr2_SAEM  = if (!inherits(nx_saem, "try-error"))
                          nx_saem$elapsed else NA_real_)
print(round(rt, 2))

## OFV: cross-tool objective values are NOT on a common scale here and must
## not be differenced. dMod models log(DV) with an additive residual (its -2LL
## lives on the log scale); nlmixr2's lnorm() likelihood lives on DV and differs
## by the change-of-variables Jacobian 2*sum(log DV); and nlmixr2's SAEM objf
## carries its own Gaussian-quadrature constant. The population estimates above
## are the meaningful cross-tool validation. The native values plus the exact
## bridge term are printed so the reader can reconcile deliberately.
jac <- 2 * sum(log(Theoph_pos$conc))
cat("\n================= Objective values (native scales) ===\n")
cat("  ! Not directly comparable across tools -- see the note below.\n")
cat(sprintf("  dMod   -2LL (log-DV scale)                 : %10.2f\n", fit_dmod$value))
cat(sprintf("  dMod   -2LL bridged to DV scale (+2*SlogDV): %10.2f\n",
            fit_dmod$value + jac))
cat(sprintf("  dMod   value_nonmem (- N_obs*log 2pi)      : %10.2f\n", fit_dmod$value_nonmem))
cat(sprintf("  nlmixr2 FOCEI objf                         : %10.2f\n",
            as.numeric(nx_focei$fit$objf)))
if (!inherits(nx_saem, "try-error"))
  cat(sprintf("  nlmixr2 SAEM  objf                         : %10.2f\n",
              as.numeric(nx_saem$fit$objf)))
cat(sprintf("\n  (Jacobian 2*sum(log DV) = %.2f; dMod models log(DV), nlmixr2\n", jac))
cat("   models DV via lnorm(), so their raw -2LLs differ by this term plus\n")
cat("   each tool's own additive OFV constant. Validate on estimates, not OFV.)\n")

invisible(list(dmod = fit_dmod, nlmixr2_focei = nx_focei$fit,
               nlmixr2_saem = if (!inherits(nx_saem, "try-error")) nx_saem$fit else NULL,
               estimates = tab, runtime = rt))
