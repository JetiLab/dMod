## WS4 capstone: clustered cell-line identification on the Bachmann JAK2/STAT5
## model (Bachmann et al., Mol Syst Biol 2011), 25 states / 34 reactions.
##
## We take the real Bachmann ODE (imported from PEtab) and manufacture 8 synthetic
## cell lines that differ ONLY in one core signalling rate, STAT5ActJAK2 (the
## STAT5-by-JAK2 phosphorylation rate, which drives the pSTAT5 readout), grouped
## into a KNOWN 3-cluster structure. We simulate pSTAT5 / CIS data for the 8 lines,
## then ask sparsify() to recover the grouping by exact-marginal model
## selection over the cluster path, exercising the adaptive quadrature (WS1, the
## anchored dimension reaches n-1 = 7) and the C++ marginal kernel (WS5).
##
## Run from the repo root (needs the sibling PEtab tree + libsbml Python):
##   Rscript bench/bench_laplace_bachmann.R
Sys.setenv(DMOD_LIBSBML_PYTHON = "/home/simon/.virtualenvs/libsbml/bin/python3")
suppressMessages(devtools::load_all(".", quiet = TRUE))
owd <- setwd(tempdir())            # keep compiled artifacts out of the repo
set.seed(20260723)

yaml <- Sys.getenv("DMOD_BACHMANN_YAML",
  "/home/simon/Documents/Projects/PEtab2dmod/BenchmarkModels/Bachmann_MSB2011/Bachmann_MSB2011.yaml")
pe <- importPEtab(yaml, solver = "CppODE", compile = TRUE, cores = 4L)
m  <- pe$odemodel

## ---- natural-scale nominal parameter values (from the Bachmann SBML) ----------
sbml <- c(
  CISEqc=432.860413, CISInh=785269991, CISRNADelay=0.14477776, CISRNAEqc=1,
  CISRNATurn=1000, CISTurn=0.0083988695, CISEqcOE=0.53026445, EpoRActJAK2=0.26730485,
  EpoRCISInh=1e6, EpoRCISRemove=5.4298069, JAK2ActEpo=633167.43, JAK2EpoRDeaSHP1=142.72332,
  SHP1ActEpoR=0.001, SHP1Dea=0.0081622049, SHP1ProOE=2.8256815, SOCS3EqcOE=0.67916552,
  SOCS3Eqc=173.64470, SOCS3Inh=10.407865, SOCS3RNADelay=1.0645845, SOCS3RNAEqc=1,
  SOCS3RNATurn=0.0083091764, SOCS3Turn=10000, STAT5ActEpoR=38.995799, STAT5ActJAK2=0.078106886,
  STAT5Exp=0.074515082, STAT5Imp=0.026886508, init_EpoRJAK2=3.9762237, init_SHP1=26.725116,
  init_STAT5=79.753640, cyt=0.4, nuc=0.275)
state_init <- c(EpoRJAK2 = 3.9762237, STAT5 = 79.753640, SHP1 = 26.725116)
EPO_DOSE   <- 2.5e-5              # strong Epo stimulation
CAND       <- "STAT5ActJAK2"
cand_nat   <- sbml[[CAND]]

## ---- clean observation of states (no PEtab observable machinery) --------------
x <- Xs(m, compile = TRUE)
states <- getParameters(x)
states <- states[!(states %in% names(sbml)) & !(states %in% c("Epo"))]  # the 25 state inits
g <- Y(eqnvec(obs_pSTAT5 = "log(pSTAT5 + npSTAT5 + 1)", obs_CIS = "log(CIS + 1)"),
       x, compile = TRUE, modelname = "bach_obs", attach.input = FALSE)
e <- Y(eqnvec(obs_pSTAT5 = "sigma", obs_CIS = "sigma"),
       g, compile = TRUE, modelname = "bach_err", attach.input = FALSE)

## ---- parameter transformation: constants everywhere except the candidate ------
inner <- unique(c(getParameters(x), getParameters(g), getParameters(e)))
trafo <- setNames(rep("0", length(inner)), inner)                 # states default 0
for (nm in inner) {
  if (nm %in% names(state_init))      trafo[nm] <- format(state_init[[nm]], digits = 12)
  else if (nm %in% names(sbml))       trafo[nm] <- format(sbml[[nm]], digits = 12)
  else if (nm %in% states)            trafo[nm] <- "0"            # other states start at 0
}
trafo["Epo"]        <- format(EPO_DOSE, digits = 12)
trafo[CAND]         <- "exp(log_STAT5ActJAK2 + eta_STAT5ActJAK2)" # the clustered candidate
trafo["sigma"]      <- "exp(log_sigma)"
trafo <- as.eqnvec(trafo)

lines <- LETTERS[1:8]
tab <- data.frame(eta_STAT5ActJAK2 = paste0("eta_STAT5ActJAK2_", lines),
                  row.names = lines, stringsAsFactors = FALSE)
p8  <- P(branch(trafo, table = tab, apply = "insert"), method = "explicit",
         modelname = "bach_p", compile = TRUE)
prd <- g * x * p8

## ---- ground truth: 3 clusters over the 8 lines --------------------------------
## {A,B,C} high, {D,E,F} baseline, {G,H} low  (log-offsets +1.2 / 0 / -1.2)
truth_eta <- c(A=1.2, B=1.2, C=1.2, D=0, E=0, F=0, G=-1.2, H=-1.2)
truth_grp <- list(c("A","B","C"), c("D","E","F"), c("G","H"))
st <- c(log_STAT5ActJAK2 = log(cand_nat), log_sigma = log(0.08))
times <- c(0, 3, 6, 10, 15, 20, 30, 45, 60, 90, 120)

pars_all <- c(st["log_STAT5ActJAK2"],
              setNames(truth_eta[lines], paste0("eta_STAT5ActJAK2_", lines)))
sim <- prd(times, c(pars_all, log_sigma = st[["log_sigma"]]), deriv = FALSE)

cat("== sanity: pSTAT5 trajectory range per line (should differ by cluster) ==\n")
for (s in lines) {
  tr <- sim[[s]]
  cat(sprintf("  line %s: obs_pSTAT5 range [%.3f, %.3f]\n", s,
              min(tr[, "obs_pSTAT5"]), max(tr[, "obs_pSTAT5"])))
}

## ---- simulate noisy data for the 8 lines --------------------------------------
dl_rows <- do.call(rbind, lapply(lines, function(s) {
  tr <- sim[[s]]
  do.call(rbind, lapply(c("obs_pSTAT5", "obs_CIS"), function(ob)
    data.frame(name = ob, time = tr[, "time"],
               value = tr[, ob] + rnorm(nrow(tr), 0, 0.08),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)))
}))
dlist <- as.datalist(dl_rows)

## ---- clustered selection: recover the grouping --------------------------------
pen  <- penaltyL1(CAND, method = "clustered", subjects = lines)
obj  <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) + constraintL1(pen)
init <- emInit(c(st, log_sigma = st[["log_sigma"]]), pen, lambda = 5)

cat("\n== regSelectCluster on the 8 Bachmann cell lines ==\n")
t0 <- Sys.time()
sel <- sparsify(obj, init, fixed = c(log_sigma = st[["log_sigma"]]),
                        verbose = TRUE)
cat(sprintf("  regSelectCluster wall-clock: %.1fs\n", as.numeric(Sys.time()-t0, units="secs")))

rec <- sel$perParam[[paste0("eta_", CAND)]]
cat("\n== recovered vs truth ==\n")
cat("  truth   : {A,B,C} {D,E,F} {G,H}\n")
cat(sprintf("  recovered: %s  (G = %d)\n",
    paste(vapply(rec$clusters, paste, "", collapse=","), collapse=" "), rec$G))
norm_grp <- function(P) sort(vapply(lapply(P, sort), paste, "", collapse=","))
ok <- identical(norm_grp(rec$clusters), norm_grp(truth_grp))
cat(sprintf("  MATCH the true 3-cluster structure: %s\n", ok))
print(rec$chain)
setwd(owd)
