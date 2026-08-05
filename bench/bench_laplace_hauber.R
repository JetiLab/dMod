## Head-to-head: marginal-likelihood clustered selection (regSelectCluster) vs the
## classical Hauber-et-al-2023 penalty scan, on the same synthetic Becker data
## (6 individuals, three ground-truth groupings; see inst/examples/EMclustered.R).
##
## Claim to demonstrate: the per-parameter marginal selector recovers the full
## per-parameter grouping (all three parameters) via the Occam factor, at a small
## FRACTION of the full-ODE fit count. The classical single-lambda scan, given
## every fair advantage (fine grid, relaxed post-lasso BIC, per-parameter
## selection), recovers only the well-identified granular parameter (kt) and
## either over-selects or misses the others, at many times the fit count.
##   ours   : one all-separate structural fit, then the whole cluster path scored
##            on the LINEARISED data (no ODE per grouping / per lambda).
##   Hauber : a lambda grid, a FULL-ODE fusion (MAP) fit at every grid point;
##            per candidate parameter, the distinct clusterings its coordinate
##            visits along the path are the models, each scored by BIC on a
##            relaxed unpenalised refit (the fair post-lasso baseline). The two
##            classical failure modes are made transparent by the printed
##            per-parameter visited-G path:
##              - kon (truth all-shared): G=1 IS visited, but BIC over-selects
##                along the sloppy kon/koff binding degeneracy (a hard parameter
##                count cannot see the near-flat direction the Occam factor does).
##              - ke  (truth two groups): the single global lambda collapses the
##                path past the balanced 2-split, so it is never a candidate.
##            The marginal selector, per parameter with the covariance-aware
##            score, recovers all three.
##
## Run from the repo root:  Rscript bench/bench_laplace_hauber.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
owd <- setwd(tempdir())            # keep compiled artefacts out of the repo
set.seed(2)

NLAMBDA <- as.integer(Sys.getenv("DMOD_HAUBER_NLAMBDA", "40"))  # scan resolution (fine: fair to Hauber)
OURFITS <- as.integer(Sys.getenv("DMOD_HAUBER_FITS", "8"))      # our multistart

## ---- Becker Epo receptor network, inline (no SBML / Python) ------------------
reactions <- NULL
reactions <- addReaction(reactions, "Epo + EpoR", "Epo_EpoR",
                         "Epo * EpoR * kon / init_Epo", "binding")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo + EpoR", "Epo_EpoR * koff", "unbinding")
reactions <- addReaction(reactions, "", "EpoR", "4 * init_Epo * init_EpoR_rel * kt", "synthesis")
reactions <- addReaction(reactions, "EpoR", "", "EpoR * kt", "turnover")
reactions <- addReaction(reactions, "Epo_EpoR", "Epo_EpoR_i", "Epo_EpoR * ke", "internalisation")
reactions <- addReaction(reactions, "Epo_EpoR_i", "Epo + EpoR", "Epo_EpoR_i * kex", "recycling")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_i", "Epo_EpoR_i * kdi", "intra degr.")
reactions <- addReaction(reactions, "Epo_EpoR_i", "dEpo_e", "Epo_EpoR_i * kde", "extra degr.")

m <- odemodel(reactions, modelname = "beckerh", backend = "cppDE",
              compile = FALSE, deriv2 = TRUE)
x <- Xs(m, compile = FALSE)
g <- Y(eqnvec(y_ext = "log(Epo + dEpo_e + 1)", y_mem = "log(Epo_EpoR + 1)",
              y_int = "log(Epo_EpoR_i + dEpo_i + 1)"),
       x, modelname = "beckerh_obs", compile = FALSE, deriv2 = TRUE, attach.input = FALSE)
e <- Y(eqnvec(y_ext = "sigma", y_mem = "sigma", y_int = "sigma"), g,
       modelname = "beckerh_err", compile = FALSE, deriv2 = TRUE, attach.input = FALSE)

dose  <- 1347.49
lines <- c("A", "B", "C", "D", "E", "F")
cand  <- c("kon", "ke", "kt")
trafo <- eqnvec(
  kon = "exp(log_kon + eta_kon)", ke = "exp(log_ke + eta_ke)",
  kt  = "exp(log_kt + eta_kt)",   koff = "exp(log_koff)",
  kex = "exp(log_kex)", kdi = "exp(log_kdi)", kde = "exp(log_kde)",
  init_Epo = "dose", init_EpoR_rel = "exp(log_EpoR_rel)",
  Epo = "dose", EpoR = "4 * dose * exp(log_EpoR_rel)",
  Epo_EpoR = "0", Epo_EpoR_i = "0", dEpo_i = "0", dEpo_e = "0", sigma = "exp(log_sigma)")
tab <- data.frame(eta_kon = paste0("eta_kon_", lines), eta_ke = paste0("eta_ke_", lines),
                  eta_kt = paste0("eta_kt_", lines), dose = dose,
                  row.names = lines, stringsAsFactors = FALSE)
p   <- P(branch(trafo, table = tab, apply = "insert"), method = "explicit",
         modelname = "beckerh_p", compile = FALSE, deriv2 = TRUE)
prd <- g * x * p
compile(prd, e, cores = 4)

## ---- ground truth: kon shared (G1), ke split (G2), kt three pairs (G3) -------
st <- c(log_kon = log(0.1513), log_koff = log(0.08069), log_kt = log(0.01601),
        log_ke = log(0.05546), log_kex = log(0.000578), log_kdi = log(0.001273),
        log_kde = log(0.01199), log_EpoR_rel = log(0.09205), log_sigma = log(0.05))
truth <- rbind(eta_kon = c(A = 0, B = 0, C = 0, D = 0, E = 0, F = 0),
               eta_ke  = c(A = 0.8, B = 0.8, C = 0.8, D = -0.8, E = -0.8, F = -0.8),
               eta_kt  = c(A = 0.7, B = 0.7, C = 0, D = 0, E = -0.7, F = -0.7))
truth_grp <- list(kon = list(lines), ke = list(c("A","B","C"), c("D","E","F")),
                  kt = list(c("A","B"), c("C","D"), c("E","F")))
eta_true <- unlist(lapply(cand, function(pp)
  stats::setNames(truth[paste0("eta_", pp), ], paste0("eta_", pp, "_", lines))))
times <- c(3, 6, 12, 20, 30, 45, 60, 90, 120, 180, 240, 300)
sim <- prd(times, c(st, eta_true), deriv = FALSE)
dl_rows <- do.call(rbind, lapply(lines, function(s) {
  pr <- sim[[s]]
  do.call(rbind, lapply(c("y_ext", "y_mem", "y_int"), function(nm)
    data.frame(name = nm, time = pr[, "time"], value = pr[, nm] + rnorm(nrow(pr), 0, 0.03),
               sigma = NA_real_, condition = s, stringsAsFactors = FALSE)))
}))
dlist <- as.datalist(dl_rows)

pen   <- penaltyL1(cand, method = "clustered", subjects = lines)
obj   <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE) + constraintL1(pen)
objD  <- normL2(dlist, prd, errmodel = e, use.bessel = FALSE)   # unpenalised, for BIC
fixed <- st[c("log_koff", "log_kex", "log_kdi", "log_kde", "log_EpoR_rel")]
init  <- emInit(c(st[paste0("log_", cand)], log_sigma = log(0.03)), pen, lambda = 5)
Ndat  <- nrow(dl_rows)

norm_grp <- function(P) sort(vapply(lapply(P, sort), paste, "", collapse = ","))
match_all <- function(rec)
  all(vapply(cand, function(pp)
    identical(norm_grp(rec[[paste0("eta_", pp)]]), norm_grp(truth_grp[[pp]])), NA))
show_grp <- function(cl)
  paste(vapply(cl, paste, "", collapse = ","), collapse = "|")

## ---- (A) ours: marginal clustered selection ---------------------------------
cat("== (A) regSelectCluster (marginal-likelihood) ==\n")
tA <- system.time(
  selA <- sparsify(obj, init, fixed = fixed, fits = OURFITS, cores = 4,
                           sd = 0.3, verbose = FALSE))[["elapsed"]]
recA <- lapply(selA$perParam, function(z) z$clusters)
fitsA <- OURFITS                                   # full-ODE structural fits used

## ---- (B) Hauber: full-ODE fusion scan, then PER-PARAMETER relaxed-refit BIC --
## Each grid point is a full nonlinear clustered MAP fit (regFit, lambda fixed).
## From the joint path we read off, PER candidate parameter, the distinct
## clusterings that parameter's coordinate visits, and BIC-select each parameter
## independently on a RELAXED (unpenalised, reference-encoded) refit: within each
## cluster the deviation is shared, cluster 1 is the reference (deviation 0), and
## clusters 2..G plus the structural log-means and log_sigma are free. Per-param
## selection is the fair mirror of regSelectCluster: a single global lambda scan
## cannot mix granularities (kon fully shared while kt splits three ways), so
## selecting per parameter gives the scan its best shot at the truth. Scoring the
## penalised fit directly over-selects (its deviations are lambda-shrunk and not
## comparable across the grid), hence the relaxed refit.
cat(sprintf("== (B) Hauber penalty scan (%d full-ODE fits) + per-param relaxed refits ==\n", NLAMBDA))
lamGrid   <- 10^seq(-2, 4, length.out = NLAMBDA)   # span all-separate to all-fused
eta_names <- as.vector(pen$subjectEtas)
structNms <- setdiff(names(init), c(names(fixed), pen$lambdaName))  # free log-means + log_sigma
nStruct   <- length(structNms)

## unpenalised reference-encoded refit of one grouping -> data -2logL + free df
relaxFit <- function(grp) {
  devNames <- unlist(lapply(cand, function(pp) {
    cl <- grp[[paste0("eta_", pp)]]
    if (length(cl) <= 1L) character(0) else paste0("dev_", pp, "_", seq_along(cl)[-1])
  }))
  toEta <- function(par) {
    ev <- setNames(numeric(length(eta_names)), eta_names)
    for (pp in cand) {
      cl <- grp[[paste0("eta_", pp)]]
      for (gi in seq_along(cl)) {
        val <- if (gi == 1L) 0 else par[[paste0("dev_", pp, "_", gi)]]
        ev[paste0("eta_", pp, "_", cl[[gi]])] <- val
      }
    }
    ev
  }
  fn <- function(par)
    objD(c(par[structNms], fixed, toEta(par)[eta_names]), deriv = FALSE)$value
  gr <- function(par) {
    g <- objD(c(par[structNms], fixed, toEta(par)[eta_names]), deriv = TRUE)$gradient
    gd <- vapply(devNames, function(dn) {
      s  <- strsplit(dn, "_", fixed = TRUE)[[1]]
      cl <- grp[[paste0("eta_", s[2])]][[as.integer(s[3])]]
      sum(g[paste0("eta_", s[2], "_", cl)])
    }, 0.0)
    c(g[structNms], gd)
  }
  par0 <- c(init[structNms], setNames(rep(0, length(devNames)), devNames))
  lo   <- c(init[structNms] - 6, setNames(rep(-6, length(devNames)), devNames))
  hi   <- c(init[structNms] + 6, setNames(rep( 6, length(devNames)), devNames))
  opt  <- optim(par0, fn, gr, method = "L-BFGS-B", lower = lo, upper = hi,
                control = list(maxit = 300))
  list(m2ll = opt$value, df = nStruct + length(devNames))
}

sepGroup <- lapply(lines, function(s) s)           # every subject its own singleton

tB <- system.time({
  ## scan: collect, per parameter, the distinct clusterings its coordinate visits
  visited <- setNames(lapply(cand, function(pp) list()), cand)
  keys    <- setNames(lapply(cand, function(pp) character(0)), cand)
  nScan   <- 0L
  for (lam in lamGrid) {
    ini <- init; ini[pen$lambdaName] <- lam
    f <- try(suppressMessages(EM(obj, ini, fixed = fixed, method = "focei",
             control = list(estimateLambda = FALSE, maxOuter = 25L), verbose = FALSE)),
             silent = TRUE)
    nScan <- nScan + 1L
    if (inherits(f, "try-error")) next
    for (pp in cand) {
      cl <- f$clusters[[paste0("eta_", pp)]]
      k  <- paste(norm_grp(cl), collapse = "|")
      if (!k %in% keys[[pp]]) {
        keys[[pp]]    <- c(keys[[pp]], k)
        visited[[pp]] <- c(visited[[pp]], list(cl))
      }
    }
  }
  ## per parameter: BIC-select over its own clusterings, other params all-separate
  nRefit <- 0L
  recB   <- setNames(vector("list", length(cand)), paste0("eta_", cand))
  for (pp in cand) {
    cands_pp <- visited[[pp]]
    bic <- vapply(cands_pp, function(cl) {
      grp <- setNames(rep(list(sepGroup), length(cand)), paste0("eta_", cand))
      grp[[paste0("eta_", pp)]] <- cl                # target param grouped, rest separate
      rf <- relaxFit(grp)
      rf$m2ll + log(Ndat) * rf$df
    }, 0.0)
    nRefit <- nRefit + length(cands_pp)
    best   <- which.min(bic)
    recB[[paste0("eta_", pp)]] <- cands_pp[[best]]
    cat(sprintf("   %-4s path G visited: {%s}  -> BIC-min G=%d\n", pp,
                paste(vapply(cands_pp, length, 0L), collapse = ","),
                length(cands_pp[[best]])))
  }
})[["elapsed"]]
fitsB <- nScan + nRefit                             # scan fits + relaxed refits

## ---- head-to-head -----------------------------------------------------------
cat("\n================= head-to-head (Becker, 6 individuals) =================\n")
cat(sprintf("%-6s | %-22s | %-22s | %-22s\n", "param", "truth", "ours (marginal)", "Hauber (scan+BIC)"))
for (pp in cand) {
  ec <- paste0("eta_", pp)
  cat(sprintf("%-6s | %-22s | %-22s | %-22s\n", pp,
              show_grp(truth_grp[[pp]]), show_grp(recA[[ec]]), show_grp(recB[[ec]])))
}
cat("------------------------------------------------------------------------\n")
cat(sprintf("recovers truth   : ours = %s    Hauber = %s\n",
            match_all(recA), match_all(recB)))
cat(sprintf("full-ODE fits    : ours = %-4d   Hauber = %-4d  (%.1fx fewer)\n",
            fitsA, fitsB, fitsB / fitsA))
cat(sprintf("wall-clock (s)   : ours = %-6.1f Hauber = %-6.1f (%.1fx faster)\n",
            tA, tB, tB / tA))
cat("========================================================================\n")
setwd(owd)
