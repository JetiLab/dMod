# What a symmetry does to a profile likelihood.
#
# P <-> pP with a log-scale observation y = log2(pP) + s. The scaling s is fixed, so
# the remaining non-identifiability is the polynomial direction symmetryDetection()
# finds: along its orbit the likelihood does not move, and the profile is flat until
# the orbit runs into {k_p = 0}, a Darboux surface the flow cannot cross -- that wall
# is where the profile turns up again. symmetryReduction() reparameterises the
# direction away, and the same profile closes.
#
# Meant to be stepped through in RStudio, top to bottom. Every plot goes to the plot
# pane; nothing is written to disk.

suppressMessages(devtools::load_all(".", quiet = TRUE))
setwd(tempdir())          # generated sources and shared objects land here

cores <- detectFreeCores()


# ---- the model and its symmetry --------------------------------------------

reactions <- eqnlist() |>
  addReaction("P",  "pP", "k_p*P",  "phosphorylation") |>
  addReaction("pP", "P",  "k_d*pP", "dephosphorylation")

# rational form for the analysis -- the engine works over Q, not over log2
res <- symmetryDetection(reactions, eqnvec(y = "s*pP"),
                         method = "observability", reconstruct = TRUE)
summary(res)

red <- symmetryReduction(res, fixed = "s")
summary(red)


# ---- the two prediction chains ---------------------------------------------

# The ODE itself is integrated in log10: log10Transform() rewrites dot x = f(x) as
# dot x_l10 = f / (log(10) * 10^x_l10). A power of ten is positive for every real
# exponent, so a state can never go negative -- which matters here, because the
# profile walks k_p towards the Darboux surface {k_p = 0} where the linear system
# would happily undershoot. The trafos below keep using exp10(), which reads better
# and is fine there; only the ODE and the observable need a differentiable form.
fLog <- log10Transform(as.eqnvec(reactions))
fLog

x <- Xs(odemodel(fLog, modelname = "x", compile = FALSE), condition = "C1",
        optionsOde  = list(atol = 1e-8, rtol = 1e-6, maxsteps = 1e7),
        optionsSens = list(atol = 1e-8, rtol = 1e-6, maxsteps = 1e7))
# the observable is now linear in the state: log2(10^pP_l10) = pP_l10*log2(10).
# Written with 10^ rather than exp10() for the same reason log10Transform() uses it:
# exp10 has no entry in R's derivatives table, so Y() could not differentiate it
g <- Y(eqnvec(y = "log2(10^(pP_l10)) + s"), f = x, attach.input = FALSE,
       modelname = "g", compile = FALSE)

states    <- names(as.eqnvec(reactions))
innerpars <- getParameters(g, x)
linpars   <- c(states, setdiff(innerpars, paste0(states, "_l10")))

# The fit lives on a log10 scale, so the chain is two transformations: exp10 back to
# the linear model parameters, then log10 again for the two initial values the
# transformed ODE asks for. The round trip costs nothing in the generated code, and
# it keeps red$trafo -- which is stated in the LINEAR coordinates -- insertable
# exactly where it belongs.
p.lin <- eqnvec() |>
  define("x~x", x = linpars) |>
  define("s~0") |>                     # the gauge symmetryReduction() fixed
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "plin", compile = FALSE)

# the same, with the symmetry direction reparameterised away
p.red <- eqnvec() |>
  define("x~x", x = linpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  define("s~0") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "pred", compile = FALSE)

# linear model parameters -> what the log10 ODE takes: rates stay, states go to log10
p.log <- eqnvec() |>
  define("x~x", x = setdiff(linpars, states)) |>
  define("x_l10 ~ log10(x)", x = states) |>
  P(modelname = "plog", compile = FALSE)

compile(x, g, p.lin, p.log, p.red, cores = cores)

prd     <- g * x * p.log * p.lin
prd.red <- g * x * p.log * p.red


# ---- data ------------------------------------------------------------------

truth <- log10(c(P = 1, pP = 0.2, k_p = 0.3, k_d = 0.05, sigma_pP = 0.3))
times <- seq(0, 10, len = 300)
nrep  <- 12

set.seed(1111)
sim <- as.data.frame(prd(seq(0, 10, by = 0.8), truth, deriv = FALSE))

data <- data.frame(
  name = "y", time = rep(sim$time, each = nrep),
  value = rep(sim$value, each = nrep) +
    rnorm(nrep * nrow(sim), 0, 10^(truth[["sigma_pP"]])),
  sigma = NA, condition = "C1") |>
  reduceReplicates() |>
  as.datalist(split.by = "condition")

plot(data)
plot(prd(times, truth, deriv = FALSE), data)


# ---- fit and profile, both chains ------------------------------------------

stepControl <- list(stepsize = 1e-8, min = 1e-8, max = Inf,
                    atol = 1e-3, rtol = 1e-3, limit = 1e3)
optControl  <- list(rinit = 0.1, rmax = 10, iterlim = 2e3)

pouter <- structure(rep(-1, length(getParameters(prd))), names = getParameters(prd))
obj    <- normL2(data, prd)
fit    <- mstrust(obj, pouter, rinit = 0.1, rmax = 10, sd = 4, fits = 100,
                  iterlim = 1e3, cores = cores, studyname = "full")

plotValues(as.parframe(fit), tol = 0.1, value < 1e4)
bestfit <- as.parvec(as.parframe(fit))
plot(prd(times, bestfit, deriv = FALSE), data)

prof <- profile(obj, bestfit, whichPar = names(bestfit), method = "optimize",
                stepControl = stepControl, optControl = optControl,
                limits = c(lower = -10, upper = 10), cores = cores)
plotProfilesAndPaths(prof, names(bestfit))


pouter.red <- structure(rep(-1, length(getParameters(prd.red))),
                        names = getParameters(prd.red))
obj.red    <- normL2(data, prd.red)
fit.red    <- mstrust(obj.red, pouter.red, rinit = 0.1, rmax = 10, sd = 4, fits = 100,
                      iterlim = 1e3, cores = cores, studyname = "reduced")

plotValues(as.parframe(fit.red), tol = 0.1, value < 1e4)
bestfit.red <- as.parvec(as.parframe(fit.red))
plot(prd.red(times, bestfit.red, deriv = FALSE), data)

prof.red <- profile(obj.red, bestfit.red, whichPar = names(bestfit.red),
                    method = "optimize",
                    stepControl = stepControl, optControl = optControl,
                    limits = c(lower = -10, upper = 10), cores = cores)
plotProfilesAndPaths(prof.red, names(bestfit.red))


# the reduction removes a direction, not a degree of freedom of the fit: the two
# optima agree, only the reduced profile is finite
cat(sprintf("optimum full %.6f, reduced %.6f, difference %.2e\n",
            obj(bestfit)$value, obj.red(bestfit.red)$value,
            abs(obj.red(bestfit.red)$value - obj(bestfit)$value)))
