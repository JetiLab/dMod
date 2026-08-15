# Profiles under a polynomial symmetry: flat along the orbit, walled at its
# reach ({k = 0} is a Darboux surface the orbit cannot cross).
# Rscript bench/symmetry_profile_flatness.R   -> symmetry_profile_flatness.pdf

suppressMessages(devtools::load_all(".", quiet = TRUE))
outdir <- normalizePath(".")
setwd(tempdir())

reactions <- eqnlist() |>
  addReaction("P", "pP", "k_p*P",  "phosphorylation") |>
  addReaction("pP", "P", "k_d*pP", "dephosphorylation")
observables <- eqnvec(y = "s*pP")   # rational form for the analysis

res <- symmetryDetection(reactions, observables, method = "observability", reconstruct = TRUE)

summary(res)

red <- reduceSymmetry(res, fixed = "s")

print(red)

m <- odemodel(reactions, modelname = "x", compile = FALSE)
x <- Xs(m, condition = "C1", optionsOde = list(atol = 1e-8, rtol = 1e-6, maxsteps = 1e7), 
        optionsSens = list(atol = 1e-8, rtol = 1e-6, maxsteps = 1e7))
g <- Y(eqnvec(y = "log2(pP) + s"), f = x, condition = NULL,
       attach.input = FALSE, modelname = "g", compile = FALSE)
e <- Y(c(y = "sigma_pP"), f = g, attach.input = FALSE, modelname = "e", compile = FALSE)
innerpars <- getParameters(g, x, e)

p <- eqnvec() |> 
  define("x~x", x = innerpars) |> 
  define("s~0") |> 
  insert("x~exp10(x)", x = .currentSymbols) |> 
  P(modelname = "p", compile = FALSE)

p.red <- eqnvec() |> 
  define("x~x", x = innerpars) |> 
  insert("x~y", x = names(red$trafo), y = red$trafo) |> 
  define("s~0") |> 
  insert("x~exp10(x)", x = .currentSymbols) |> 
  P(modelname = "pred", compile = FALSE)

compile(x, g, e, p, p.red, cores = 6)
prd <- g*x*p

truth <- log10(c(P = 1, pP = 0.2, k_p = 0.3, k_d = 0.05, sigma_pP = 0.1))

timesD <- seq(0, 10, by = 0.8)
nrep <- 6

set.seed(5555)

sim <- prd(timesD, truth, deriv = FALSE) |> as.data.frame()

data <- as.datalist(data.frame(
  name = "y", time = rep(sim$time, each = nrep),
  value = rep(sim$value, each = nrep) +
    rnorm(nrep * nrow(sim), 0, 10^(truth[["sigma_pP"]])),
  sigma = NA, condition = "C1"), split.by = "condition")

times <- seq(0, 10, len = 300)
plot(prd(times, truth, deriv = FALSE), data)

pouter <- structure(rep(-1, length(getParameters(prd))), names = getParameters(prd))

obj <- normL2(data, prd, errmodel = e)
fit <- mstrust(obj, pouter, rinit = 0.1, rmax = 10, sd = 4, 
               studyname = "msfit", iterlim = 1e3, fits = 100, cores = 20)

plotValues(as.parframe(fit), tol = 0.1, value < 1e4)

bestfit <- as.parframe(fit) |> as.parvec()

plot(prd(times, bestfit, deriv = FALSE), data)

prof <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate",
                algoControl = list(reoptimize = TRUE),
                stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-3, rtol = 1e-3, limit = 1e3),
                limits = c(lower = -100, upper = 100), cores = length(fit$argument))

plotProfilesAndPaths(prof, names(bestfit))


prd <- g*x*p.red

pouter <- structure(rep(-1, length(getParameters(prd))), names = getParameters(prd))

obj <- normL2(data, prd, errmodel = e)
fit <- mstrust(obj, pouter, rinit = 0.1, rmax = 10, sd = 4, 
               studyname = "msfit", iterlim = 1e3, fits = 100, cores = 20)


plotValues(as.parframe(fit), tol = 0.1, value < 1e4)

bestfit <- as.parframe(fit) |> as.parvec()

plot(prd(times, bestfit, deriv = FALSE), data)


prof.red <- profile(obj, bestfit, whichPar = names(bestfit), method = "integrate",
                algoControl = list(reoptimize = TRUE),
                limits = c(lower = -2, upper = 2), cores = length(fit$argument))

plotProfilesAndPaths(prof.red, names(bestfit))
