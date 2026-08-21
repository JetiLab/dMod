# What a symmetry does to a profile likelihood.
#
# P <-> pP observed through y = log2(s*pP), the log of the rational observable the
# engine analyses. First the two directions are flown both ways, as in the
# invariance bench: the states move, the observable does not. Then they are profiled
# -- the scaling is gauged away by fixing s, the curved direction is flat until its
# orbit hits the Darboux surface {k_p = 0}, and symmetryReduction() closes that
# profile.
#
# Step through in RStudio, top to bottom. Nothing is written to disk.

suppressMessages(devtools::load_all(".", quiet = TRUE))
library(ggplot2)
setwd(tempdir())

cores <- detectFreeCores()


# ---- the model and its symmetries -------------------------------------------

reactions <- eqnlist() |>
  addReaction("P",  "pP", "k_p*P",  "phosphorylation") |>
  addReaction("pP", "P",  "k_d*pP", "dephosphorylation")

observables <- eqnvec(y = "s*pP")          # rational: the engine works over Q

res <- symmetryDetection(reactions, observables,
                         method = "observability", reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates

red <- symmetryReduction(res, fixed = "s") # s is the gauge the fit fixes too
summary(red)


# ---- the prediction chain ---------------------------------------------------

# the ODE is integrated as it stands; only the parameters live on log10, which keeps
# every rate positive and puts {k_p = 0} at -Inf, where the profile walks
x <- Xs(odemodel(reactions, modelname = "x", compile = FALSE), condition = "C1",
        optionsOde  = list(atol = 1e-10, rtol = 1e-10, maxsteps = 1e7),
        optionsSens = list(atol = 1e-10, rtol = 1e-10, maxsteps = 1e7))
# attach.input: the states ride along for the flow picture; only the observable
# decides the verdict, and the objective ignores what the data has no entry for.
# log2 rather than exp10's counterpart: Y() differentiates R's table only.
g <- Y(eqnvec(y = "log2(s*pP)"), f = x, attach.input = TRUE,
       modelname = "g", compile = FALSE)

innerpars <- getParameters(g, x)

# the flow check needs every coordinate free, in the _l10 coordinates the flow is
# integrated in -- output and input meet without conversion
p.flow <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x ~ exp10(x_l10)", x = .currentSymbols) |>
  P(modelname = "pflow", compile = FALSE)

p.lin <- eqnvec() |>
  define("x~x", x = innerpars) |>
  define("s~1") |>                         # gauge: the log2 offset is 0
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "plin", compile = FALSE)

# red$trafo is stated in the linear coordinates, so it goes in before exp10. Every
# carrier here is certified positive (red$blocks[[1]]$carrierDomain); a
# "[real-valued]" one would have to stay linear.
p.red <- eqnvec() |>
  define("x~x", x = innerpars) |>
  insert("x~y", x = names(red$trafo), y = red$trafo) |>
  define("s~1") |>
  insert("x~exp10(x)", x = .currentSymbols) |>
  P(modelname = "pred", compile = FALSE)

# one flow per direction, with a sign switch: dz/deps = sgn*xi(z) runs the same
# compiled model both ways, so eps stays a forward integration either way. Its
# states ARE the moved _l10 coordinates; log10Transform() renames only those, so the
# rest enter linearly and keep their nominal value.
flow <- lapply(seq_along(res$symmetries), function(i) {
  fl <- log10Transform(res$symmetries[[i]]$generator)
  Xf(odemodel(as.eqnvec(setNames(paste0("sgn*(", fl, ")"), names(fl))),
              deriv = FALSE, modelname = paste0("flow", i), compile = FALSE),
     condition = "flow",
     optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
})

do.call(compile, c(list(x, g, p.flow, p.lin, p.red), flow, list(cores = cores)))

prd     <- g * x * p.flow                  # every coordinate free
prd.lin <- g * x * p.lin                   # s gauged
prd.red <- g * x * p.red                   # s gauged, curved direction removed


# ============================================================================
# the directions, flown
# ============================================================================

truth  <- c(P = 1, pP = 0.2, k_p = 0.3, k_d = 0.05, s = 1)
truthL <- setNames(log10(truth), paste0(names(truth), "_l10"))

times <- seq(0, 10, len = 300)
pred0 <- as.data.frame(prd(times, truthL, deriv = FALSE))   # the unmoved reference

direction <- 2                            # click through 1 .. length(flow)
# the two ways go unequally far: one runs a coordinate to zero at finite eps and the
# flow stops there, the other runs towards {k_p = 0} and only reaches it at
# infinity. Which side is which follows the sign of the returned generator --
# shorten whichever range stops.
eps  <- seq(0, 5,    len = 50)            # forward
epsb <- seq(0, 0.1541, len = 50)            # backward

mv  <- names(res$symmetries[[direction]]$generator)
mvL <- paste0(mv, "_l10")
z0  <- c(truthL[mvL], truth[setdiff(coords, mv)])

zf <- flow[[direction]](eps, c(z0, sgn =  1), deriv = FALSE)$flow
zb <- flow[[direction]](epsb, c(z0, sgn = -1), deriv = FALSE)$flow

# one orbit through the truth, indexed by signed eps. The axis is taken from the
# returned time columns, so unequal ranges -- or a run the solver cut short -- stay
# labelled correctly.
z <- rbind(zb[rev(seq_len(nrow(zb))), ], zf[-1, ])
z[, "time"] <- c(-rev(zb[, "time"]), zf[-1, "time"])

ggplot(wide2long(z), aes(time, value)) +
  geom_line() + geom_vline(xintercept = 0, linetype = 3) +
  facet_wrap(~ name, scales = "free_y") +
  labs(x = expression(epsilon), y = NULL) + theme_dMod(base_size = 9)

theta <- lapply(seq_len(nrow(z)), function(j) replace(truthL, mvL, z[j, mvL]))
pred  <- lapply(theta, function(th) as.data.frame(prd(times, th, deriv = FALSE)))

# the verdict: max deviation of the observable, relative to the RMS of the
# unmoved curve
obs <- pred0$name %in% names(observables)
rms <- ave(pred0$value[obs], pred0$name[obs], FUN = function(v) sqrt(mean(v^2)))

sapply(pred, function(d) max(abs(d$value[obs] - pred0$value[obs]) / rms))


# ---- the picture: the states move, the observable does not ------------------

long <- do.call(rbind, Map(function(d, e) transform(d, eps = e), pred, z[, "time"]))

# mark the observable and move it to the end
lev <- levels(long$name)
levels(long$name) <- ifelse(lev %in% names(observables), paste0(lev, " [observed]"), lev)
long$name <- factor(long$name, levels(long$name)[order(lev %in% names(observables))])

ggplot(long, aes(time, value, group = eps, colour = eps)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  scale_color_dMod_div(name = expression(epsilon)) +
  labs(title = paste0("P/pP, X", direction), x = "time", y = NULL) +
  theme_dMod(base_size = 9)


# ---- every direction, both ways, one line each ------------------------------

for (i in seq_along(flow)) {
  mv  <- names(res$symmetries[[i]]$generator)
  mvL <- paste0(mv, "_l10")
  z0  <- c(truthL[mvL], truth[setdiff(coords, mv)])
  for (sgn in c(1, -1)) {
    e  <- if (sgn > 0) eps else epsb
    zi <- try(flow[[i]](e, c(z0, sgn = sgn), deriv = FALSE)$flow, silent = TRUE)
    if (inherits(zi, "try-error")) {
      cat(sprintf("X%d  %-8s  eps %+.2f  escapes\n", i,
                  res$symmetries[[i]]$type, sgn * max(e)))
      next
    }
    d <- sapply(seq_len(nrow(zi)), function(j) {
      q <- as.data.frame(prd(times, replace(truthL, mvL, zi[j, mvL]), deriv = FALSE))
      max(abs(q$value[obs] - pred0$value[obs]) / rms)
    })
    cat(sprintf("X%d  %-8s on %-20s  eps %+.2f  moved %.2f  deviates %.2e\n", i,
                res$symmetries[[i]]$type,
                paste(res$symmetries[[i]]$support, collapse = ","),
                sgn * max(zi[, "time"]),
                max(abs(10^(zi[nrow(zi), mvL] - truthL[mvL]) - 1)), max(d)))
  }
}


# ============================================================================
# the same directions, profiled
# ============================================================================

# ---- data ------------------------------------------------------------------

truth.fit <- log10(truth[getParameters(prd.lin)])
sigma     <- 0.3                           # noise sd of the log2 observation
nrep      <- 12

set.seed(1111)
sim <- as.data.frame(prd.lin(seq(0, 10, by = 0.8), truth.fit, deriv = FALSE))
sim <- sim[sim$name == "y", ]              # attach.input carries the states too

data <- data.frame(
  name = "y", time = rep(sim$time, each = nrep),
  value = rep(sim$value, each = nrep) + rnorm(nrep * nrow(sim), 0, sigma),
  sigma = NA, condition = "C1") |>
  reduceReplicates() |>
  as.datalist(split.by = "condition")

plot(data)
plot(prd.lin(times, truth.fit, deriv = FALSE), data)


# ---- fit and profile, both chains ------------------------------------------

stepControl <- list(stepsize = 1e-8, min = 1e-8, max = Inf,
                    atol = 1e-3, rtol = 1e-3, limit = 1e3)
optControl  <- list(rinit = 0.1, rmax = 10, iterlim = 2e3)

pouter <- structure(rep(-1, length(getParameters(prd.lin))),
                    names = getParameters(prd.lin))
obj    <- normL2(data, prd.lin)
fit    <- mstrust(obj, pouter, rinit = 0.1, rmax = 10, sd = 4, fits = 100,
                  iterlim = 1e3, cores = cores, studyname = "full")

plotValues(as.parframe(fit), tol = 0.1, value < 1e4)
bestfit <- as.parvec(as.parframe(fit))
plot(prd.lin(times, bestfit, deriv = FALSE), data)

prof <- profile(obj, bestfit, whichPar = names(bestfit), method = "optimize",
                stepControl = stepControl, optControl = optControl,
                limits = c(lower = -10, upper = 10), cores = cores)
plotProfilesAndPaths(prof, names(bestfit))


# the fit in the new coordinates
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


# the reduction removes a direction, not a degree of freedom: the optima agree,
# only the reduced profile is finite
cat(sprintf("optimum full %.6f, reduced %.6f, difference %.2e\n",
            obj(bestfit)$value, obj.red(bestfit.red)$value,
            abs(obj.red(bestfit.red)$value - obj(bestfit)$value)))
