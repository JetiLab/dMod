# What a symmetry generator claims, checked numerically.
#
# The verification inside symmetryDetection() is algebraic: Taylor jet at t0, over
# GF(p), stopped at a saturation guard. Here the claim itself is tested -- flow the
# coordinates along X and the prediction g*x does not move -- over R, on the compiled
# model. It can only falsify: a surviving direction may still break at a larger eps,
# and a MISSED direction is invisible by construction.
#
# Two passes: the published NFkB model (Merkt, Timmer & Kaschek, Phys Rev E 92,
# 012920 (2015), supplement S.2), and the same check on an EGF cascade held at steady
# state by a `trafo`, with a dose event.
#
# Meant to be stepped through in RStudio, top to bottom. Every plot goes to the plot
# pane; nothing is written to disk.

library(dMod)
library(ggplot2)
setwd(tempdir())

# Throughout, coordinates live on a log10 scale and carry an "_l10" suffix. The
# generator flow is integrated in those coordinates -- log10Transform() rewrites
# dz/deps = xi(z) as dz_l10/deps = xi / (log(10) * 10^z_l10) -- and the biological
# system takes them back to the model with exp10() in its parameter transformation.
# Output and input therefore meet in the same coordinates, with nothing converted in
# between, and a power of ten is positive for every real exponent, so no coordinate
# can come out negative.



# ============================================================================
# NFkB
# ============================================================================

# The input u = pIKK is held constant, so it is simply a parameter -- carrying it as
# a state with u' = 0 is a detour that changes nothing (same rank, same coordinates,
# same generators). What matters is that u IS a coordinate: that is what puts the
# paper's Moebius transformation u -> u/(1 - eps*u), k0 -> k0 + eps into the space.
# Initial values are left free -- the honest choice for an identifiability analysis,
# and it makes every one of them a coordinate the flow has to move.
nfkb <- eqnvec(
  x1  = "k11*x10 - (k1*u/(1 + k0*u) + k1_0)*x1",
  x2  = "(k1*u/(1 + k0*u) + k1_0)*x1 - k2*x2",
  x3  = "k2*x2 - k3*x3",
  x4  = "k2*x2 - k4*x4",
  x5  = "k3*rho_vol*x3 - k5*x5",
  x6  = "k5*x5 - k10*x9*x6",
  x7  = "k6*x6 - k7*x7",
  x8  = "k8*x7 - k9*x8",
  x9  = "k9*rho_vol*x8 - k10*x9*x6",
  x10 = "k10*x9*x6 - k11*rho_vol*x10")

observables <- eqnvec(y1 = "s1*(x1 + x2 + x3) + I0cyt",
                      y2 = "s2*(x10 + x5 + x6) + I0nuc",
                      y3 = "s3*(x2 + x3)",
                      y4 = "s4*(x2 + x4)")

res <- symmetryDetection(nfkb, observables, method = "observability",
                         reconstruct = TRUE)
summary(res)

coords <- res$info$coordinates


# ---- the prediction the generators are a claim about -----------------------

x <- Xf(odemodel(nfkb, deriv = FALSE, modelname = "x", compile = FALSE),
        condition = "C1",
        optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7))
# attach.input = TRUE: one prediction carries the states next to the observables.
# Only the observables decide the verdict -- the states are meant to move.
g <- Y(observables, f = x, attach.input = TRUE, modelname = "g", compile = FALSE)

# the parameters arrive on a log10 scale, so p is where they come back to the model
p <- eqnvec() |>
  define("x~x", x = getParameters(g, x)) |>
  insert("x ~ exp10(x_l10)", x = .currentSymbols) |>
  P(modelname = "p", compile = FALSE)

# one flow per direction, integrated in the same coordinates: its states ARE the
# moved _l10 coordinates. Coordinates a direction leaves alone are not states at
# all; log10Transform() renames only the states, so those enter the right-hand side
# linearly and keep their nominal value.
flow <- lapply(seq_along(res$symmetries), function(i)
  Xf(odemodel(log10Transform(res$symmetries[[i]]$generator), deriv = FALSE,
              modelname = paste0("flow", i), compile = FALSE), condition = "flow",
     optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)))

do.call(compile, c(list(x, g, p), flow, list(cores = 10)))

prd <- g * x * p


# ---- flow one direction ----------------------------------------------------

truth <- c(x1 = 1, x2 = 0.1, x3 = 0.05, x4 = 0.05, x5 = 0.02, x6 = 0.03,
           x7 = 0.04, x8 = 0.2, x9 = 0.06, x10 = 0.07, u = 0.5,
           k0 = 0.4, k1 = 0.7, k1_0 = 0.05, k2 = 0.3, k3 = 0.25, k4 = 0.2,
           k5 = 0.35, k6 = 0.15, k7 = 0.1, k8 = 0.45, k9 = 0.12, k10 = 0.8,
           k11 = 0.22, rho_vol = 1.5, s1 = 1.1, s2 = 0.9, s3 = 1.3, s4 = 0.7,
           I0cyt = 0.2, I0nuc = 0.15)

truthL <- setNames(log10(truth), paste0(names(truth), "_l10"))

times <- seq(0, 60, len = 500)
eps   <- seq(0, 5, len = 100)

direction <- 2                      # click through this with 1 .. length(flow)
mv  <- names(res$symmetries[[direction]]$generator)
mvL <- paste0(mv, "_l10")
z   <- flow[[direction]](eps, c(truthL[mvL], truth[setdiff(coords, mv)]),
                         deriv = FALSE)$flow

plot(z) + xlab("eps")

theta <- lapply(seq_len(nrow(z)), function(j) replace(truthL, mvL, z[j, mvL]))
pred  <- lapply(theta, function(th) as.data.frame(prd(times, th, deriv = FALSE)))

# the verdict: max deviation of the OBSERVABLES, relative to the RMS of the
# unperturbed curve of the same observable -- one scale per curve, so a small
# observable is not drowned by a large one
obs <- pred[[1]]$name %in% names(observables)
rms <- ave(pred[[1]]$value[obs], pred[[1]]$name[obs], FUN = function(v) sqrt(mean(v^2)))

sapply(pred, function(d) max(abs(d$value[obs] - pred[[1]]$value[obs]) / rms))


# ---- the picture: the states move, the observables do not ------------------

long <- do.call(rbind, Map(function(d, e) transform(d, eps = e), pred, z[, "time"]))

# wide2long() levels `name` in the prediction's own column order; mark the
# observables and move them to the end
lev <- levels(long$name)
levels(long$name) <- ifelse(lev %in% names(observables), paste0(lev, " [observed]"), lev)
long$name <- factor(long$name, levels(long$name)[order(lev %in% names(observables))])

ggplot(long, aes(time, value, group = eps, colour = eps)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  scale_color_dMod_c(name = expression(epsilon), direction = -1, palette = "grey") +
  labs(title = paste0("NFkB, X", direction), x = "time", y = NULL) +
  theme_dMod(base_size = 9)


# ---- every direction, one line each ----------------------------------------

for (i in seq_along(flow)) {
  mv  <- names(res$symmetries[[i]]$generator)
  mvL <- paste0(mv, "_l10")
  # a polynomial direction can escape to infinity at finite eps -- a property of the
  # direction, not a failure. Lower eps and click through that one on its own.
  z   <- try(flow[[i]](eps, c(truthL[mvL], truth[setdiff(coords, mv)]),
                       deriv = FALSE)$flow, silent = TRUE)
  if (inherits(z, "try-error")) {
    cat(sprintf("X%d  %-8s escapes before eps = %.2f\n", i,
                res$symmetries[[i]]$type, max(eps)))
    next
  }
  d <- sapply(seq_len(nrow(z)), function(j) {
    q <- as.data.frame(prd(times, replace(truthL, mvL, z[j, mvL]), deriv = FALSE))
    max(abs(q$value[obs] - pred[[1]]$value[obs]) / rms)
  })
  cat(sprintf("X%d  %-8s on %-40s  eps %.2f  moved %.2f  deviates %.2e\n", i,
              res$symmetries[[i]]$type,
              paste(res$symmetries[[i]]$support, collapse = ","),
              max(z[, "time"]),
              max(abs(10^(z[nrow(z), mvL] - truthL[mvL]) - 1)), max(d)))
}


# ============================================================================
# EGF cascade at steady state, with a dose
# ============================================================================
# The coordinates are no longer the model's own parameters: `trafo` holds the
# network at rest, so the twin is g * x * p with that same transformation. A state
# symbol on the right of a trafo entry is an INITIAL VALUE, not the running state --
# reading it as the state is what used to collapse f to zero and report the model as
# far more identifiable than it is.

egf <- eqnlist() |>
  addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR") |>
  addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR") |>
  addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK / (Km_MEK * (1 + k_inh_ppERK * ppERK) + MEK)") |>
  addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK / (Km_pMEK + pMEK)") |>
  addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK / (Km_ERK + ERK)") |>
  addReaction("pERK", "ERK", "k_dephos_pERK * pERK / (Km_dpERK + ppERK)") |> 
  addReaction("pERK", "ppERK", "k_phos_pERK * pMEK * pERK / (Km_pERK + pERK)") |>
  addReaction("ppERK", "pERK", "k_dephos_ppERK * ppERK / (Km_dppERK + ppERK)")
observables2 <- eqnvec(tMEK_obs = "MEK + pMEK", ppERK_obs = "ppERK")

# steadyStates(egf) with the four moieties as custom totals, transcribed so the
# script does not depend on the representative sympy's solve happens to return
ss <- steadyStates(egf) |> as.eqnvec()

dose <- eventlist(var = "EGF", time = 0, value = "1", method = "add")

res2 <- symmetryDetection(egf, observables2, trafo = ss, events = dose,
                          method = "observability", reduceCQ = FALSE,
                          reconstruct = TRUE)
summary(res2)

coords2 <- res2$info$coordinates


x2 <- Xf(odemodel(egf, deriv = FALSE, events = dose, modelname = "x2",
                  compile = FALSE), condition = "C1",
         optionsOde = list(atol = 1e-11, rtol = 1e-11, maxsteps = 1e7))
g2 <- Y(observables2, f = x2, attach.input = TRUE, modelname = "g2", compile = FALSE)

# the steady state first, then the same exp10 step: the outer coordinates of p2 are
# the trafo's right-hand-side symbols, on a log10 scale
p2 <- eqnvec() |>
  define("x~x", x = getParameters(g2, x2)) |>
  insert("x~y", x = names(ss), y = ss) |>
  insert("x ~ exp10(x_l10)", x = .currentSymbols) |>
  P(modelname = "p2", compile = FALSE)

flow2 <- lapply(seq_along(res2$symmetries), function(i)
  Xf(odemodel(log10Transform(res2$symmetries[[i]]$generator), deriv = FALSE,
              modelname = paste0("flow2", i), compile = FALSE), condition = "flow",
     optionsOde = list(atol = 1e-12, rtol = 1e-12, maxsteps = 1e7)))

do.call(compile, c(list(x2, g2, p2), flow2, list(cores = 10)))

prd2 <- g2 * x2 * p2

truth2 <- c(EGFR = 2.1, EGF_EGFR = 0.7, ERK = 1.3, MEK = 0.9, pMEK = 0.4,
            pERK = 0.6, ppERK = 0.35,
            k_bind = 0.8, k_unbind = 0.5, k_dephos_MEK = 0.3, k_phos_pERK = 0.45,
            k_inh_ppERK = 0.6, r_pERK_1 = 1.2,
            Km_MEK = 1.1, Km_pMEK = 0.7, Km_ERK = 1.4, Km_pERK = 0.9,
            Km_dpERK = 0.8, Km_dppERK = 0.65)

truth2L <- setNames(log10(truth2), paste0(names(truth2), "_l10"))

times2 <- seq(0, 50, len = 300)
# the coordinates stay positive well past this, but EGF + EGF_EGFR is conserved and
# further out the flow pushes EGF itself negative -- legal for the ODE, distracting
# in a figure
eps2   <- seq(0, 0.1103, len = 100)

mv2  <- names(res2$symmetries[[1]]$generator)
mv2L <- paste0(mv2, "_l10")
z2   <- flow2[[1]](eps2, c(truth2L[mv2L], truth2[setdiff(coords2, mv2)]),
                   deriv = FALSE)$flow

plot(z2)

pred2 <- lapply(seq_len(nrow(z2)), function(j)
  as.data.frame(prd2(times2, replace(truth2L, mv2L, z2[j, mv2L]), deriv = FALSE)))

obs2 <- pred2[[1]]$name %in% names(observables2)
rms2 <- ave(pred2[[1]]$value[obs2], pred2[[1]]$name[obs2],
            FUN = function(v) sqrt(mean(v^2)))

sapply(pred2, function(d) max(abs(d$value[obs2] - pred2[[1]]$value[obs2]) / rms2))


long2 <- do.call(rbind, Map(function(d, e) transform(d, eps = e), pred2, z2[, "time"]))
lev2 <- levels(long2$name)
levels(long2$name) <- ifelse(lev2 %in% names(observables2),
                             paste0(lev2, " [observed]"), lev2)
long2$name <- factor(long2$name,
                     levels(long2$name)[order(lev2 %in% names(observables2))])

ggplot(long2, aes(time, value, group = eps, colour = eps)) +
  geom_line() +
  facet_wrap(~ name, scales = "free_y") +
  scale_color_dMod_c(name = expression(epsilon), direction = -1, palette = "grey") +
  labs(title = "EGF at steady state with a dose, X1", x = "time", y = NULL) +
  theme_dMod(base_size = 9)
