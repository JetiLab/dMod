## Load Libraries -----
rm(list=ls())
library(dMod)
library(dplyr)
library(ggplot2)
theme_set(theme_dMod())
setwd(tempdir()) # change if you want to keep the files that dMod creates internally
set.seed(5555)

## Set up the ODE model -----------------------------------------------------------------------

# Define via reactions

reactions <- eqnlist() %>%
  addReaction("TCA_buffer", "TCA_cell",  rate = "k_import*TCA_buffer", description = "Uptake") %>%
  addReaction("TCA_cell", "TCA_buffer",  rate = "k_export_sinus*TCA_cell", description = "Sinusoidal export") %>%
  addReaction("TCA_cell", "TCA_cana",    rate = "k_export_cana*TCA_cell", description = "Canalicular export") %>%
  addReaction("TCA_cana", "TCA_buffer",  rate = "k_reflux*TCA_cana", description = "Reflux into the buffer")

# Define via ODE system (equivalent)

# f <- eqnvec(TCA_buffer = "-k_import * TCA_buffer + k_export_sinus * TCA_cell + k_reflux * TCA_cana",
#             TCA_cana = "k_export_cana * TCA_cell - k_reflux * TCA_cana",
#             TCA_cell = "k_import * TCA_buffer - k_export_sinus * TCA_cell - k_export_cana * TCA_cell")

# Translate reactions into ODE model object
mymodel <- odemodel(reactions, modelname = "bamodel", compile = F, solver = "boost", deriv2 = T)
# Generate trajectories for the default condition
x <- Xs(mymodel)

# Define observables buffer and cellular
observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
g <- Y(observables, f = x, condition = NULL, compile = F, modelname = "obsfn_bamodel", attach.input = F, deriv2 = T)

# Define parameter transformations using define(), insert() and branch(). Old function repar also avaiable!
innerpars <- getParameters(x,g)
trafo <- NULL %>%
  define("x~x", x = innerpars) %>% # identity
  define("TCA_buffer~0") %>%
  insert("x~10^y", x = .currentSymbols, y = toupper(.currentSymbols))


# # # Explicit trafo (this is equivalent to the lines above)
# trafo <- eqnvec(TCA_buffer = "0",
#                 TCA_cell = "exp(log(10)*TCA_cell)",
#                 TCA_cana = "exp(log(10)*TCA_cana)",
#                 k_import = "exp(log(10)*k_import)",
#                 k_export_sinus = "exp(log(10)*k_export_sinus)",
#                 k_export_cana = "exp(log(10)*k_export_cana)",
#                 k_reflux = "exp(log(10)*k_reflux)",
#                 s = "exp(log(10)*s)")

p <- P(trafo, condition = "closed", compile = F, deriv2 = T)


### Compile the objects
compile(g, x, p, cores = 5) # Compile C/C++ output of odemodel in parallel

## Use simulate data to calibrate outer model parameters -------------------------------------------------------------
data(badata)
data <- badata %>% subset(condition == "closed") %>% as.datalist()
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

p(pouter[!grepl("K_IMPORT", outerpars)], fixed = c(K_IMPORT = -2))


times <- seq(0, 45, len = 300)

p(pouter)
(x*p)(times, pouter)

# # One Fit
obj <- normL2(data, g * x * p)
obj(pouter)

# Fit on time (starting from pouter)
myfit <- trust(obj, pouter, rinit = 0.1, rmax = 5, iterlim = 500, printIter = T)
mypred <- (g * x * p)(times, myfit$argument)
plot(mypred, data)

obj(myfit$argument)

system.time({obj(myfit$argument, deriv2 = F)})
system.time({obj(myfit$argument, deriv2 = T)})
## Handling different experimental conditions

# Define reflux and open condition according to "Dynamic Modelling in R p. 19-20"
pars["k_reflux"] <- 1e3
outD_open <- (g * x)(timesD, pars, conditions = "open")
datasheet <- subset(as.data.frame(outD_open), time %in% timesD & name %in% names(observables))
datasheet <- within(datasheet, {
  sigma <- sqrt(value + 1)
  value <- rnorm(length(value), value, sigma)
})
data <- data + as.datalist(datasheet)
plotData(data)

# Parameter Trafo, usage of "+" operator for trafo functions (output of P())
trafo <- getEquations(p, conditions = "closed") %>% 
  insert("K_REFLUX~K_REFLUX_OPEN")
p <- p + P(trafo, condition = "open", compile = T, deriv2 = T)

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*p)(times, pouter),data)

datasheet <- as.data.frame(data)
p(pouter, deriv2 = T)

# Objective function
obj <- normL2(data, g * x * p) + constraintL2(pouter, sigma = 4)

# Evaluation of obj at pouter
obj(pouter, deriv2 = T)

system.time({obj(pouter, deriv2 = F)})
system.time({obj(pouter, deriv2 = T)})

myfit <- trust(obj, pouter, rinit = 0.1, rmax = 5, iterlim = 500, printIter = T, deriv2 = F)

# Fit 50 times, sample with sd=4 around pouter
out_frame <- mstrust(obj, pouter, sd = 4, studyname = "bamodel", cores=8, fits=50, iterlim = 200, deriv2 = F)
out_frame <- as.parframe(out_frame)
plotValues(out_frame) # Show "Waterfall" plot
# plotPars(out_frame) # Show parameter plot
bestfit <- as.parvec(myfit$argument)


# Plot predictions along data
plot((g * x * p)(times, bestfit), data)
# 
# # Plot sensis
plot(getDerivs((g * x * p)(times, bestfit)))
# 
# Calculate Parameter Profiles and plot different contributions (for identifiablility only "data" is of interest)
profiles <- profile(obj, bestfit, whichPar = c("K_REFLUX"), method = "optimize", cores = 1, stepControl = list(stop = "data"))
# plotProfile(profiles)
plotProfile(profiles, mode %in% c("data", "prior"))
# plotPaths(profiles, whichPar = "TCA_cana")
# plotPaths(profiles, whichPar = "export_cana")
# plotPaths(profiles, whichPar = "s")
# 
# # Tighten model assumptions with steady state constraint
# pSS <- NULL
# equations <- getEquations(p) # getEquations retrieves the "trafo" from p 
# conditions <- names(equations)
# 
# # Set up steady state constraint
# for (n in conditions) {
#   equations[[n]]["TCA_cana"] <- "exp(log(10)*export_cana)*exp(log(10)*TCA_cell)/exp(log(10)*reflux)"
#   pSS <- pSS + P(equations[[n]], condition = n)
# }
# outerpars <- getParameters(pSS)
# pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
# plot((g*x*pSS)(times, pouter),data)
# 
# # Objective function
# obj <- normL2(data, g * x * pSS, attr.name = "data") + constraintL2(pouter, sigma = 5, attr.name = "prior")
# 
# # Fit 50 times again
# out_frame <- mstrust(obj,pouter,studyname = "bamodel", cores=10,fits=50)
# out_frame <- as.parframe(out_frame)
# plotValues(out_frame) # Show "Waterfall" plot
# plotPars(out_frame) # Show parameter plot
# bestfit <- as.parvec(out_frame)
# plot((g * x * pSS)(times, bestfit), data)
# 
# 
# # Calculate Parameter Profiles
# profiles <- profile(obj, bestfit, names(bestfit), cores = 10, limits = c(lower = -10, upper = 10),
#                     stepControl = list(stepsize = 1e-4, min = 1e-4, max = Inf, atol = 1e-2, rtol = 1e-2, limit = 100, stop = "data"))
# plotProfile(profiles,mode %in% c("data", "prior"))
# plotPaths(profiles, whichPar = "s")
# 
# 
# # s and initial value TCA_cell render model structurally non identifiable
# #  --> fix s to 1 (on log scale to 0)
# 
# trafo <- getEquations(pSS)
# trafo <- repar("x~0", x = "s", trafo)
# p <- P(trafo)
# 
# 
# # Objective function
# outerpars <- getParameters(p)
# pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
# obj <- normL2(data, g * x * p, attr.name = "data") + constraintL2(pouter, sigma = 20, attr.name = "prior")
# 
# # Fit 50 times again
# out_frame <- mstrust(obj,pouter,studyname = "bamodel", cores=10,fits=50)
# out_frame <- as.parframe(out_frame)
# plotValues(out_frame) # Show "Waterfall" plot
# plotPars(out_frame) # Show parameter plot
# bestfit <- as.parvec(out_frame)
# plot((g * x * p)(times, bestfit), data)
# 
# # remove the prior
# obj <- normL2(data, g * x * p, attr.name = "data")
# 
# # refit
# refit <- trust(obj, bestfit, rinit = 1, rmax = 10)
# plot((g * x * p)(times, refit$argument), data)
# 
# profiles <- profile(obj, refit$argument, names(refit$argument), cores = 10, limits = c(lower = -5, upper = 5), method = "optimize",
#                     stepControl = list(stop = "data"))
# plotProfile(profiles, mode == "data")
# plotPaths(profiles, whichPar = "reflux_open")
# 
# # Note: The parameter "reflux_open" is still practical non-identifiable 
# # The profile is open to the left -> A possible reduction of the model would be the assumption of an immediate reflux, i.e. the limit case reflux_open -> infinity
# # An pragmatic but ugly way to circumvent the reformulation of the ODE system is to fix the reflux_open parameter to a high value. E.g. 1e3
# 
# # One could also check the models ability to produce reliable predictions, by the calculation of prediction uncertainty with profile likelihood
# # The calculation of prediction confidence intervals is done at next
# 
# ## Prediction uncertainty taken from validation profile --------------------------------------------------------------------------
# 
# # choose sigma below 1 percent of the prediction in order to pull the prediction strongly towards d1
# obj.validation <- normL2(data, g * x * p, times = c(10), attr.name = "data") +
#   datapointL2(name = "TCA_cell", time = 10, value = "v", sigma = 1, attr.name = "validation", condition = "closed")
# 
# # If sigma is not known, and you therefore decide to calculate prediction confidence intervals, just choose a very small sigma, in order to "pull strongly" on the trajectory
# 
# # refit
# myfit <- trust(obj.validation, parinit = c(v = 180, bestfit), rinit = 1, rmax = 10, iterlim = 1000)
# 
# # Calculate profile
# validation_profile <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "optimize",
#                               optControl = list(rinit = .1, rmax = 10, iterlim = 100))
# 
# 
# # plotProfile(validation_profile) # This also plos the prediction colums, which is a bug in the code.
# plotProfile(validation_profile, mode %in% c("validation", "data")) # Plots only the two contributions validation and data, along with the sum (total)
# # Is the contribution of validation small?? # If yes: The total aligns with a prediction profile
# 
# # Confidence Interval of the prediction
# confint(validation_profile, val.column = "value")
# 
# 
# 
# ## Prediction band (prediction uncertainty for several time points) --------------------------------------------------------------
# # Here we calculate a prediction CI for different timepoints. In the end we interpolate to a "prediction band"
# prediction_band <- do.call(rbind, lapply(seq(10, 50, 10), function(t) {
#   
#   cat("Computing prediction profile for t =", t, "\n")
#   
#   obj.validation <- normL2(data, g * x * p, times = c(t), attr.name = "data") + 
#     datapointL2(name = "TCA_cell", time = t, value = "v", sigma = 1, attr.name = "validation", condition = "closed")
#   
#   refit <- trust(obj.validation, parinit = c(v = 180, bestfit), rinit = 1, rmax = 10, iterlim = 1000)
#   
#   profile_prediction <- profile(obj.validation, myfit$argument, "v", cores = 4, method = "optimize")
#   
#   d1 <- confint(profile_prediction, val.column = "value")
#   
#   # Output
#   data.frame(time = t, condition = "closed", name = "TCA_cell",  d1[-1])
#   
# }))
# 
# prediction <- (g * x * p)(times, refit$argument) %>%
#   as.data.frame() 
# 
# 
# prediction_band_spline <- data.frame(
#   time = prediction$time[prediction$time>=10],
#   value = prediction$value[prediction$time>=10],
#   condition = "closed",
#   name = "TCA_cell",
#   lower = spline(prediction_band$time, prediction_band$lower, xout = prediction$time[prediction$time>=10])$y,
#   upper = spline(prediction_band$time, prediction_band$upper, xout = prediction$time[prediction$time>=10])$y
# )
# 
# # Create the ggplot
# ggplot(prediction, aes(x = time, y = value, color = condition)) +
#   geom_line() +  # Line connecting the points for each condition
#   geom_ribbon(data = prediction_band_spline, aes(x = time, ymin = lower, ymax = upper, fill = condition), 
#               lty = 0, alpha = .3, show.legend = F) +  # Show ribbon in the legend
#   geom_point(data = prediction_band, aes(x = time, y = lower, color = condition), shape = 4, show.legend = F) + 
#   geom_point(data = prediction_band, aes(x = time, y = upper, color = condition), shape = 4, show.legend = F) + 
#   facet_wrap(~ name, scales = "free_y") +  # Facet by 'name' column
#   labs(
#     x = "Time",
#     y = "Value",
#     color = "Condition"
#   ) +
#   dMod::theme_dMod() +  # Apply dMod theme
#   dMod::scale_color_dMod() +  # Apply dMod color scale to lines
#   dMod::scale_fill_dMod()   # Apply the same color scale to the fill
# 
# 
# ## Alternative implementation of the Steady state via implicit parameter transformation of the steady state  -----------------------------------------------
# 
# # Redefine reactions in order to control the standard and open condition by events
# reactions <- eqnlist() %>% 
#   addReaction("TCA_buffer", "TCA_cell", rate = "import*TCA_buffer", description = "Uptake")%>% 
#   addReaction("TCA_cell", "TCA_buffer", rate = "export_sinus*TCA_cell", description = "Sinusoidal export")%>% 
#   addReaction("TCA_cell", "TCA_cana", rate = "export_cana*TCA_cell", description = "Canalicular export")%>% 
#   addReaction("TCA_cana", "TCA_buffer", rate = "(reflux*(1-switch) + reflux_open*switch)*TCA_cana", description = "Reflux into the buffer")%>% 
#   addReaction("0", "switch", rate = "0", description = "Create a switch")
# 
# events <- NULL
# events <- addEvent(events, var = "TCA_buffer", time = 0, value = 0)
# events <- addEvent(events, var = "switch" , time = 0, value = "OnOff")
# mymodel <- odemodel(reactions, modelname = "bamodel2", events = events)
# x <- Xs(mymodel)
# 
# 
# 
# # Replace one reaction with a analytical expression for the conserved quantity: TCR_tot
# f <- as.eqnvec(reactions)[c("TCA_buffer", "TCA_cana", "TCA_cell")]
# f["TCA_cell"] <- "TCA_buffer + TCA_cana + TCA_cell - TCA_tot"
# pSS <- P(f, method = "implicit", compile = TRUE, modelname = "pfn") 
# 
# observables <- eqnvec(buffer = "s*TCA_buffer", cellular = "s*(TCA_cana + TCA_cell)")
# 
# innerpars <- unique(c(getParameters(mymodel), getSymbols(observables), getSymbols(f)))
# trafo <- repar("x~x" , x = innerpars)
# trafo <- repar("x~0" , x = reactions$states, trafo)
# 
# trafo <- repar("x~exp(log(10)*x)", x = setdiff(innerpars, "OnOff"), trafo)
# p <- P(repar("OnOff~0", trafo), condition = "closed") + P(repar("OnOff~1", trafo), condition = "open")
# g <- Y(observables, f = x, compile = TRUE, modelname = "obsfn2")
# 
