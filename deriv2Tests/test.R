rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(dplyr)
library(ggplot2)

reactions <- eqnvec() %>% 
  addReaction("", "A", "k_p") %>% 
  addReaction("A", "B", "k1 * A") %>% 
  addReaction("B", "A", "k2 * B") %>% 
  addReaction("B", "", "k_d * B")


events <- eventlist() %>%
  addEvent(var = "A", time = 0, value = "dose", method = "add")


myopts = list(atol = 1e-4, rtol = 1e-4)

x <- odemodel(reactions, events = events, compile = F, deriv2 = T, solver = "boost", modelname = "toymodel") %>% 
  Xs(condition = "cond1",
     optionsOde = myopts,
     optionsSens = myopts,
     optionsSens2 = myopts)

p.SS <- P(reactions, attach.input = T, deriv2 = T, compile = F)
pouter <- c(k_p = 0.3, k1 = 0.1, k2 = 0.2, k_d = 0.4, dose = 1)
p.SS(pouter, deriv2 = T)

p.SS(pouter, deriv2 = T) %>% getDerivs()


observables <- eqnvec(
  B_obs =  "scale_B * B + offset_B"
)
g <- Y(observables, f = x, compile = F, modelname = "toymodel_obsfun", attach.input = F)
innerpars <- setdiff(getParameters(x, p.SS, g), reactions$states)
trafo <- eqnvec() %>% 
  define("x~x", x = innerpars) %>% 
  insert("x~exp(log(10)*y)", x = setdiff(.currentSymbols, "dose"), y = toupper(setdiff(.currentSymbols, "dose")))
p.logtrafo <- P(trafo, deriv2 = T, condition = "cond1", compile = F)

getParameters(p.logtrafo)
pouter <- structure(rep(-1, length(getParameters(p.logtrafo))), names = getParameters(p.logtrafo))
pouter["dose"] <- 1
poutLog <- p.logtrafo(pouter)
poutLog
poutLog %>% getDerivs()
pout <- p.SS(poutLog$cond1)
pout %>% getDerivs()

p <- p.SS*p.logtrafo

times <- seq(-10, 100, len = 300)

p(pouter)
p(pouter) %>% getDerivs()


compile(g,x,p, output = "toymodelAll", cores = 8)


prd <- g*x*p
# debugonce(g)
out <- prd(times, pouter, deriv2 = T)
plot(out)
out %>% getDerivs() %>% plot()
out %>% getDerivs2() %>% plot()
system.time({prd(times, pouter, deriv2 = T)})
system.time({prd(times, pouter, deriv2 = F)})
