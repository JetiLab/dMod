rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(dplyr)

reactions <- eqnvec() %>% 
  addReaction("", "A", "k_p") %>% 
  addReaction("A", "B", "k1 * A") %>% 
  addReaction("B", "A", "k2 * B") %>% 
  addReaction("B", "", "k_d * B")


events <- eventlist() %>%
  addEvent(var = "A", time = 0, value = "dose", method = "add")


optionsOde = list(atol = 1e-4, rtol = 1e-4)
optionsSens = list(atol = 1e-4, rtol = 1e-4)

x.ds <- odemodel(reactions, events = events, compile = F, deriv2 = F, solver = "deSolve", modelname = "AB_deSolve") %>% 
  Xs(condition = "cond1", 
     optionsOde = optionsOde,
     optionsSens = optionsSens)

x.bt <- odemodel(reactions, events = events, compile = F, deriv2 = T, solver = "boost", modelname = "AB_boost") %>% 
  Xs(condition = "cond1",
     optionsOde = optionsOde,
     optionsSens = optionsSens)

compile(x.ds, x.bt, cores = 8)
times <- seq(-10, 100, len = 300)
# debugonce(x)
out.ds <- x.ds(times, c(A=4.5, B=0.75, k_p = 0.3, k_d = 0.4, k1=0.1, k2=0.2, dose = 1))
out.bt <- x.bt(times, c(A=4.5, B=0.75, k_p = 0.3, k_d = 0.4, k1=0.1, k2=0.2, dose = 1), deriv2 = T)

plot(out.ds)
plot(out.bt)
getDerivs(out.ds) %>% plot()
getDerivs(out.bt) %>% plot()
# getDerivs2(out) %>% plot()

# innerpars <- getParameters(x)
# trafo <- eqnvec() %>% 
#   define("x~x", x = innerpars) %>% 
#   insert("x~exp(log(10)*y)", x = setdiff(innerpars, "dose"), y = toupper(setdiff(innerpars, "dose")))
# 
# p <- P(trafo, compile = F, condition = "cond1")
# 
# outerpars <- getParameters(p)
# pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
# pouter["dose"] <- 1
# 
# pinner <- p(pouter)
# 
# prd <- x*p
# out <- prd(times, pouter, deriv2 = T)
# plot(out)

# derivs <- getDerivs(out)
# plot(derivs)
# derivs2 <- getDerivs2(out, full = F)
# plot(derivs2)

reactions

p.impl <- Pimpl(reactions, attach.input = T, deriv2 = T, condition = "cond1", compile = F)
p.impl
pouter <- c(k_p = 0.3, k1 = 0.1, k2 = 0.2, k_d = 0.4, dose = 1)

p.impl(pouter, deriv2 = T)

innerpars <- getParameters(x.ds)
ss.expl <- c(A="k_p * (k2 + k_d) / (k1*k_d)", B = "k_p/k_d")

trafo.expl <- eqnvec() %>% 
  define("x~x", x = innerpars) %>% 
  insert("x~y", x = names(ss.expl), y = ss.expl)

p.expl <- P(trafo.expl, condition = "cond1", deriv2 = T, compile = F)


prd.impl <- x.bt*p.impl
prd.expl.ds <- x.ds*p.expl
prd.expl.bt <- x.bt*p.expl
# debugonce(prd)
out.impl <- prd.impl(times, pouter, deriv2 = F)
out.expl.ds <- prd.expl.ds(times,pouter)
out.expl.bt <- prd.expl.bt(times,pouter, deriv2 = T)
plot(out.expl.ds)
plot(out.expl.bt)
getDerivs(out.expl.ds) %>% plot()
getDerivs(out.expl.bt) %>% plot()
getDerivs(out.impl) %>% plot()

compile(p.impl, p.expl)

system.time({prd.impl(times,pouter, deriv2 = T)})