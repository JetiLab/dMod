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


optionsOde = list(atol = 1e-4, rtol = 1e-4)
optionsSens = list(atol = 1e-4, rtol = 1e-4)
optionsSens2 = list(atol = 1e-4, rtol = 1e-4)

x.ds <- odemodel(reactions, events = events, compile = F, deriv2 = F, solver = "deSolve", modelname = "AB_deSolve") %>% 
  Xs(condition = "cond1", 
     optionsOde = optionsOde,
     optionsSens = optionsSens)

x.bt <- odemodel(reactions, events = events, compile = F, deriv2 = T, solver = "boost", modelname = "AB_boost") %>% 
  Xs(condition = "cond1",
     optionsOde = optionsOde,
     optionsSens = optionsSens,
     optionsSens2 = optionsSens2)



times <- seq(-10, 100, len = 300)
# debugonce(x)
pars <- c(A=4.5, B=0.75, k_p = 0.3, k_d = 0.4, k1=0.1, k2=0.2, dose = 1)
# out.ds <- x.ds(times, pars)
# out.bt <- x.bt(times, pars, deriv2 = F)
# 
# outdsframe <- out.ds %>% as.data.frame()
# outbtframe <- out.bt %>% as.data.frame()
# 
# plot(out.ds)
# plot(out.bt)
# getDerivs(out.ds) %>% plot()
# getDerivs(out.bt) %>% plot()
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

p.expl(pouter) %>% getDerivs()
p.impl(pouter) %>% getDerivs()

compile(x.ds, x.bt, p.impl, p.expl, output = "dMod2test", cores = 8)

prd.impl.ds <- x.ds*p.impl
prd.impl.bt <- x.bt*p.impl
prd.expl.ds <- x.ds*p.expl
prd.expl.bt <- x.bt*p.expl
# debugonce(prd)
out.impl.ds <- prd.impl.ds(times, pouter)
out.expl.ds <- prd.expl.ds(times,pouter)
out.impl.bt <- prd.impl.bt(times, pouter, deriv2 = T)
out.expl.bt <- prd.expl.bt(times,pouter, deriv2 = T)
plot(out.expl.ds)


plot(out.impl.bt)

derivs.bt.impl <- getDerivs(out.expl.bt) %>% as.data.frame() %>% mutate(parfun = "impl")
derivs.bt.expl <- getDerivs(out.impl.bt) %>% as.data.frame() %>% mutate(parfun = "expl")

derivs.bt <- rbind(derivs.bt.expl, derivs.bt.impl)

derivs2.bt.impl <- getDerivs2(out.expl.bt, full = T) %>% as.data.frame() %>% 
  mutate(parfun = "impl")
derivs2.bt.expl <- getDerivs2(out.impl.bt, full = T) %>% as.data.frame() %>% 
  mutate(parfun = "expl")

derivs2.bt <- rbind(derivs2.bt.expl, derivs2.bt.impl)

ggplot(derivs.bt, aes(x=time, y=value, linetype = parfun, color = parfun)) +
  geom_line() + 
  facet_wrap(~name, scales = "free") + 
  theme_dMod() +
  scale_color_dMod()

ggplot(derivs2.bt, aes(x=time, y=value, linetype = parfun, color = parfun)) +
  geom_line() + 
  facet_wrap(~name, scales = "free") + 
  theme_dMod() +
  scale_color_dMod()

system.time({prd.impl.bt(times, pouter, deriv2 = T)})
system.time({prd.expl.bt(times, pouter, deriv2 = T)})


reactions2 <- eqnvec() %>% 
  addReaction("ERK", "ppERK", "k1 * ERK") %>% 
  addReaction("ppERK", "ERK", "k2 * ppERK")

reactions2
p.impl2 <- Pimpl(reactions2, attach.input = T, deriv2 = T, condition = "cond1", compile = F)
