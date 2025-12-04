rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(ggplot2)
library(dplyr)


# Set up reactions
r <- eqnvec() %>%
  addReaction("", "Stim", "0") %>%
  addReaction("", "R", "k_act_R_bas + k_act_R * Stim", "stim act R") %>%
  addReaction("R", "", "k_deact_R * R", "deact R") %>% 
  addReaction("A", "pA", "k1 * A * R",  "phos A") %>%
  addReaction("pA", "A", "k2 * pA",     "dephos pA")

mysteadies <- steadyStates(r)

p.expl <- eqnvec() %>% 
  define("x~x", x = getParameters(r)) %>% 
  insert("x~y", x = names(mysteadies), y = mysteadies) %>% 
  insert("A~y", y = "totA/(1+k1*k_act_R_bas/(k2*k_deact_R))") %>% 
  P()

p.impl <- P(r, positive = T, nstart = 200,
            optionsRootSolve = list(atol = 1e-10, rtol = 1e-10))

outerpars <- getParameters(p.impl)
pouter <- structure(rep(0.1, length(outerpars)), names = outerpars)
pouter["Stim"] = 0

p.impl(pouter)
p.expl(pouter)
p.impl(pouter) %>% getDerivs()
p.expl(pouter) %>% getDerivs()

getEquations(p.expl)
getEquations(p.impl)