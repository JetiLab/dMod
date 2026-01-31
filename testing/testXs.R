rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
unlink(list.files(".", pattern = "\\.(cpp|c|o|so|dll)$", full.names = TRUE), force = TRUE)

set.seed(5555)

library(dMod)
library(ggplot2)
library(dplyr)


# Set up reactions
eqns <- eqnvec(x = "-k*x")

x <- odemodel(eqns, modelname = "test_Xs") %>% Xs()

p <- eqnvec() %>% 
  define("x~x", x = getParameters(x)) %>% 
  insert("x~1") %>% P(compile = T)
