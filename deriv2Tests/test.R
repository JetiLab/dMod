rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(5555)

library(dMod)
library(dplyr)

reactions <- eqnvec() %>% 
  addReaction("A", "B", "k1 * A") %>% 
  addReaction("B", "A", "k2 * B")


events <- eventlist() %>% 
  addEvent(var = "A", time = 0, value = "dose", method = "replace")

x <- odemodel(reactions, events = events, compile = T, deriv2 = T, solver = "boost") %>% Xs(condition = "cond1")
# loadDLL(x)
times <- seq(-10, 100, len = 400)
# debugonce(x)
out <- x(times, c(A=0, B=0, k1=0.1, k2=0.2, dose = 1), deriv = T)

plot(out)

innerpars <- getParameters(x)
trafo <- eqnvec() %>% 
  define("x~x", x = innerpars) %>% 
  insert("x~exp(log(10)*y)", x = setdiff(innerpars, "dose"), y = toupper(setdiff(innerpars, "dose")))

p <- P(trafo, compile = F, condition = "cond1")

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
pouter["dose"] <- 1

pinner <- p(pouter)

prd <- x*p
out <- prd(times, pouter, deriv2 = T)
plot(out)

# derivs <- getDerivs(out)
# plot(derivs)
# derivs2 <- getDerivs2(out, full = F)
# plot(derivs2)

f <- c(A = "-k1*A + k2*B",
       B = "k1*A - k2*B")
P.steadyState <- Pimpl(f, "A")

p.outerValues <- c(k1 = 1, k2 = 0.1, A = 10, B = 1)
P.steadyState(p.outerValues)


