## Microbenchmark: the regFit / regClustered C++ marginal kernels vs the pure-R
## reference. Run from the repo root: Rscript bench/bench_laplace_kernel.R
## The C++ path is the default (options(dMod.regUseCpp = TRUE)); the R path is the
## bit-comparable reference (agreement asserted in tests/testthat/test-EMlaplace.R).
suppressMessages(devtools::load_all(".", quiet = TRUE))
if (!requireNamespace("microbenchmark", quietly = TRUE))
  stop("install 'microbenchmark' to run this bench.")
library(microbenchmark)

cat("== normalLaplace (K coordinates) ==\n")
set.seed(1); a <- runif(200, 0.5, 50); m <- rnorm(200); lam <- 2
print(microbenchmark(R   = normalLaplace(a, m, lam),
                     Cpp = normalLaplaceCpp(a, m, lam), times = 200L))

cat("\n== clustered marginal (one parameter, G groups; the AGHQ node loop) ==\n")
gh <- .gaussHermite(24L)
for (G in c(4, 5, 6)) {
  set.seed(G); H <- runif(G, 5, 40); mm <- rnorm(G, 0, 0.7)
  fR <- function() { old <- options(dMod.regUseCpp = FALSE); on.exit(options(old))
    .clusterParamMarginal(H, mm, 3, gh, rule = "auto") }
  fC <- function() .clusterParamMarginal(H, mm, 3, gh, rule = "auto")
  cat(sprintf("[ G = %d (anchored d = %d, nq = %d) ]\n", G, G - 1L, .adaptiveNq(G - 1L)))
  print(microbenchmark(R = fR(), Cpp = fC(), times = if (G >= 5) 20L else 50L))
}

cat("\n== grouping score (a full cluster-path candidate marginal) ==\n")
Hd <- rep(30, 6); m6 <- c(0.6, 0.55, 0.5, -0.5, -0.55, -0.6)
P  <- list(1:3, 4:6)
fR <- function() { old <- options(dMod.regUseCpp = FALSE); on.exit(options(old))
  .regClusterGroupScore(P, Hd, m6, gh, rule = "auto") }
fC <- function() .regClusterGroupScore(P, Hd, m6, gh, rule = "auto")
print(microbenchmark(R = fR(), Cpp = fC(), times = 30L))
