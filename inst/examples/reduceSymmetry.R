\dontrun{
  library(dMod)

  ## A gene cascade observed only through the protein: transcription and
  ## translation cannot be told apart (only ktx*ktl is identifiable)
  f <- eqnvec(m = "ktx - dm*m", p = "ktl*m - dp*p")
  g <- eqnvec(y = "p")
  res <- symmetryDetection(f, g, method = "observability", reconstruct = TRUE)

  ## constructive reduction: invariant monomials, a certified transversal, and
  ## the full admissible family (the gauge choice stays yours)
  red <- reduceSymmetry(res)
  red
  red$trafo                       # ready for P() or symmetryDetection(trafo = )

  ## the emitted trafo removes the direction
  symmetryDetection(f, g, method = "observability", trafo = red$trafo)$identifiable

  ## pre-gauge with your own fixing instead: the reduction then only handles
  ## what your choice leaves over, and reports redundant fixings
  reduceSymmetry(res, fixed = "ktl")

  ## a curved trade-off weighted by a known input: no monomial invariant exists
  ## (certified), the polynomial stage finds k1 + u*k2 and solves it rationally
  f2 <- eqnvec(A = "k1 + u*k2 - kdeg*A")
  res2 <- symmetryDetection(f2, eqnvec(y = "A"), method = "observability",
                            fixed = "u", reconstruct = TRUE)
  red2 <- reduceSymmetry(res2)
  red2$trafo                      # k1 = "k1 - k2*u": the influx lump keeps k1's name

  ## a block beyond the search caps stays honest: negative certificates say what
  ## was excluded instead of silently giving up
  reduceSymmetry(res2, dPoly = 0L, dDarboux = 0L, dExp = 0L)

  ## -- harder curved blocks ---------------------------------------------------

  ## A scaling ENTANGLED with an affine direction: A <-> B observed on a gain.
  ## The two directions share coordinates, so they merge into one curved block;
  ## its three invariants (A*alpha, k1 + k2, alpha*k2*(A + B)) are found across
  ## the monomial and polynomial stages and solved jointly -- pinning either
  ## direction alone would contradict the other.
  f3 <- eqnvec(A = "-k1*A + k2*B", B = "k1*A - k2*B")
  g3 <- eqnvec(y = "alpha*A")
  res3 <- symmetryDetection(f3, g3, method = "observability", reconstruct = TRUE)
  red3 <- reduceSymmetry(res3)
  red3
  symmetryDetection(f3, g3, method = "observability",
                    trafo = red3$trafo)$identifiable

  ## A degree-2 polynomial direction whose invariant is RATIONAL: only the
  ## difference of inverse rates 1/a - 1/b drives the decay, so the generator is
  ## a^2 d/da + b^2 d/db and no monomial or polynomial invariant exists (both
  ## certified negative). The Darboux stage factors the extactic determinant
  ## into a, b, a - b and assembles (a - b)/(a*b) from the cofactor kernel.
  f4 <- eqnvec(x = "-(b - a)/(a*b)*x")
  res4 <- symmetryDetection(f4, eqnvec(y = "x"), method = "observability",
                            reconstruct = TRUE)
  red4 <- reduceSymmetry(res4, dDarboux = 1L)
  red4$blocks[[1]]$invariants     # "a^-1*b^-1*(a - b)"
  symmetryDetection(f4, eqnvec(y = "x"), method = "observability",
                    trafo = red4$trafo)$identifiable

  ## -- the Liouvillian stages ------------------------------------------------

  ## The SAME direction with the Darboux stage capped away: the escalation
  ## continues into exponential factors (Prelle-Singer: every Liouvillian first
  ## integral is a product of Darboux factors and exp(g/h) terms). The stage
  ## certifies exp((b - a)/(a*b)) exactly and the solve becomes a log-solve:
  ## b = -1/(log(b) - 1). Entries with log/exp are fine for P(), but NOT for the
  ## symmetryDetection(trafo = ) round trip, which requires rational entries.
  red4e <- reduceSymmetry(res4, dDarboux = 0L)
  red4e$blocks[[1]]$invariants    # "exp((-a + b)/(a*b))"
  red4e$trafo

  ## A quadratic invariant solved by a root: the rotation generator b d/da -
  ## a d/db leaves a^2 + b^2 invariant; the carrier is solved as a fractional
  ## power (positive branch, noted in the gauge note).
  ## (Directions like this reach reduceSymmetry from symmetryDetection results;
  ## the block then reads:
  ##   invariants: a^2 + b^2      trafo: a = 1, b = sqrt(b - 1))

  ## An integrating factor when no closed-form invariant exists in the rational
  ## or exponential class: the damped-oscillator direction y d/dx - (x + y) d/dy
  ## admits the Darboux-form integrating factor M = 1/(x^2 + x*y + y^2) (the
  ## quadratic IS a Darboux polynomial, cofactor -1, and div X = -1). The first
  ## integral -log(x^2 + x*y + y^2)/2 - sqrt(3)*atan(...)/3 is found by
  ## quadrature and proven exactly, but reported invariantOnly: atan is outside
  ## the trafo language, so no reparametrisation is claimed.

  ## The EGF/EGFR -> MEK/ERK cascade of the vignette, observed only through the
  ## two phospho-forms on unknown gains: four directions. The MEK- and ERK-module
  ## unit scalings reduce as clean scaling blocks; the receptor-module scaling is
  ## ENTANGLED with the genuinely rational confounder and both merge into one
  ## curved block whose four invariants read biologically -- the effective fluxes
  ## EGF_EGFR*k_bind and EGF_EGFR*k_phos_MEK, the relaxation rate
  ## k_unbind + k_bind*(EGF + EGFR), and totalEGF*totalEGFR/EGF_EGFR^2.
  reactions <- eqnlist() |>
    addReaction("EGF + EGFR", "EGF_EGFR", "k_bind * EGF * EGFR")   |>
    addReaction("EGF_EGFR", "EGF + EGFR", "k_unbind * EGF_EGFR")   |>
    addReaction("MEK", "pMEK", "k_phos_MEK * EGF_EGFR * MEK")      |>
    addReaction("pMEK", "MEK", "k_dephos_MEK * pMEK")              |>
    addReaction("ERK", "pERK", "k_phos_ERK * pMEK * ERK")          |>
    addReaction("pERK", "ERK", "k_dephos_ERK * pERK")
  observables <- eqnvec(pMEK_obs = "scale_pMEK * pMEK",
                        pERK_obs = "scale_pERK * pERK")
  egf <- symmetryDetection(reactions, observables, method = "observability",
                           reduceCQ = FALSE, reconstruct = TRUE)
  redEgf <- reduceSymmetry(egf)
  redEgf
  symmetryDetection(reactions, observables, method = "observability",
                    reduceCQ = FALSE, trafo = redEgf$trafo)$identifiable

  ## The Michaelis-Menten depletion assay on an arbitrary readout scale: the
  ## turnover and enzyme amount lump into Vmax = kcat*Etot, and the molar units
  ## of S, Km, Etot trade against the scale s. Both directions land in one
  ## scaling block whose invariants are exactly the classic identifiable set
  ## {S/Km, s*Km, kcat*Etot/Km}.
  f5 <- eqnvec(S = "-kcat*Etot*S/(Km + S)")
  g5 <- eqnvec(y = "s*S")
  res5 <- symmetryDetection(f5, g5, method = "observability", reconstruct = TRUE,
                            reduceCQ = FALSE)
  red5 <- reduceSymmetry(res5)
  red5
  symmetryDetection(f5, g5, method = "observability", reduceCQ = FALSE,
                    trafo = red5$trafo)$identifiable
}
