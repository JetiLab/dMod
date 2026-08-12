## Backend comparison for steadyStates() on three models: what 1.2 fixes
## over 1.1, and what 1.3 fixes over 1.2. Ranking is 1.3 >= 1.2 > 1.1.
## Standalone script, not part of the man page.

library(dMod)
setwd(tempdir())

## Largest ODE residual at a random positive point.
residual <- function(reactions, ss, seed = 1) {
  odes <- as.eqnvec(reactions)
  defs <- ss[trimws(ss) != names(ss)]
  symbols <- unique(unlist(lapply(c(as.character(odes), ss),
                                  function(e) all.vars(parse(text = e)))))
  set.seed(seed)
  free <- setdiff(symbols, names(defs))
  env <- setNames(as.list(runif(length(free), 0.3, 2)), free)
  for (nm in names(defs)) env[[nm]] <- eval(parse(text = defs[[nm]]), env)
  max(abs(vapply(odes, function(e) eval(parse(text = e), env), numeric(1))))
}

## ------ 1.1 -> 1.2: hidden sinks --------------------------------------
## TGFb + C leak jointly through k_dg_C. No single ODE shows it, so both
## being zero is only visible to the sink-cluster search (v1.2+).

m1 <- eqnlist()
m1 <- addReaction(m1, "", "R", "k_pr_R", "synthesis")
m1 <- addReaction(m1, "R", "", "k_dg_R*R", "degradation")
m1 <- addReaction(m1, "TGFb + R", "C", "k_on*TGFb*R", "binding")
m1 <- addReaction(m1, "C", "TGFb + R", "k_off*C", "dissociation")
m1 <- addReaction(m1, "C", "", "k_dg_C*C", "complex degradation")

## v1.1 misses the sink: TGFb stays free and k_off is spent on a
## division by zero, yet the version's own test claims success.
steadyStates(m1, version = "1.1")

## v1.2+ report 'These states are zero: TGFb, C'.
ss1 <- steadyStates(m1)
stopifnot(ss1[["TGFb"]] == "0", ss1[["C"]] == "0",
          residual(m1, ss1) < 1e-10)

## ------ 1.1 -> 1.2 -> 1.3: positivity ---------------------------------
## Basal plus induced synthesis and a saturating sink. A steady-state
## trafo must stay non-negative for every parameter choice, otherwise
## optimization can start states at negative values.

m2 <- eqnlist()
m2 <- addReaction(m2, "", "Src", "k_pr_Src", "source")
m2 <- addReaction(m2, "Src", "Lig", "k_sec*Src", "secretion")
m2 <- addReaction(m2, "", "Rec", "(k_pr_Rec + k_basal_Rec)", "two-arm synthesis")
m2 <- addReaction(m2, "Rec", "", "k_dg_Rec*Rec", "degradation")
m2 <- addReaction(m2, "Lig + Rec", "Cpx", "k_form*Lig*Rec/(Km + Rec)", "binding")
m2 <- addReaction(m2, "Cpx", "", "k_dg_Cpx*Cpx", "complex degradation")

## v1.1 returns k_dg_Rec as a difference, negative for some parameters.
ss2_11 <- steadyStates(m2, version = "1.1")
stopifnot(any(grepl("- ", ss2_11, fixed = TRUE)))

## v1.2 has the positivity contract but substitutes solutions eagerly.
## By the time Rec is reached its sink is a constant, no pivot is
## sign-definite, and it refuses with a diagnosis.
stopifnot(identical(steadyStates(m2, version = "1.2"), 0))

## v1.3 keeps fluxes symbolic. The two synthesis arms are still intact
## and are spent as a ratio pivot: no subtraction anywhere.
ss2 <- steadyStates(m2)
stopifnot(!any(grepl("-", ss2, fixed = TRUE)), residual(m2, ss2) < 1e-10)

## ------ 1.2 -> 1.3: runtime --------------------------------------------
## TGFb-like network, 20 states: receptor turnover, complex cycle, a
## shared saturable Smad7 pool, Smad trimer with conserved totals, and
## downstream mRNA/protein chains. v1.2 substitutes every solution into
## all fluxes immediately, so expressions explode and cancel()/GCD
## dominates. v1.3 records solutions and resolves once at output time.

m3 <- eqnlist()
m3 <- addReaction(m3, "", "TGFbmRNA", "k_pr_TGFbmRNA_basal", "basal transcription")
m3 <- addReaction(m3, "", "TGFbmRNA", "k_pr_TGFbmRNA * C3", "autoinduction")
m3 <- addReaction(m3, "TGFbmRNA", "", "k_dg_TGFbmRNA * TGFbmRNA", "decay")
m3 <- addReaction(m3, "", "TGFb", "k_pr_TGFb * TGFbmRNA", "secretion")
m3 <- addReaction(m3, "TGFb", "", "k_depl_TGFb * TGFb", "depletion")
m3 <- addReaction(m3, "", "R1mRNA", "k_pr_R1mRNA * (1 + k_act_FB * FB)", "feedback transcription")
m3 <- addReaction(m3, "R1mRNA", "", "k_dg_R1mRNA * R1mRNA", "decay")
m3 <- addReaction(m3, "", "R1", "k_pr_R1 * R1mRNA", "translation")
m3 <- addReaction(m3, "R1", "", "k_dg_R1 * R1", "turnover")
m3 <- addReaction(m3, "", "R2", "k_pr_R2", "synthesis")
m3 <- addReaction(m3, "R2", "", "k_dg_R2 * R2", "turnover")
m3 <- addReaction(m3, "TGFb + R2", "R2_TGFb", "k_form_R2_TGFb * TGFb * R2", "capture")
m3 <- addReaction(m3, "R2_TGFb", "", "k_dg_R2_TGFb * R2_TGFb", "turnover")
m3 <- addReaction(m3, "R1 + R2_TGFb", "C", "k_form_C * R1 * R2_TGFb", "signaling complex")
m3 <- addReaction(m3, "C", "C_int", "k_int * C", "internalization")
m3 <- addReaction(m3, "C_int", "C", "k_rec * C_int", "recycling")
m3 <- addReaction(m3, "C", "", "k_dg_C * C", "turnover")
m3 <- addReaction(m3, "R1 + S7", "", "k_dg_R1_S7 * R1 * S7 / (Km + R1 + C_int)", "shared S7 pool")
m3 <- addReaction(m3, "C_int + S7", "", "k_dg_C_S7 * C_int * S7 / (Km + R1 + C_int)", "shared S7 pool")
m3 <- addReaction(m3, "Smad2", "pSmad2", "k_p2 * C_int * Smad2 / (1 + k_inh * S7)", "phosphorylation")
m3 <- addReaction(m3, "pSmad2", "Smad2", "k_dp2 * pSmad2", "dephosphorylation")
m3 <- addReaction(m3, "Smad3", "pSmad3", "k_p3 * C_int * Smad3 / (1 + k_inh * S7)", "phosphorylation")
m3 <- addReaction(m3, "pSmad3", "Smad3", "k_dp3 * pSmad3", "dephosphorylation")
m3 <- addReaction(m3, "pSmad2 + pSmad3 + Smad4", "C3", "k_form_C3 * pSmad2 * pSmad3 * Smad4", "trimer")
m3 <- addReaction(m3, "C3", "pSmad2 + pSmad3 + Smad4", "k_dec_C3 * C3", "trimer decay")
m3 <- addReaction(m3, "", "S7mRNA", "k_pr_S7mRNA_basal", "basal transcription")
m3 <- addReaction(m3, "", "S7mRNA", "k_pr_S7mRNA * C3", "induced transcription")
m3 <- addReaction(m3, "S7mRNA", "", "k_dg_S7mRNA * S7mRNA", "decay")
m3 <- addReaction(m3, "", "S7", "k_pr_S7 * S7mRNA", "translation")
m3 <- addReaction(m3, "S7", "", "k_dg_S7 * S7", "turnover")
m3 <- addReaction(m3, "", "FBmRNA", "k_pr_FBmRNA * C3", "transcription")
m3 <- addReaction(m3, "FBmRNA", "", "k_dg_FBmRNA * FBmRNA", "decay")
m3 <- addReaction(m3, "", "FB", "k_pr_FB * FBmRNA", "translation")
m3 <- addReaction(m3, "FB", "", "k_dg_FB * FB", "turnover")
m3 <- addReaction(m3, "", "ROmRNA", "k_pr_ROmRNA_basal", "basal transcription")
m3 <- addReaction(m3, "", "ROmRNA", "k_pr_ROmRNA * C3", "induced transcription")
m3 <- addReaction(m3, "ROmRNA", "", "k_dg_ROmRNA * ROmRNA", "decay")
m3 <- addReaction(m3, "", "RO", "k_pr_RO * ROmRNA", "translation")
m3 <- addReaction(m3, "RO", "", "k_dg_RO * RO", "turnover")
m3 <- customTotals(m3, list(totalS2 = "Smad2 + pSmad2 + C3",
                            totalS3 = "Smad3 + pSmad3 + C3",
                            totalS4 = "Smad4 + C3"))

t12 <- system.time(ss3_12 <- steadyStates(m3, version = "1.2"))["elapsed"]
t13 <- system.time(ss3 <- steadyStates(m3))["elapsed"]
## one measurement: v1.2 218 s, v1.3 2.6 s
print(c(v1.2 = t12, v1.3 = t13))
stopifnot(residual(m3, ss3) < 1e-10)
