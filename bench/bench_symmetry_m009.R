# M009-scale benchmark for symmetryDetection(method = "observability").
#
# The full TGFb 2-receptor model M009 (H1975, capacity-limited Smad7 pool) with
# the resolved steadyStates(simplify = "full") fixed point as per-condition
# trafo: 34 conditions x 3 event segments = 102 tapes whose pivot solves
# (k_pr_S7, k_pr_R1, ...) are very large rational expressions. This is the
# setup that motivated the 2026-08 memo/dedup optimization of the observability
# compile and scaling peel — before it, this call did not finish; after it,
# steadyStates takes ~20s and the SI check ~3-4 min (cores = 10).
#
# The condition grid approximates the H1975 data layout (TC arms with the full
# observable panel, DR dose ladders with the DR panel, ligand-depletion assay
# with the TGFb panel); the verdict on the real data layout may differ.
#
# Usage (from the repo root):
#   Rscript bench/bench_symmetry_m009.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(dMod.sym.verbose = FALSE)
setwd(tempdir())   # steadyStates writes its csv into the working directory

addRC <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "Cell",      ...)
addRE <- function(eq, from, to, rate, ...) addReaction(eq, from, to, rate, compartment = "extraCell", ...)

reactions <- eqnlist() |>
  addRC("", "bool_ActD",  "0") |> addRC("", "bool_CHX",   "0") |>
  addRC("", "bool_MG132", "0") |> addRC("", "bool_R1OE",  "0") |>
  addRC("", "bool_R2OE",  "0") |> addRC("", "bool_R1Knd", "0") |>
  addRC("", "bool_R2Knd", "0") |>
  addRC("", "TGFB1mRNA", "k_pr_TGFB1mRNA_basal * (1 - bool_ActD)") |>
  addRC("", "TGFB1mRNA", "k_pr_TGFB1mRNA * C3 * (1 - bool_ActD)") |>
  addRC("TGFB1mRNA", "", "k_dg_TGFB1mRNA * TGFB1mRNA") |>
  addRE("TGFb", "", "k_depl_TGFb_basal * TGFb") |>
  addRE("", "TGFb", "k_pr_TGFb * TGFB1mRNA * (1 - bool_CHX)", rateCompartment = "Cell") |>
  addRC("", "R1mRNA", "(k_pr_R1mRNA * (1 + k_act_FB3 * FB3) + k_pr_R1mRNA_OE * bool_R1OE) * (1 - bool_ActD)") |>
  addRC("", "R2mRNA", "(k_pr_R2mRNA + k_pr_R2mRNA_OE * bool_R2OE) * (1 - bool_ActD)") |>
  addRC("R1mRNA", "", "(k_dg_R1mRNA + k_si_R1mRNA * bool_R1Knd) * R1mRNA") |>
  addRC("R2mRNA", "", "(k_dg_R2mRNA + k_si_R2mRNA * bool_R2Knd + k_dg_R2mRNA_FB4 * FB4) * R2mRNA") |>
  addRC("", "R1", "k_pr_R1 * R1mRNA * (1 - bool_CHX)") |>
  addRC("", "R2", "k_pr_R2 * R2mRNA * (1 - bool_CHX)") |>
  addRC("R1", "", "k_dg_R1 * R1 * (1 - bool_MG132)") |>
  addRC("R2", "", "k_dg_R2 * R2 * (1 - bool_MG132)") |>
  addRC("R2 + TGFb", "R2_TGFb", "k_form_R2_TGFb * TGFb * R2", rateCompartment = "Cell") |>
  addRC("R2_TGFb", "",          "k_dg_R2_TGFb * R2_TGFb * (1 - bool_MG132)") |>
  addRC("R1 + R2_TGFb", "R1_R2_TGFb",       "k_form_R1_R2_TGFb * R1 * R2_TGFb") |>
  addRC("R1_R2_TGFb",     "R1_R2_TGFb_int", "k_int_R1_R2_TGFb * R1_R2_TGFb") |>
  addRC("R1_R2_TGFb_int", "R1_R2_TGFb",     "k_rec_R1_R2_TGFb_int * R1_R2_TGFb_int") |>
  addRC("R1_R2_TGFb",     "",               "k_dg_R1_R2_TGFb * R1_R2_TGFb * (1 - bool_MG132)") |>
  addRC("R1 + Smad7", "",             "k_dg_R1_S7 * R1 * Smad7 / (Km_S7_R1 + R1 + R1_R2_TGFb + R1_R2_TGFb_int) * (1 - bool_MG132)") |>
  addRC("R1_R2_TGFb + Smad7", "",     "k_dg_R1_R2_TGFb_S7 * R1_R2_TGFb * Smad7 / (Km_S7_R1 + R1 + R1_R2_TGFb + R1_R2_TGFb_int) * (1 - bool_MG132)") |>
  addRC("R1_R2_TGFb_int + Smad7", "", "k_dg_R1_R2_TGFb_int_S7 * R1_R2_TGFb_int * Smad7 / (Km_S7_R1 + R1 + R1_R2_TGFb + R1_R2_TGFb_int)* (1 - bool_MG132)") |>
  addRC("Smad2", "pSmad2", "k_pS2_R1_R2_TGFb * R1_R2_TGFb * Smad2 / (1 + k_inh_pS2_S7 * Smad7)") |>
  addRC("Smad2", "pSmad2", "k_pS2_R1_R2_TGFb_int * R1_R2_TGFb_int * Smad2 / (1 + k_inh_pS2_S7 * Smad7)") |>
  addRC("pSmad2", "Smad2", "k_dp_S2 * pSmad2") |>
  addRC("Smad3", "pSmad3", "k_pS3_R1_R2_TGFb * R1_R2_TGFb * Smad3 / (1 + k_inh_pS3_S7 * Smad7)") |>
  addRC("Smad3", "pSmad3", "k_pS3_R1_R2_TGFb_int * R1_R2_TGFb_int * Smad3 / (1 + k_inh_pS3_S7 * Smad7)") |>
  addRC("pSmad3", "Smad3", "k_dp_S3 * pSmad3") |>
  addRC("pSmad2 + pSmad3 + Smad4", "C3", "k_form_C3 * pSmad2 * pSmad3 * Smad4") |>
  addRC("C3", "pSmad2 + pSmad3 + Smad4", "k_dec_C3 * C3") |>
  addRC("", "Smad7mRNA", "k_pr_S7mRNA_basal * (1 - bool_ActD)") |>
  addRC("", "Smad7mRNA", "k_pr_S7mRNA * C3 * (1 - bool_ActD)") |>
  addRC("Smad7mRNA", "", "k_dg_S7mRNA * Smad7mRNA") |>
  addRC("", "Smad7",     "k_pr_S7 * Smad7mRNA * (1 - bool_CHX)") |>
  addRC("Smad7", "",     "k_dg_S7 * Smad7 * (1 - bool_MG132)") |>
  addRC("", "FB3mRNA", "k_pr_FB3mRNA * C3 * (1 - bool_ActD)") |>
  addRC("", "FB4mRNA", "k_pr_FB4mRNA * C3 * (1 - bool_ActD)") |>
  addRC("FB3mRNA", "", "k_dg_FB3mRNA * FB3mRNA") |>
  addRC("FB4mRNA", "", "k_dg_FB4mRNA * FB4mRNA") |>
  addRC("", "FB3", "k_pr_FB3 * FB3mRNA * (1 - bool_CHX)") |>
  addRC("", "FB4", "k_pr_FB4 * FB4mRNA * (1 - bool_CHX)") |>
  addRC("FB3", "", "k_dg_FB3 * FB3 * (1 - bool_MG132)") |>
  addRC("FB4", "", "k_dg_FB4 * FB4 * (1 - bool_MG132)") |>
  addRC("", "SERPINE1mRNA", "k_pr_SERPINE1mRNA_basal * (1 - bool_ActD)") |>
  addRC("", "SERPINE1mRNA", "k_pr_SERPINE1mRNA * C3 * (1 - bool_ActD)") |>
  addRC("SERPINE1mRNA", "", "k_dg_SERPINE1mRNA * SERPINE1mRNA") |>
  addRC("", "SERPINE1",     "k_pr_SERPINE1 * SERPINE1mRNA * (1 - bool_CHX)") |>
  addRC("SERPINE1", "",     "k_dg_SERPINE1 * SERPINE1 * (1 - bool_MG132)")

reactions$compartments$Cell$volume      <- "volumeC"
reactions$compartments$extraCell$volume <- "volumeEC"
reactions <- customTotals(reactions, list(
  totalSmad2 = "Smad2 + pSmad2 + C3",
  totalSmad3 = "Smad3 + pSmad3 + C3",
  totalSmad4 = "Smad4 + C3"))

events <- eventlist() |>
  addEvent(var = "TGFb",       time = 0,     value = "var_TGFb_dose",  method = "add") |>
  addEvent(var = "bool_CHX",   time = -30,   value = "var_bool_CHX",   method = "replace") |>
  addEvent(var = "bool_MG132", time = -30,   value = "var_bool_MG132", method = "replace") |>
  addEvent(var = "bool_ActD",  time = -30,   value = "var_bool_ActD",  method = "replace") |>
  addEvent(var = "bool_R1OE",  time = -2880, value = "var_bool_R1OE",  method = "replace") |>
  addEvent(var = "bool_R2OE",  time = -2880, value = "var_bool_R2OE",  method = "replace") |>
  addEvent(var = "bool_R1Knd", time = -2880, value = "var_bool_R1Knd", method = "replace") |>
  addEvent(var = "bool_R2Knd", time = -2880, value = "var_bool_R2Knd", method = "replace")

.forcings <- c("bool_ActD", "bool_CHX", "bool_MG132",
               "bool_R1OE", "bool_R2OE", "bool_R1Knd", "bool_R2Knd")

cat("== steadyStates(simplify = 'full') ==\n")
t_ss <- system.time(
  mysteadies <- steadyStates(reactions, forcings = .forcings, simplify = "full"))
cat(sprintf("steadyStates: %.1fs, %d equations\n", t_ss["elapsed"], length(mysteadies)))

# rationalised observables (log2/offset form reduced to scale * state sums)
.R1tot <- "R1+R1_R2_TGFb+R1_R2_TGFb_int"
.R2tot <- "R2+R2_TGFb+R1_R2_TGFb+R1_R2_TGFb_int"
observables_rational <- eqnvec(
  pSmad2_obs         = "scale_pSmad2*(pSmad2+C3)",
  pSmad2_obs_DR1     = "scale_pSmad2_DR1*(pSmad2+C3)",
  pSmad3_obs         = "scale_pSmad3*(pSmad3+C3)",
  pSmad3_obs_DR1     = "scale_pSmad3_DR1*(pSmad3+C3)",
  SERPINE1_mRNA_obs_RNAseq = "scale_SERPINE1_mRNA_RNAseq*(SERPINE1mRNA)",
  SERPINE1_prot_obs  = "scale_SERPINE1_prot*(SERPINE1)",
  Smad4_CoIP_obs     = "scale_Smad4_CoIP*(C3)",
  Smad4_CoIP_obs_DR1 = "scale_Smad4_CoIP_DR1*(C3)",
  Smad4_CoIP_obs_DR2 = "scale_Smad4_CoIP_DR2*(C3)",
  SMAD7_mRNA_obs        = "scale_SMAD7_mRNA*(Smad7mRNA)",
  SMAD7_mRNA_obs_RNAseq = "scale_SMAD7_mRNA_RNAseq*(Smad7mRNA)",
  TGFb_obs           = "scale_TGFb*(TGFb)",
  TGFb_obs_CHX       = "scale_TGFb_CHX*(TGFb)",
  TGFb_obs_MG132     = "scale_TGFb_MG132*(TGFb)",
  TGFB1_mRNA_obs_RNAseq  = "scale_TGFB1_mRNA_RNAseq*(TGFB1mRNA)",
  TGFBR1_mRNA_obs        = "scale_TGFBR1_mRNA*(R1mRNA)",
  TGFBR1_mRNA_obs_ActD   = "scale_TGFBR1_mRNA_ActD*(R1mRNA)",
  TGFBR1_mRNA_obs_RNAseq = "scale_TGFBR1_mRNA_RNAseq*(R1mRNA)",
  TGFBR1_prot_obs    = paste0("(", .R1tot, ")"),
  TGFBR1_prot_obs_MS = paste0("(", .R1tot, ")"),
  TGFBR2_mRNA_obs        = "scale_TGFBR2_mRNA*(R2mRNA)",
  TGFBR2_mRNA_obs_ActD   = "scale_TGFBR2_mRNA_ActD*(R2mRNA)",
  TGFBR2_mRNA_obs_RNAseq = "scale_TGFBR2_mRNA_RNAseq*(R2mRNA)",
  TGFBR2_prot_obs    = paste0("(", .R2tot, ")"),
  TGFBR2_prot_obs_MS = paste0("(", .R2tot, ")"),
  TSmad2_obs         = "(Smad2+pSmad2+C3)",
  TSmad2_obs_DIA     = "(Smad2+pSmad2+C3)",
  TSmad3_obs         = "(Smad3+pSmad3+C3)",
  TSmad3_obs_DIA     = "(Smad3+pSmad3+C3)",
  TSmad4_obs         = "(Smad4+C3)",
  TSmad4_obs_DIA     = "(Smad4+C3)")

# --- 34 conditions approximating the data layout: 8 TC arms (full panel),
# DR dose ladders (sparse panel), ligand-depletion assay (TGFb panel) ---
mkrow <- function(pert, dose, assay) data.frame(
  Pertubation = pert, init_TGFb = dose, Assay = assay,
  var_bool_ActD  = as.integer(pert == "ActD"),
  var_bool_CHX   = as.integer(pert == "CHX"),
  var_bool_MG132 = as.integer(pert == "MG132"),
  var_bool_R1OE  = as.integer(pert == "OE_R1"),
  var_bool_R2OE  = as.integer(pert == "OE_R2"),
  var_bool_R1Knd = as.integer(pert == "R1Knd"),
  var_bool_R2Knd = as.integer(pert == "R2Knd"),
  var_TGFb_dose  = dose / 25, stringsAsFactors = FALSE)

tcPerts <- c("Ctrl", "ActD", "CHX", "MG132", "R1Knd", "R2Knd", "OE_R1", "OE_R2")
rows <- c(
  lapply(tcPerts, mkrow, dose = 1, assay = "main"),                       #  8 TC
  lapply(c(0.008, 0.04, 0.2, 1, 5, 25), function(d) mkrow("Ctrl", d, "main")),   #  6 DR1
  lapply(c(0.04, 0.2, 1, 5, 25), function(d) mkrow("R1Knd", d, "main")),  #  5 DR
  lapply(c(0.04, 0.2, 1, 5, 25), function(d) mkrow("R2Knd", d, "main")),  #  5 DR
  lapply(c(0.2, 1, 5), function(d) mkrow("OE_R1", d, "main")),            #  3 DR
  lapply(c(0.2, 1, 5), function(d) mkrow("OE_R2", d, "main")),            #  3 DR
  lapply(c("Ctrl", "CHX", "MG132"), mkrow, dose = 1, assay = "depl"),     #  3 depl
  list(mkrow("Ctrl", 5, "depl")))                                         #  1 depl
cond.grid <- do.call(rbind, rows)
rownames(cond.grid) <- sprintf("c%02d_%s_%g_%s", seq_len(nrow(cond.grid)),
                               cond.grid$Pertubation, cond.grid$init_TGFb,
                               cond.grid$Assay)
cat(sprintf("conditions: %d\n", nrow(cond.grid)))

tcObs   <- setdiff(names(observables_rational),
                   c(grep("_DR", names(observables_rational), value = TRUE),
                     "TGFb_obs", "TGFb_obs_CHX", "TGFb_obs_MG132"))
drObs   <- c("pSmad2_obs_DR1", "pSmad3_obs_DR1", "Smad4_CoIP_obs_DR1", "Smad4_CoIP_obs_DR2")
deplObs <- c("TGFb_obs", "TGFb_obs_CHX", "TGFb_obs_MG132")
gSI <- lapply(seq_len(nrow(cond.grid)), function(i) {
  a <- cond.grid$Assay[i]
  isDR <- cond.grid$init_TGFb[i] != 1 || i > 8
  o <- if (a == "depl") deplObs else if (isDR) drObs else tcObs
  as.eqnvec(observables_rational[o])
})
names(gSI) <- rownames(cond.grid)

trafoSI <- lapply(cond.grid$Assay, function(a) as.eqnvec(c(
  mysteadies, volumeEC = paste0("volumeEC_", a))))

.fixToOneSI <- c("k_pr_FB3", "k_pr_FB3mRNA", "k_pr_FB4", "k_pr_FB4mRNA", "volumeC",
                 paste0("scale_", c("TGFB1_mRNA_RNAseq", "TGFBR1_mRNA_RNAseq",
                                    "TGFBR2_mRNA_RNAseq", "SMAD7_mRNA_RNAseq",
                                    "SERPINE1_mRNA_RNAseq", "SERPINE1_prot")))

cat("== symmetryDetection (M009 SI check, 34 conditions, reconstruct) ==\n")
p0 <- proc.time()
idResult <- symmetryDetection(
  f = reactions, g = gSI, trafo = trafoSI, fixed = .fixToOneSI,
  method = "observability", events = events, conditions = cond.grid,
  forcings = .forcings, equilibrate = FALSE, reduceCQ = FALSE,
  reconstruct = TRUE, cores = 10, verbose = FALSE)
p1 <- proc.time() - p0
cat(sprintf("symmetryDetection: %.1fs elapsed | rank %d/%d, %d directions\n",
            p1["elapsed"], idResult$rank, idResult$dim,
            length(idResult$nonIdentifiable)))
print(summary(idResult))
