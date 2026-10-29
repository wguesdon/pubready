#!/usr/bin/env Rscript
# Write the stress datasets of the render audit. Each dataset pushes one recipe
# toward a layout limit that the example data does not reach: long labels, many
# groups, small and large n, extreme values. The data is seeded, so a rerun
# writes the same files. Run inside the pinned image from the repository root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:<version> \
#     Rscript tools/render_audit/make_stress_data.R pubready_output/render_audit/data
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "pubready_output/render_audit/data"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(20261029)

write_csv <- function(df, name) {
  utils::write.csv(df, file.path(out, name), row.names = FALSE)
}

# One correlated variable: a scaled latent factor plus noise.
cv_column <- function(i, latent) {
  latent[, (i %% 3) + 1] * stats::runif(1, 0.3, 1) + stats::rnorm(nrow(latent))
}

# One best response, with a better response more likely at a higher arm index.
draw_response <- function(k) {
  sample(c("complete response", "partial response", "stable disease", "progressive disease"),
         1, prob = c(0.05 * k, 0.1 + 0.04 * k, 0.4, 0.45 - 0.09 * k))
}

# The lesion measurements of one patient, from three to eight visits.
spider_patient <- function(i, arms) {
  weeks <- seq(0, 56, by = 8)[seq_len(sample(3:8, 1))]
  slope <- c(control = 0.5, "drug A" = -0.4, "drug A + B" = -0.9)[[arms[i]]] +
    stats::rnorm(1, 0, 0.4)
  base <- stats::runif(1, 20, 120)
  noise <- stats::rnorm(length(weeks), 0, 0.03)
  data.frame(patient = sprintf("PT-%03d", i), arm = arms[i], week = weeks,
             lesion_mm = round(pmax(0, base * (1 + slope * weeks / 100 + noise)), 1))
}

# One falling dose-response curve in triplicate.
dose_curve <- function(i, cmpds, ec50, doses) {
  d <- rep(doses, each = 3)
  data.frame(compound = cmpds[i], concentration_um = d,
             viability_pct = round(5 + 95 / (1 + (d / ec50[i])^1.1) +
                                     stats::rnorm(length(d), 0, 4), 1))
}

# The tumour volumes of one mouse, eight mice to an arm.
tumour_mouse <- function(i, arms, days) {
  a <- arms[(i - 1) %/% 8 + 1]
  g <- c(0.11, 0.09, 0.07, 0.04)[match(a, arms)]
  data.frame(mouse = sprintf("M%02d", i), arm = a, day = days,
             volume_mm3 = round(100 * exp(g * days) * exp(stats::rnorm(length(days), 0, 0.08)), 1))
}

# Groups with labels of the length a scientist types into a spreadsheet.
long_groups <- c("Vehicle control (DMSO 0.1%)", "Compound 7b, 10 uM (72 h)")
write_csv(data.frame(
  group = rep(long_groups, each = 8),
  value = round(c(rnorm(8, 100, 12), rnorm(8, 72, 10)), 1)),
  "two_group_long_labels.csv")

write_csv(data.frame(group = rep(c("ctrl", "trt"), each = 3),
                     value = c(4.1, 5.3, 4.8, 7.9, 9.2, 8.4)),
          "two_group_small_n.csv")

write_csv(data.frame(group = rep(c("wild type", "knockout"), each = 400),
                     value = round(c(rnorm(400, 50, 9), rnorm(400, 53, 11)), 2)),
          "two_group_large_n.csv")

write_csv(data.frame(group = rep(c("baseline", "stimulated"), each = 12),
                     value = round(c(rnorm(12, -2e5, 6e4), rnorm(12, 1.5e5, 9e4)))),
          "two_group_large_values.csv")

# Eight groups: the number of significant pairs, and so of brackets, is large.
mg_levels <- c("WT", "Het", "KO", "KO + rescue", "KO + vehicle",
               "Double KO", "Overexpression", "Empty vector")
mg_means <- c(10, 11, 16, 10.5, 15.8, 19, 7, 10.2)
write_csv(data.frame(
  genotype = rep(mg_levels, each = 7),
  expression = round(rnorm(56, rep(mg_means, each = 7), 1.1), 2)),
  "multi_group_eight.csv")

fa <- expand.grid(rep = 1:6,
                  treatment = c("Vehicle", "Low dose (1 mg/kg)", "High dose (10 mg/kg)"),
                  strain = c("C57BL/6J", "BALB/c", "NOD scid gamma", "FVB/N"),
                  stringsAsFactors = FALSE)
fa$response <- round(20 + 3 * (fa$strain == "NOD scid gamma") +
                       4 * (fa$treatment != "Vehicle") * (fa$strain != "FVB/N") +
                       rnorm(nrow(fa), 0, 2), 2)
write_csv(fa[, c("strain", "treatment", "response")], "factorial_four_by_three.csv")

# Four arms with long names and a long follow-up.
arms <- c("Placebo", "Drug A 5 mg once daily", "Drug A 10 mg once daily",
          "Drug A 10 mg + standard of care")
rates <- c(0.06, 0.045, 0.035, 0.025)
n_arm <- 70
km_t <- unlist(lapply(rates, stats::rexp, n = n_arm))
km_c <- stats::runif(4 * n_arm, 6, 60)
km <- data.frame(
  months = round(pmin(km_t, km_c), 1),
  event = as.integer(km_t <= km_c),
  arm = rep(arms, each = n_arm),
  age = round(stats::rnorm(4 * n_arm, 63, 9)),
  sex = sample(c("female", "male"), 4 * n_arm, replace = TRUE),
  ecog = sample(c("0", "1", "2"), 4 * n_arm, replace = TRUE, prob = c(0.5, 0.4, 0.1)),
  histology = sample(c("adenocarcinoma", "squamous", "large cell", "other", "mixed"),
                     4 * n_arm, replace = TRUE, prob = c(0.45, 0.3, 0.12, 0.1, 0.03)),
  ldh_u_l = round(stats::rlnorm(4 * n_arm, log(220), 0.4)))
write_csv(km, "survival_four_arms.csv")

# A large heatmap with three groups: no stars, many labels.
n_gene <- 80; n_per <- 12
grp3 <- rep(c("control", "low dose", "high dose"), each = n_per)
samples <- sprintf("%s_%02d", rep(c("ctl", "low", "high"), each = n_per), rep(seq_len(n_per), 3))
mat <- matrix(stats::rnorm(n_gene * 3 * n_per, 8, 1.2), nrow = n_gene)
shift <- stats::rnorm(n_gene, 0, 1.5)
mat[, grp3 == "high dose"] <- mat[, grp3 == "high dose"] + shift
mat[, grp3 == "low dose"] <- mat[, grp3 == "low dose"] + shift / 2
hm <- data.frame(gene = sprintf("ENSG%011d", 1e6 + seq_len(n_gene)), round(mat, 2))
names(hm) <- c("gene", samples)
write_csv(hm, "heatmap_large_matrix.csv")
write_csv(data.frame(sample = samples, group = grp3), "heatmap_large_annotation.csv")

cx <- stats::rnorm(500, 50, 10)
write_csv(data.frame(dose_mg_per_kg = round(cx, 2),
                     plasma_concentration_ng_ml = round(0.8 * cx + stats::rnorm(500, 0, 12), 2)),
          "correlation_large_n.csv")

# Twenty variables with long names for the correlation heatmap.
latent <- matrix(stats::rnorm(120 * 3), ncol = 3)
cv_names <- c("interleukin_6_pg_ml", "tumour_necrosis_factor_pg_ml", "c_reactive_protein_mg_l",
              "interleukin_1_beta_pg_ml", "interleukin_10_pg_ml", "tgf_beta_ng_ml",
              "interleukin_4_pg_ml", "interferon_gamma_pg_ml", "ferritin_ng_ml",
              "procalcitonin_ng_ml", "d_dimer_ug_ml", "fibrinogen_g_l",
              "lymphocyte_count", "neutrophil_count", "platelet_count",
              "albumin_g_l", "ldh_u_l", "creatinine_umol_l", "troponin_ng_l", "lactate_mmol_l")
cv <- sapply(seq_along(cv_names), cv_column, latent = latent)
colnames(cv) <- cv_names
write_csv(as.data.frame(round(cv, 3)), "correlation_heatmap_twenty.csv")

# Seven sets and 300 elements for the UpSet plot.
sets <- c("apoptosis", "inflammation", "metabolism", "cell_cycle", "dna_repair",
          "angiogenesis", "immune_response")
um <- as.data.frame(matrix(stats::rbinom(300 * 7, 1, 0.3), ncol = 7))
names(um) <- sets
write_csv(data.frame(gene = sprintf("gene_%03d", 1:300), um), "upset_seven_sets.csv")

# A DE table with p = 0, very large fold changes and 20000 genes.
ng <- 20000
lfc <- stats::rnorm(ng, 0, 0.6)
hit <- sample(ng, 600)
lfc[hit] <- lfc[hit] + sample(c(-1, 1), 600, replace = TRUE) * stats::rexp(600, 0.4)
z <- abs(lfc) / 0.25 + stats::rnorm(ng, 0, 0.5)
pv <- 2 * stats::pnorm(-abs(z))
pv[order(-abs(lfc))[1:15]] <- 0
de <- data.frame(gene = sprintf("Gene%05d", seq_len(ng)), log2FoldChange = round(lfc, 4),
                 pvalue = signif(pv, 4), padj = signif(stats::p.adjust(pv, "BH"), 4))
write_csv(de, "volcano_extreme.csv")

go_terms <- c(
  "regulation of transcription by RNA polymerase II",
  "positive regulation of cytokine production involved in inflammatory response",
  "negative regulation of extrinsic apoptotic signaling pathway via death domain receptors",
  "antigen processing and presentation of exogenous peptide antigen via MHC class II",
  "mitochondrial electron transport, NADH to ubiquinone",
  "extracellular matrix organization", "T cell receptor signaling pathway",
  "response to unfolded protein", "DNA damage response, signal transduction by p53 class mediator",
  "Wnt signaling pathway, planar cell polarity pathway", "angiogenesis",
  "fatty acid beta-oxidation using acyl-CoA dehydrogenase", "glycolytic process",
  "type I interferon signaling pathway", "neutrophil degranulation",
  "cellular response to tumor necrosis factor", "collagen fibril organization",
  "ribosome biogenesis", "mRNA splicing, via spliceosome", "oxidative phosphorylation",
  "complement activation, classical pathway", "cell cycle G2/M phase transition",
  "autophagosome assembly", "regulation of TOR signaling", "cholesterol biosynthetic process")
cnt <- sample(8:140, length(go_terms))
write_csv(data.frame(Description = go_terms, GeneRatio = sprintf("%d/1500", cnt), Count = cnt,
                     p.adjust = signif(10^-stats::runif(length(go_terms), 2, 18), 3)),
          "enrichment_long_terms.csv")

pid <- sprintf("S%03d", 1:120)
pre <- stats::rnorm(120, 100, 15)
write_csv(data.frame(subject = rep(pid, 2), condition = rep(c("baseline", "week 12"), each = 120),
                     value = round(c(pre, pre - 4 + stats::rnorm(120, 0, 8)), 1)),
          "paired_many_subjects.csv")

parms <- c("Placebo", "Low dose", "Medium dose", "High dose", "High dose + adjuvant")
pr <- data.frame(arm = sample(parms, 400, replace = TRUE))
pr$response <- vapply(match(pr$arm, parms), draw_response, character(1))
write_csv(pr, "proportions_five_by_four.csv")

pc_groups <- c("healthy", "mild", "moderate", "severe", "critical", "recovered")
pcg <- rep(pc_groups, each = 20)
pcm <- matrix(stats::rnorm(120 * 40), nrow = 120)
pcm[, 1:5] <- pcm[, 1:5] + match(pcg, pc_groups) * 0.7
write_csv(data.frame(sample = sprintf("S%03d", 1:120), group = pcg, round(pcm, 3)),
          "pca_six_groups.csv")

# Sixty patients over eight visits in three arms for the spider and waterfall plots.
sp_arm <- rep(c("control", "drug A", "drug A + B"), each = 20)
sp <- do.call(rbind, lapply(seq_along(sp_arm), spider_patient, arms = sp_arm))
write_csv(sp, "spider_sixty_patients.csv")

doses <- c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
cmpds <- c("Cmpd-001", "Cmpd-002", "Cmpd-003", "Cmpd-004", "Cmpd-005", "Inactive analogue")
ec50 <- c(0.02, 0.08, 0.4, 1.5, 6, 1e4)
dr <- do.call(rbind, lapply(seq_along(cmpds), dose_curve, cmpds = cmpds, ec50 = ec50,
                            doses = doses))
write_csv(dr, "dose_response_six.csv")

ry <- rep(c("case", "control"), c(90, 110))
rcase <- ry == "case"
write_csv(data.frame(
  outcome = ry,
  troponin = round(stats::rnorm(200, ifelse(rcase, 2, 0), 1), 3),
  nt_probnp = round(stats::rnorm(200, ifelse(rcase, 1, 0), 1), 3),
  crp = round(stats::rnorm(200, ifelse(rcase, 0.3, 0), 1), 3),
  inverse_marker = round(stats::rnorm(200, ifelse(rcase, -0.8, 0), 1), 3)),
  "roc_four_markers.csv")

ba_ref <- stats::runif(300, 5, 500)
write_csv(data.frame(sample = sprintf("X%03d", 1:300), reference = round(ba_ref, 1),
                     candidate = round(ba_ref * 1.06 + 2 + stats::rnorm(300, 0, 0.04 * ba_ref), 1)),
          "bland_altman_large.csv")

ms_n <- 25
ms_hr <- exp(stats::rnorm(ms_n, log(0.8), 0.25))
ms_se <- stats::runif(ms_n, 0.08, 0.6)
ms_se[c(3, 17)] <- c(1.2, 1.5)
ms_names <- sprintf("%s et al. (%d)", sample(c("Smith", "Garcia-Fernandez", "Nakamura", "O'Brien",
                                                "van der Berg", "Kowalczyk", "Li", "Okonkwo"),
                                              ms_n, replace = TRUE), sample(1998:2024, ms_n, replace = TRUE))
write_csv(data.frame(study = make.unique(ms_names, sep = " "), hr = round(ms_hr, 3),
                     lower = round(ms_hr * exp(-1.96 * ms_se), 3),
                     upper = round(ms_hr * exp(1.96 * ms_se), 3)),
          "forest_meta_twenty_five.csv")

sw_n <- 70
sw_dur <- round(stats::rexp(sw_n, 1 / 14) + 0.5, 1)
sw_resp <- ifelse(stats::runif(sw_n) < 0.5, round(stats::runif(sw_n, 0.2, 0.6) * sw_dur, 1), NA)
sw_prog <- ifelse(stats::runif(sw_n) < 0.6, sw_dur, NA)
write_csv(data.frame(patient = sprintf("PT-%03d", seq_len(sw_n)),
                     arm = sample(c("Arm A", "Arm B", "Arm C"), sw_n, replace = TRUE),
                     months = sw_dur, ended = as.integer(stats::runif(sw_n) < 0.75),
                     response_month = sw_resp, progression_month = sw_prog),
          "swimmer_seventy.csv")

tc_days <- c(0, 2, 4, 7, 9, 11, 14, 16, 18, 21)
tc_arms <- c("vehicle", "drug 1 mg/kg", "drug 3 mg/kg", "drug 10 mg/kg")
tc <- do.call(rbind, lapply(seq_len(32), tumour_mouse, arms = tc_arms, days = tc_days))
write_csv(tc, "timecourse_four_arms.csv")

cat(sprintf("Wrote the stress datasets to %s\n", out))
