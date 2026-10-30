# Benchmarks of the v3 numbers audit: 32 cases on public datasets, each checked
# against a reference value that an independent call to the canonical R function
# computed (benchmarks/reference_values.csv, sections 40 to 62 of
# tools/make_benchmarks/make_benchmark_data.R). `figkit test` mounts the
# benchmarks folder; the tests skip when it is absent. testthat runs each file in
# its own environment, so the helpers of test-benchmarks.R are repeated here.

bench_dir <- function() Sys.getenv("PUBREADY_BENCHMARKS", "")

bench_skip <- function() {
  skip_if(!nzchar(bench_dir()) || !dir.exists(bench_dir()), "benchmarks not mounted")
}

bench_data <- function(name) read_tidy(file.path(bench_dir(), "data", paste0(name, ".csv")))

bench_spec <- function(recipe, ...) {
  opt <- utils::modifyList(list(recipe = recipe, data = "bench.csv", test = "auto"), list(...))
  spec_from_opt(opt)
}

bench_refs <- function() {
  utils::read.csv(file.path(bench_dir(), "reference_values.csv"), stringsAsFactors = FALSE)
}

expect_ref <- function(actual, case, quantity) {
  refs <- bench_refs()
  r <- refs[refs$case == case & refs$quantity == quantity, , drop = FALSE]
  expect_equal(nrow(r), 1, info = paste(case, quantity))
  expect_equal(as.numeric(actual), r$value, tolerance = r$rel_tol, info = paste(case, quantity))
}

pair_row <- function(stats, a, b) {
  stats[(stats$group1 == a & stats$group2 == b) | (stats$group1 == b & stats$group2 == a), ]
}

# The count that a "Dropped <n> row(s) with <reason>" cleaning step names.
dropped_rows <- function(steps, reason) {
  hit <- grep(reason, steps, value = TRUE, fixed = TRUE)
  if (length(hit) != 1) return(NA_real_)
  as.numeric(sub("^Dropped ([0-9]+) row.*$", "\\1", hit))
}

# --- two_group_compare -------------------------------------------------------

test_that("cats: Student's t-test and the pooled-SD Cohen's d match t.test(var.equal = TRUE)", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("cats"),
                                  bench_spec("two_group_compare", x = "Sex", y = "Hwt", test = "student_t"))
  expect_equal(res$test_meta$name, "Student's two-sample t-test")
  expect_equal(res$stats$group1, "F")
  expect_ref(res$stats$statistic, "cats_student", "statistic")
  expect_ref(res$stats$df, "cats_student", "df")
  expect_ref(res$stats$p_value, "cats_student", "p_value")
  expect_ref(res$stats$effect_size, "cats_student", "cohens_d")
})

test_that("cats: the Mann-Whitney U test matches wilcox.test(correct = TRUE)", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("cats"),
                                  bench_spec("two_group_compare", x = "Sex", y = "Hwt", test = "wilcoxon"))
  expect_equal(res$test_meta$name, "Mann-Whitney U test")
  expect_ref(res$stats$statistic, "cats_mann_whitney", "statistic")
  expect_ref(res$stats$p_value, "cats_mann_whitney", "p_value")
  expect_ref(res$stats$effect_size, "cats_mann_whitney", "rank_biserial")
})

test_that("PlantGrowth ctrl and trt2: the exact Mann-Whitney U test matches wilcox.test(exact = TRUE)", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("plant_growth_ctrl_trt2"),
                                  bench_spec("two_group_compare", x = "group", y = "weight", test = "wilcoxon"))
  expect_ref(res$stats$statistic, "plant_growth_mw_exact", "statistic")
  expect_ref(res$stats$p_value, "plant_growth_mw_exact", "p_value")
  expect_ref(res$stats$effect_size, "plant_growth_mw_exact", "rank_biserial")
})

test_that("sleep through two_group_compare --paired: the paired t-test and dz match t.test(paired = TRUE)", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("sleep"),
                                  bench_spec("two_group_compare", x = "group", y = "extra", paired = TRUE,
                                             test = "t"))
  expect_equal(res$test_meta$name, "paired t-test")
  expect_equal(res$stats$effect_size_name, "Cohen's dz")
  expect_ref(res$stats$statistic, "sleep_two_group_paired", "statistic")
  expect_ref(res$stats$df, "sleep_two_group_paired", "df")
  expect_ref(res$stats$p_value, "sleep_two_group_paired", "p_value")
  expect_ref(res$stats$effect_size, "sleep_two_group_paired", "cohens_dz")
})

test_that("sleep through two_group_compare --paired: auto checks the differences and picks the signed-rank test", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("sleep"),
                                  bench_spec("two_group_compare", x = "group", y = "extra", paired = TRUE))
  expect_ref(res$stats$shapiro_p_diff, "sleep_two_group_paired_auto", "shapiro_p_diff")
  expect_ref(res$test_meta$assumptions$normality_of_differences$p, "sleep_two_group_paired_auto",
             "shapiro_p_diff")
  expect_equal(res$test_meta$name, "Wilcoxon signed-rank test")
  expect_ref(res$stats$statistic, "sleep_two_group_paired_auto", "statistic")
  expect_ref(res$stats$p_value, "sleep_two_group_paired_auto", "p_value")
})

test_that("nitrofen 80 and 160: the numeric doses sort as numbers and Welch's t-test matches t.test", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("nitrofen_80_160"),
                                  bench_spec("two_group_compare", x = "conc", y = "total", test = "welch_t"))
  expect_equal(c(res$stats$group1, res$stats$group2), c("80", "160"))
  expect_ref(res$stats$statistic, "nitrofen_two_group", "statistic")
  expect_ref(res$stats$p_value, "nitrofen_two_group", "p_value")
})

test_that("Cars93: USA is group 1 and the Mann-Whitney U test matches wilcox.test", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("cars93"),
                                  bench_spec("two_group_compare", x = "Origin", y = "MPG.highway",
                                             test = "wilcoxon"))
  expect_equal(c(res$stats$group1, res$stats$group2), c("USA", "non-USA"))
  expect_ref(res$stats$statistic, "cars93_origin", "statistic")
  expect_ref(res$stats$p_value, "cars93_origin", "p_value")
  expect_ref(res$stats$effect_size, "cars93_origin", "rank_biserial")
})

# --- feature_panel -----------------------------------------------------------

test_that("cats: the feature panel Welch tests, the pooled-SD d and BH match t.test and p.adjust", {
  bench_skip()
  res <- recipe_feature_panel(bench_data("cats"), bench_spec("feature_panel", x = "Sex", test = "welch_t"))
  for (f in c("Bwt", "Hwt")) {
    r <- res$stats[res$stats$feature == f, ]
    expect_ref(r$statistic, "cats_feature_panel", paste0("statistic_", f))
    expect_ref(r$df1, "cats_feature_panel", paste0("df_", f))
    expect_ref(r$p_value, "cats_feature_panel", paste0("p_", f))
    expect_ref(r$effect_size, "cats_feature_panel", paste0("cohens_d_", f))
    expect_ref(r$p_adj, "cats_feature_panel", paste0("p_adj_", f))
  }
})

# --- multi_group_compare -----------------------------------------------------

test_that("PlantGrowth: every Tukey pair and eta squared match TukeyHSD and anova", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("plant_growth"),
                                    bench_spec("multi_group_compare", x = "group", y = "weight", test = "anova"))
  for (pr in list(c("ctrl", "trt1"), c("ctrl", "trt2"), c("trt1", "trt2"))) {
    r <- pair_row(res$stats, pr[1], pr[2])
    key <- paste0(pr[2], "_", pr[1])
    expect_ref(r$estimate, "plant_growth_tukey", paste0("estimate_", key))
    expect_ref(r$p_adj, "plant_growth_tukey", paste0("p_", key))
  }
  expect_ref(res$test_meta$omnibus$effect_size$value, "plant_growth_tukey", "eta_squared")
})

test_that("chickwts: every Games-Howell pair matches its definition with ptukey", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("chickwts"),
                                    bench_spec("multi_group_compare", x = "feed", y = "weight", test = "welch"))
  expect_equal(res$test_meta$post_hoc, "the Games-Howell test")
  pairs <- utils::combn(c("casein", "horsebean", "linseed", "meatmeal", "soybean", "sunflower"), 2)
  for (k in seq_len(ncol(pairs))) {
    r <- pair_row(res$stats, pairs[1, k], pairs[2, k])
    key <- paste0(pairs[2, k], "_", pairs[1, k])
    expect_ref(r$estimate, "chickwts_games_howell", paste0("estimate_", key))
    expect_ref(r$p_adj, "chickwts_games_howell", paste0("p_", key))
  }
})

insect_pairs <- utils::combn(c("A", "B", "C", "D", "E", "F"), 2)

test_that("InsectSprays: the Dunn BH P values and eta squared (H) match the reference", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("insect_sprays"),
                                    bench_spec("multi_group_compare", x = "spray", y = "count", test = "kruskal"))
  for (k in seq_len(ncol(insect_pairs))) {
    r <- pair_row(res$stats, insect_pairs[1, k], insect_pairs[2, k])
    expect_ref(r$p_adj, "insect_sprays_dunn", sprintf("p_adj_%s_%s", insect_pairs[1, k], insect_pairs[2, k]))
  }
  expect_ref(res$test_meta$omnibus$effect_size$value, "insect_sprays_dunn", "eta_squared_h")
})

test_that("InsectSprays: the Dunn z matches the reference", {
  bench_skip()
  skip("insect_sprays_dunn: the Dunn table has no z column (N-11)")
  res <- recipe_multi_group_compare(bench_data("insect_sprays"),
                                    bench_spec("multi_group_compare", x = "spray", y = "count", test = "kruskal"))
  for (k in seq_len(ncol(insect_pairs))) {
    r <- res$stats[res$stats$group1 == insect_pairs[1, k] & res$stats$group2 == insect_pairs[2, k], ]
    expect_ref(r$statistic, "insect_sprays_dunn", sprintf("z_%s_%s", insect_pairs[1, k], insect_pairs[2, k]))
  }
})

test_that("nitrofen: the doses sort as numbers, and the ANOVA and Tukey match anova and TukeyHSD", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("nitrofen"),
                                    bench_spec("multi_group_compare", x = "conc", y = "total", test = "anova"))
  expect_equal(unique(res$stats$group1), c("0", "80", "160", "235"))
  expect_ref(res$test_meta$omnibus$statistic, "nitrofen_multi", "statistic")
  expect_ref(res$test_meta$omnibus$p_value, "nitrofen_multi", "p_value")
  r <- pair_row(res$stats, "80", "160")
  expect_ref(r$estimate, "nitrofen_multi", "estimate_160_80")
  expect_ref(r$p_adj, "nitrofen_multi", "p_160_80")
  r <- pair_row(res$stats, "0", "80")
  expect_ref(r$estimate, "nitrofen_multi", "estimate_80_0")
  expect_ref(r$p_adj, "nitrofen_multi", "p_80_0")
})

# --- factorial_anova ---------------------------------------------------------

quine_terms <- c("Eth", "Sex", "Lrn", "Eth:Sex", "Eth:Lrn", "Sex:Lrn", "Eth:Sex:Lrn")

quine_spec <- function(test) {
  bench_spec("factorial_anova", x = "Eth", fill = "Sex", facet = "Lrn", y = "Days", test = test)
}

test_that("quine: the unbalanced three-way ANOVA matches car::Anova(type = 2)", {
  bench_skip()
  res <- recipe_factorial_anova(bench_data("quine"), quine_spec("anova"))
  for (term in quine_terms) {
    r <- res$stats[res$stats$effect == term, ]
    expect_ref(r$statistic, "quine_threeway_anova", paste0("F_", term))
    expect_ref(r$p_value, "quine_threeway_anova", paste0("p_", term))
    expect_ref(r$effect_size, "quine_threeway_anova", paste0("pes_", term))
  }
})

test_that("quine: the three-way aligned rank transform matches ARTool", {
  bench_skip()
  res <- recipe_factorial_anova(bench_data("quine"), quine_spec("art"))
  for (term in quine_terms) {
    r <- res$stats[res$stats$effect == term, ]
    expect_ref(r$statistic, "quine_threeway_art", paste0("F_", term))
    expect_ref(r$p_value, "quine_threeway_art", paste0("p_", term))
  }
})

test_that("quine: auto checks the residuals and takes the aligned rank transform", {
  bench_skip()
  res <- recipe_factorial_anova(bench_data("quine"), quine_spec("auto"))
  expect_ref(res$test_meta$residual_normality$p, "quine_threeway_auto", "residual_shapiro_p")
  expect_ref(res$stats$residual_shapiro_p[1], "quine_threeway_auto", "residual_shapiro_p")
  expect_equal(res$test_meta$type, "aligned_rank_transform")
  r <- res$stats[res$stats$effect == "Eth:Sex:Lrn", ]
  expect_ref(r$statistic, "quine_threeway_art", "F_Eth:Sex:Lrn")
  expect_ref(r$p_value, "quine_threeway_art", "p_Eth:Sex:Lrn")
})

test_that("Higgins1990Table1: the two-way aligned rank transform matches ARTool", {
  bench_skip()
  res <- recipe_factorial_anova(bench_data("higgins1990_table1"),
                                bench_spec("factorial_anova", x = "Row", fill = "Column", y = "Response",
                                           test = "art"))
  for (term in c("Row", "Column", "Row:Column")) {
    r <- res$stats[res$stats$effect == term, ]
    expect_ref(r$statistic, "higgins_art", paste0("F_", term))
    expect_ref(r$p_value, "higgins_art", paste0("p_", term))
  }
})

# --- estimation_plot and paired_compare --------------------------------------

test_that("shoes: the paired estimation matches t.test(paired = TRUE) and the DABEST Hedges' g", {
  bench_skip()
  res <- recipe_estimation_plot(bench_data("shoes"),
                                bench_spec("estimation_plot", x = "material", y = "wear", paired = TRUE,
                                           id = "boy", ref = "A"))
  r <- res$stats[res$stats$group == "B", ]
  expect_equal(nrow(r), 1)
  expect_ref(r$mean_diff, "shoes_paired_estimation", "mean_diff")
  expect_ref(r$statistic, "shoes_paired_estimation", "statistic")
  expect_ref(r$df, "shoes_paired_estimation", "df")
  expect_ref(r$p_value, "shoes_paired_estimation", "p_value")
  expect_ref(r$hedges_g, "shoes_paired_estimation", "hedges_g")
})

test_that("sleep through paired_compare: dz and the mean difference are condition 2 minus condition 1", {
  bench_skip()
  res <- recipe_paired_compare(bench_data("sleep"),
                               bench_spec("paired_compare", x = "group", y = "extra", id = "ID", test = "t"))
  expect_ref(res$stats$effect_size, "sleep_paired_effect", "cohens_dz")
  expect_ref(res$stats$mean_diff, "sleep_paired_effect", "mean_diff")
})

# --- timecourse --------------------------------------------------------------

test_that("Oxboys: the one-way repeated-measures ANOVA matches aov(Error()) and car::Anova", {
  bench_skip()
  res <- recipe_timecourse(bench_data("oxboys"),
                           bench_spec("timecourse", x = "Occasion", y = "height", id = "Subject"))
  r <- res$stats[res$stats$effect == "Occasion", ]
  expect_equal(nrow(r), 1)
  expect_ref(r$F, "oxboys_rm", "F_time")
  expect_ref(r$p_value, "oxboys_rm", "p_time_uncorrected")
  expect_ref(r$partial_eta2, "oxboys_rm", "pes_time")
  expect_ref(r$gg_epsilon, "oxboys_rm", "gg_epsilon")
  expect_ref(r$p_gg, "oxboys_rm", "p_time_gg")
  expect_ref(res$test_meta$sphericity$p_value, "oxboys_rm", "mauchly_p")
})

# --- survival ----------------------------------------------------------------

test_that("veteran by arm: the medians, their CIs and the log-rank test match survfit and survdiff", {
  bench_skip()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_survival_km(bench_data("veteran"),
                            bench_spec("survival_km", time = "time", event = "status", x = "trt"))
  for (arm in c("standard", "test")) {
    r <- res$stats[res$stats$group == arm, ]
    expect_ref(r$median_survival, "veteran_trt_km", paste0("median_", arm))
    expect_ref(r$ci_low, "veteran_trt_km", paste0("median_ci_low_", arm))
    expect_ref(r$ci_high, "veteran_trt_km", paste0("median_ci_high_", arm))
  }
  expect_ref(res$test_meta$log_rank$chisq, "veteran_trt_km", "log_rank_chisq")
  expect_ref(res$test_meta$log_rank$p_value, "veteran_trt_km", "log_rank_p")
})

test_that("colon (death): the Cox model with the three-level rx matches survival::coxph", {
  bench_skip()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_cox_forest(bench_data("colon_death"),
                           bench_spec("cox_forest", time = "time", event = "status",
                                      covariates = "rx,sex,age"))
  keys <- c("rxLev+5FU" = "rxLev5FU", rxObs = "rxObs", sex = "sex", age = "age")
  for (term in names(keys)) {
    r <- res$stats[res$stats$term == term, ]
    expect_ref(r$HR, "colon_cox", paste0("HR_", keys[[term]]))
    expect_ref(r$ci_low, "colon_cox", paste0("ci_low_", keys[[term]]))
    expect_ref(r$ci_high, "colon_cox", paste0("ci_high_", keys[[term]]))
    expect_ref(r$p_value, "colon_cox", paste0("p_", keys[[term]]))
  }
  expect_ref(res$test_meta$proportional_hazards$p_value, "colon_cox", "ph_global_p")
})

test_that("cgd: X-linked is the reference of inherit and the Cox model matches survival::coxph", {
  bench_skip()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_cox_forest(bench_data("cgd_first"),
                           bench_spec("cox_forest", time = "time", event = "status",
                                      covariates = "treat,inherit,age"))
  keys <- c("treatrIFN-g" = "treatrIFNg", inheritautosomal = "inheritautosomal", age = "age")
  for (term in names(keys)) {
    r <- res$stats[res$stats$term == term, ]
    expect_ref(r$HR, "cgd_cox", paste0("HR_", keys[[term]]))
    expect_ref(r$p_value, "cgd_cox", paste0("p_", keys[[term]]))
  }
  expect_ref(res$test_meta$proportional_hazards$p_value, "cgd_cox", "ph_global_p")
})

# --- diagnostic and agreement ------------------------------------------------

test_that("birthwt: the AUC, the DeLong CI and the Youden point of lwt match pROC", {
  bench_skip()
  res <- recipe_roc_curve(bench_data("birthwt"),
                          bench_spec("roc_curve", y = "low", x = "lwt", positive = "1"))
  r <- res$stats[res$stats$marker == "lwt", ]
  expect_equal(r$direction, ">")
  expect_ref(r$auc, "birthwt_roc", "auc")
  expect_ref(r$auc_ci_low, "birthwt_roc", "auc_ci_low")
  expect_ref(r$auc_ci_high, "birthwt_roc", "auc_ci_high")
  expect_ref(r$youden_threshold, "birthwt_roc", "youden_threshold")
  expect_ref(r$sensitivity, "birthwt_roc", "youden_sensitivity")
  expect_ref(r$specificity, "birthwt_roc", "youden_specificity")
})

test_that("Davis: the Bland-Altman limits, their CIs and the proportional bias match the reference", {
  bench_skip()
  res <- recipe_bland_altman(bench_data("davis"), bench_spec("bland_altman", x = "weight", y = "repwt"))
  for (q in c("n", "bias", "bias_ci_low", "bias_ci_high", "sd_diff", "loa_low", "loa_high",
              "loa_low_ci_low", "loa_high_ci_high", "prop_bias_slope", "prop_bias_p")) {
    expect_ref(res$stats[[q]], "davis_bland_altman", q)
  }
})

# --- curve fits --------------------------------------------------------------

test_that("Nitrendipene: the falling 4PL curves, their CIs and the shared IC50 test match nls(SSfpl)", {
  bench_skip()
  res <- recipe_dose_response(bench_data("nitrendipene"),
                              bench_spec("dose_response", x = "NIF", y = "activity", group = "Tissue"))
  expect_ref(dropped_rows(res$clean_steps, "with a dose of zero or below"), "nitrendipene_ic50",
             "n_dropped_zero_dose")
  cols <- c(bottom = "bottom", top = "top", ic50 = "midpoint", hill = "hill",
            ic50_ci_low = "midpoint_ci_low", ic50_ci_high = "midpoint_ci_high", rss = "rss")
  for (t in paste0("tissue", 1:4)) {
    r <- res$stats[res$stats$group == t, ]
    for (q in names(cols)) expect_ref(r[[cols[[q]]]], "nitrendipene_ic50", paste0(q, "_", t))
  }
  expect_ref(res$test_meta$shared_ec50_test$statistic, "nitrendipene_ic50", "shared_F")
  expect_ref(res$test_meta$shared_ec50_test$p_value, "nitrendipene_ic50", "shared_p")
})

test_that("DNase run 1: the standard error and the 95% CI of the EC50 match nls(SSfpl)", {
  bench_skip()
  res <- recipe_dose_response(bench_data("dnase_run1"), bench_spec("dose_response", x = "conc", y = "density"))
  expect_ref(res$stats$se_log_midpoint, "dnase_ec50_ci", "se_log_ec50")
  expect_ref(res$stats$midpoint_ci_low, "dnase_ec50_ci", "ec50_ci_low")
  expect_ref(res$stats$midpoint_ci_high, "dnase_ec50_ci", "ec50_ci_high")
})

test_that("Formaldehyde: the linear standard curve and its inverse match lm", {
  bench_skip()
  res <- recipe_standard_curve(bench_data("formaldehyde_standard_curve"),
                               bench_spec("standard_curve", x = "carb", y = "optden", id = "sample",
                                          test = "linear"))
  expect_ref(res$test_meta$parameters$intercept, "formaldehyde_linear", "intercept")
  expect_ref(res$test_meta$parameters$slope, "formaldehyde_linear", "slope")
  expect_ref(res$test_meta$r_squared, "formaldehyde_linear", "r_squared")
  expect_ref(res$test_meta$n_in_range, "formaldehyde_linear", "n_in_range")
  expect_ref(res$test_meta$n_below, "formaldehyde_linear", "n_below")
  refs <- bench_refs()
  conc_q <- refs$quantity[refs$case == "formaldehyde_linear" & grepl("^conc_", refs$quantity)]
  expect_length(conc_q, 5)
  for (q in conc_q) {
    expect_ref(res$stats$concentration[res$stats$id == sub("^conc_", "", q)], "formaldehyde_linear", q)
  }
})

test_that("Puromycin: the 95% CIs of Vmax and Km match nls", {
  bench_skip()
  res <- recipe_enzyme_kinetics(bench_data("puromycin"),
                                bench_spec("enzyme_kinetics", x = "conc", y = "rate", group = "state"))
  for (st in c("treated", "untreated")) {
    r <- res$stats[res$stats$group == st, ]
    for (q in c("vmax_ci_low", "vmax_ci_high", "km_ci_low", "km_ci_high")) {
      expect_ref(r[[q]], "puromycin_ci", paste0(q, "_", st))
    }
  }
})

# --- composition, PCA and correlation ----------------------------------------

test_that("housing by influence: the per-category Kruskal-Wallis tests match kruskal.test and BH", {
  bench_skip()
  res <- recipe_composition(bench_data("housing"),
                            bench_spec("composition", id = "sample", x = "Sat", y = "Freq", group = "Infl"))
  for (s in c("Low", "Medium", "High")) {
    r <- res$stats[res$stats$category == s, ][1, ]
    expect_ref(r$statistic, "housing_kruskal", paste0("H_", s))
    expect_ref(r$p_value, "housing_kruskal", paste0("p_", s))
    expect_ref(r$p_adj, "housing_kruskal", paste0("p_adj_", s))
  }
})

test_that("airway and iris: the PERMANOVA pseudo-F and the exact airway P match the distance matrix", {
  bench_skip()
  ap <- recipe_pca(bench_data("airway_pca"), bench_spec("pca", group = "dex"))
  expect_equal(ap$test_meta$permanova$method, "exact")
  expect_ref(ap$test_meta$permanova$statistic_F, "airway_pca_permanova", "F")
  expect_ref(ap$test_meta$permanova$p_value, "airway_pca_permanova", "p_exact")
  ir <- recipe_pca(bench_data("iris"), bench_spec("pca", group = "Species"))
  expect_ref(ir$test_meta$permanova$statistic_F, "iris_pca_permanova", "F")
})

test_that("mtcars: the 95% CI of the Pearson r matches cor.test", {
  bench_skip()
  res <- recipe_correlation(bench_data("mtcars"),
                            bench_spec("correlation", x = "wt", y = "mpg", test = "pearson"))
  expect_ref(res$stats$conf_low, "mtcars_pearson_ci", "conf_low")
  expect_ref(res$stats$conf_high, "mtcars_pearson_ci", "conf_high")
})
