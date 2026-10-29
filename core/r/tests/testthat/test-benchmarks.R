# Benchmarks: every recipe on a public dataset, checked against a reference value
# that an independent call to the canonical R function computed
# (benchmarks/reference_values.csv, written by tools/make_benchmarks). `figkit
# test` mounts the benchmarks folder; the tests skip when it is absent.

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

var_pair_row <- function(stats, a, b) {
  stats[(stats$var1 == a & stats$var2 == b) | (stats$var1 == b & stats$var2 == a), ]
}

test_that("sleep (Student 1908): the paired t-test matches stats::t.test", {
  bench_skip()
  res <- recipe_paired_compare(bench_data("sleep"),
                               bench_spec("paired_compare", x = "group", y = "extra", id = "ID", test = "t"))
  expect_ref(res$stats$statistic, "sleep_paired", "statistic")
  expect_ref(res$stats$df, "sleep_paired", "df")
  expect_ref(res$stats$p_value, "sleep_paired", "p_value")
})

test_that("sleep: the Welch t-test matches stats::t.test", {
  bench_skip()
  res <- recipe_two_group_compare(bench_data("sleep"),
                                  bench_spec("two_group_compare", x = "group", y = "extra", test = "welch_t"))
  expect_ref(res$stats$statistic, "sleep_welch", "statistic")
  expect_ref(res$stats$df, "sleep_welch", "df")
  expect_ref(res$stats$p_value, "sleep_welch", "p_value")
})

test_that("PlantGrowth: auto picks the one-way ANOVA and Tukey matches TukeyHSD", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("plant_growth"),
                                    bench_spec("multi_group_compare", x = "group", y = "weight"))
  expect_equal(res$test_meta$name, "one-way ANOVA")
  expect_ref(res$test_meta$omnibus$statistic, "plant_growth_anova", "statistic")
  expect_ref(res$test_meta$omnibus$p_value, "plant_growth_anova", "p_value")
  expect_ref(pair_row(res$stats, "trt1", "trt2")$p_adj, "plant_growth_anova", "tukey_p_trt1_trt2")
})

test_that("InsectSprays: the Kruskal-Wallis test matches stats::kruskal.test", {
  bench_skip()
  res <- recipe_multi_group_compare(bench_data("insect_sprays"),
                                    bench_spec("multi_group_compare", x = "spray", y = "count", test = "kruskal"))
  expect_ref(res$test_meta$omnibus$statistic, "insect_sprays_kruskal", "statistic")
  expect_ref(res$test_meta$omnibus$p_value, "insect_sprays_kruskal", "p_value")
})

test_that("ToothGrowth: the two-way ANOVA matches car::Anova at full precision", {
  bench_skip()
  res <- recipe_factorial_anova(bench_data("tooth_growth"),
                                bench_spec("factorial_anova", x = "dose", fill = "supp", y = "len", test = "anova"))
  for (term in c("dose", "supp", "dose:supp")) {
    r <- res$stats[res$stats$effect == term, ]
    expect_ref(r$statistic, "tooth_growth_twoway", paste0("F_", term))
    expect_ref(r$p_value, "tooth_growth_twoway", paste0("p_", term))
  }
  expect_true(all(res$stats$effect_size_name == "partial eta squared"))
})

test_that("lung (NCCTG): the log-rank test matches survival::survdiff", {
  bench_skip()
  res <- recipe_survival_km(bench_data("lung"),
                            bench_spec("survival_km", time = "time", event = "event", x = "sex"))
  expect_ref(res$test_meta$log_rank$chisq, "lung_km", "log_rank_chisq")
  expect_ref(res$test_meta$log_rank$p_value, "lung_km", "log_rank_p")
})

test_that("lung: the Cox hazard ratios match survival::coxph", {
  bench_skip()
  # survminer::ggforest measures text on the open device; the test mount is
  # read-only, so a null device stands in for Rplots.pdf.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_cox_forest(bench_data("lung"),
                           bench_spec("cox_forest", time = "time", event = "event",
                                      covariates = "age,sex,ph.ecog"))
  for (term in c("age", "sexmale", "ph.ecog")) {
    r <- res$stats[res$stats$term == term, ]
    expect_ref(r$HR, "lung_cox", paste0("HR_", term))
    expect_ref(r$p_value, "lung_cox", paste0("p_", term))
  }
})

test_that("mtcars: the Pearson correlation and the correlation matrix match stats::cor", {
  bench_skip()
  res <- recipe_correlation(bench_data("mtcars"),
                            bench_spec("correlation", x = "wt", y = "mpg", test = "pearson"))
  expect_ref(res$stats$estimate, "mtcars_correlation", "estimate")
  expect_ref(res$stats$p_value, "mtcars_correlation", "p_value")
  hm <- recipe_correlation_heatmap(bench_data("mtcars"), bench_spec("correlation_heatmap"))
  expect_ref(var_pair_row(hm$stats, "mpg", "wt")$r, "mtcars_correlation_heatmap", "r_mpg_wt")
  expect_ref(var_pair_row(hm$stats, "disp", "cyl")$r, "mtcars_correlation_heatmap", "r_disp_cyl")
})

test_that("mtcars: the per-feature Welch tests of the heatmap match stats::t.test", {
  bench_skip()
  spec <- bench_spec("heatmap", annotation = file.path(bench_dir(), "data", "mtcars_annotation.csv"))
  res <- recipe_heatmap(bench_data("mtcars_matrix"), spec)
  for (f in c("mpg", "disp", "hp", "drat", "wt", "qsec")) {
    expect_ref(res$stats$p_value[res$stats$feature == f], "mtcars_heatmap", paste0("p_", f))
  }
})

test_that("airquality: the Spearman correlation drops the missing ozone values", {
  bench_skip()
  res <- recipe_correlation(bench_data("airquality"),
                            bench_spec("correlation", x = "Temp", y = "Ozone", test = "spearman"))
  expect_ref(res$stats$estimate, "airquality_spearman", "estimate")
  expect_ref(res$stats$n, "airquality_spearman", "n")
})

test_that("iris: the PCA variance matches stats::prcomp", {
  bench_skip()
  res <- recipe_pca(bench_data("iris"), bench_spec("pca", group = "Species"))
  expect_ref(res$test_meta$pc1_variance_pct, "iris_pca", "pc1_variance_pct")
  expect_ref(res$test_meta$pc2_variance_pct, "iris_pca", "pc2_variance_pct")
})

test_that("Titanic and UCBAdmissions: the chi-squared tests match stats::chisq.test", {
  bench_skip()
  ti <- recipe_proportions(bench_data("titanic"),
                           bench_spec("proportions", x = "Class", y = "Survived", test = "chisq"))
  expect_ref(ti$test_meta$statistic, "titanic_proportions", "statistic")
  expect_ref(ti$test_meta$p_value, "titanic_proportions", "p_value")
  ucb <- recipe_proportions(bench_data("ucb_admissions"),
                            bench_spec("proportions", x = "Gender", y = "Admit", test = "chisq"))
  expect_ref(ucb$test_meta$statistic, "ucb_proportions", "statistic")
  expect_ref(ucb$test_meta$p_value, "ucb_proportions", "p_value")
  expect_ref(ucb$test_meta$odds_ratio, "ucb_proportions", "odds_ratio")
})

test_that("DNase: the 4PL fit matches the self-starting SSfpl model", {
  bench_skip()
  res <- recipe_dose_response(bench_data("dnase_run1"),
                              bench_spec("dose_response", x = "conc", y = "density"))
  for (q in c("bottom", "top", "ec50", "hill", "rss")) expect_ref(res$stats[[q]], "dnase_4pl", q)
  sh <- recipe_dose_response(bench_data("dnase"),
                             bench_spec("dose_response", x = "conc", y = "density", group = "Run"))
  expect_ref(sh$test_meta$shared_ec50_test$statistic, "dnase_shared", "statistic")
  expect_ref(sh$test_meta$shared_ec50_test$p_value, "dnase_shared", "p_value")
})

test_that("aSAH: the AUCs, the DeLong CI and the DeLong test match pROC", {
  bench_skip()
  res <- recipe_roc_curve(bench_data("asah"),
                          bench_spec("roc_curve", y = "outcome", x = "s100b,ndka"))
  s <- res$stats[res$stats$marker == "s100b", ]
  n <- res$stats[res$stats$marker == "ndka", ]
  expect_ref(s$auc, "asah_roc", "auc_s100b")
  expect_ref(s$auc_ci_low, "asah_roc", "auc_low_s100b")
  expect_ref(s$auc_ci_high, "asah_roc", "auc_high_s100b")
  expect_ref(n$auc, "asah_roc", "auc_ndka")
  expect_ref(n$delong_p_vs_first, "asah_roc", "delong_p")
})

test_that("PEFR (Bland and Altman 1986): the bias and the SD match the paper", {
  bench_skip()
  res <- recipe_bland_altman(bench_data("pefr"),
                             bench_spec("bland_altman", x = "wright_1", y = "mini_1"))
  expect_ref(res$stats$bias, "pefr_bland_altman", "bias_published")
  expect_ref(res$stats$sd_diff, "pefr_bland_altman", "sd_published")
  expect_ref(res$stats$bias, "pefr_bland_altman", "bias")
  expect_ref(res$stats$loa_low, "pefr_bland_altman", "loa_low")
  expect_ref(res$stats$loa_high, "pefr_bland_altman", "loa_high")
})

test_that("BCG trials: the pooled estimates and the heterogeneity match metafor", {
  bench_skip()
  res <- recipe_forest_meta(bench_data("bcg"),
                            bench_spec("forest_meta", y = "risk_ratio", ci_low = "ci_low",
                                       ci_high = "ci_high", label = "trial"))
  expect_equal(res$test_meta$measure, "ratio")
  re <- res$stats[res$stats$label == "Random effects (DerSimonian-Laird)", ]
  fe <- res$stats[res$stats$label == "Fixed effect", ]
  expect_ref(re$estimate, "bcg_meta", "random_estimate")
  expect_ref(re$ci_low, "bcg_meta", "random_ci_low")
  expect_ref(re$ci_high, "bcg_meta", "random_ci_high")
  expect_ref(re$tau2, "bcg_meta", "tau2")
  expect_ref(re$i2, "bcg_meta", "i2")
  expect_ref(re$q, "bcg_meta", "q")
  expect_ref(fe$estimate, "bcg_meta", "fixed_estimate")
})

check_mixed <- function(res, case, group, time) {
  g <- res$stats[res$stats$effect == group, ]
  t <- res$stats[res$stats$effect == time, ]
  i <- res$stats[res$stats$effect == paste0(group, ":", time), ]
  expect_ref(g$F, case, "F_group")
  expect_ref(g$p_value, case, "p_group")
  expect_ref(t$F, case, "F_time")
  expect_ref(i$F, case, "F_interaction")
  expect_ref(i$p_value, case, "p_interaction")
  expect_ref(t$gg_epsilon, case, "gg_epsilon")
  expect_ref(res$test_meta$sphericity$p_value, case, "mauchly_p")
}

test_that("Orthodont: the mixed ANOVA matches aov and car::Anova", {
  bench_skip()
  res <- recipe_timecourse(bench_data("orthodont"),
                           bench_spec("timecourse", x = "age", y = "distance", id = "Subject", group = "Sex"))
  check_mixed(res, "orthodont_timecourse", "Sex", "age")
  expect_false(res$test_meta$correction_applied)
})

test_that("ChickWeight: the mixed ANOVA drops the incomplete chicks and corrects sphericity", {
  bench_skip()
  res <- recipe_timecourse(bench_data("chick_weight"),
                           bench_spec("timecourse", x = "Time", y = "weight", id = "Chick", group = "Diet"))
  expect_ref(res$test_meta$n_subjects, "chick_weight_timecourse", "n_complete_chicks")
  check_mixed(res, "chick_weight_timecourse", "Diet", "Time")
  expect_true(res$test_meta$correction_applied)
})

test_that("aml (Miller 1997): the swimmer arm summaries match the data", {
  bench_skip()
  res <- recipe_swimmer(bench_data("aml"),
                        bench_spec("swimmer", id = "patient", time = "time", event = "status",
                                   group = "x"))
  expect_ref(res$test_meta$arms$Maintained$median_duration, "aml_swimmer", "median_maintained")
  expect_ref(res$test_meta$arms$Nonmaintained$median_duration, "aml_swimmer", "median_nonmaintained")
  expect_ref(res$test_meta$arms$Maintained$ongoing, "aml_swimmer", "ongoing_maintained")
  expect_ref(res$test_meta$arms$Nonmaintained$ongoing, "aml_swimmer", "ongoing_nonmaintained")
})

test_that("airway (Himes 2014): the volcano counts match the DESeq2 results", {
  bench_skip()
  res <- recipe_volcano(bench_data("airway_de"), bench_spec("volcano", y = "padj"))
  for (q in c("n_up", "n_down", "n_total")) expect_ref(res$test_meta[[q]], "airway_volcano", q)
})

test_that("airway: the per-gene Welch tests of the heatmap match stats::t.test", {
  bench_skip()
  spec <- bench_spec("heatmap", annotation = file.path(bench_dir(), "data", "airway_annotation.csv"))
  res <- recipe_heatmap(bench_data("airway_heatmap"), spec)
  for (g in c("SPARCL1", "CACNB2", "DUSP1", "SAMHD1", "MAOA")) {
    expect_ref(res$stats$p_value[res$stats$feature == g], "airway_heatmap", paste0("p_", g))
  }
})

test_that("airway: the PCA variance matches stats::prcomp and the PERMANOVA is exact", {
  bench_skip()
  res <- recipe_pca(bench_data("airway_pca"), bench_spec("pca", group = "dex"))
  expect_ref(res$test_meta$pc1_variance_pct, "airway_pca", "pc1_variance_pct")
  expect_ref(res$test_meta$pc2_variance_pct, "airway_pca", "pc2_variance_pct")
  # Eight samples in two groups of four give 70 labelings, so the test is exact.
  # The observed split and its mirror are the two most extreme labelings: p = 2/70.
  expect_equal(res$test_meta$permanova$method, "exact")
  expect_equal(res$test_meta$permanova$p_value, 2 / 70)
})

test_that("airway: the top GO term of the enrichment dot plot matches phyper", {
  bench_skip()
  res <- recipe_enrichment_dot(bench_data("airway_enrichment"), bench_spec("enrichment_dot"))
  expect_ref(res$stats$count[1], "airway_enrichment", "top_count")
  expect_ref(res$stats$p_adjust[1], "airway_enrichment", "top_p_adjust")
  expect_ref(res$stats$gene_ratio[1], "airway_enrichment", "top_gene_ratio")
})

test_that("airway: the UpSet intersections match the membership patterns", {
  bench_skip()
  res <- recipe_upset(bench_data("airway_upset"), bench_spec("upset"))
  for (q in c("n_elements", "largest_intersection", "n_intersections")) {
    expect_ref(res$test_meta[[q]], "airway_upset", q)
  }
})

test_that("anorexia (Hand 1994): the family therapy paired t-test matches stats::t.test", {
  bench_skip()
  res <- recipe_paired_compare(bench_data("anorexia_ft"),
                               bench_spec("paired_compare", x = "condition", y = "weight_lbs",
                                          id = "subject", test = "t"))
  expect_ref(res$stats$statistic, "anorexia_paired", "statistic")
  expect_ref(res$stats$p_value, "anorexia_paired", "p_value")
})

test_that("veteran: the log-rank test by cell type matches survival::survdiff", {
  bench_skip()
  res <- recipe_survival_km(bench_data("veteran"),
                            bench_spec("survival_km", time = "time", event = "status", x = "celltype"))
  expect_ref(res$test_meta$log_rank$chisq, "veteran_km", "log_rank_chisq")
  expect_ref(res$test_meta$log_rank$p_value, "veteran_km", "log_rank_p")
})

test_that("veteran: the Cox hazard ratios match survival::coxph", {
  bench_skip()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_cox_forest(bench_data("veteran"),
                           bench_spec("cox_forest", time = "time", event = "status",
                                      covariates = "trt,celltype,karno,age"))
  for (term in c("trttest", "celltypelarge", "celltypesmallcell", "celltypesquamous", "karno", "age")) {
    r <- res$stats[res$stats$term == term, ]
    expect_ref(r$HR, "veteran_cox", paste0("HR_", term))
    expect_ref(r$p_value, "veteran_cox", paste0("p_", term))
  }
})

test_that("state.x77: the correlation matrix keeps the column names with a space", {
  bench_skip()
  res <- recipe_correlation_heatmap(bench_data("state_x77"), bench_spec("correlation_heatmap"))
  lm <- var_pair_row(res$stats, "Murder", "Life Exp")
  expect_ref(lm$r, "state_correlation_heatmap", "r_murder_life_exp")
  expect_ref(lm$p_adj, "state_correlation_heatmap", "padj_murder_life_exp")
  expect_ref(var_pair_row(res$stats, "Illiteracy", "Frost")$r, "state_correlation_heatmap",
             "r_illiteracy_frost")
})

test_that("mgus2 (Kyle 2002): Gray's test and the cumulative incidence match cmprsk", {
  bench_skip()
  # ggpubr::ggarrange measures the plots on the open device; the test mount is
  # read-only, so a null device stands in for Rplots.pdf.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- recipe_cumulative_incidence(bench_data("mgus2"),
                                     bench_spec("cumulative_incidence", time = "etime",
                                                event = "event", x = "sex"))
  for (cause in c("1", "2")) {
    r <- res$stats[res$stats$cause == cause, ]
    expect_ref(r$gray_statistic[1], "mgus2_cif", paste0("gray_stat_", cause))
    expect_ref(r$gray_p[1], "mgus2_cif", paste0("gray_p_", cause))
    for (g in c("F", "M")) {
      cv <- res$curves[res$curves$cause == cause & res$curves$group == g, ]
      for (tt in c(120, 240)) {
        expect_ref(.cif_value_at(tt, cv$time, cv$est), "mgus2_cif",
                   sprintf("cif_%s_%s_%d", g, cause, tt))
      }
    }
  }
})

test_that("Pima.te (Smith 1988): the calibration statistics match hoslem.test and glm", {
  bench_skip()
  res <- recipe_calibration(bench_data("pima_te"),
                            bench_spec("calibration", y = "type", x = "predicted", positive = "Yes"))
  for (q in c("hl_statistic", "hl_df", "hl_p", "brier", "slope", "intercept")) {
    expect_ref(res$stats[[q]], "pima_calibration", q)
  }
})

test_that("pbcseq (Murtaugh 1994): Bowker's test of edema matches stats::mcnemar.test", {
  bench_skip()
  res <- recipe_sankey(bench_data("pbcseq_edema"),
                       bench_spec("sankey", id = "id", x = "year", y = "edema"))
  expect_equal(res$test_meta$name, "Bowker's test of symmetry")
  for (q in c("statistic", "df", "p_value", "n_paired")) {
    expect_ref(res$test_meta[[q]], "pbcseq_sankey", q)
  }
  st <- res$stats
  expect_ref(st$n[st$from_state == "0" & st$to_state == "0.5"], "pbcseq_sankey", "n_0_to_0.5")
  expect_ref(st$n[st$from_state == "0.5" & st$to_state == "0"], "pbcseq_sankey", "n_0.5_to_0")
})
