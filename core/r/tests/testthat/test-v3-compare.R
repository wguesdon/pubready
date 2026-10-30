# The v3 fixes of the comparison recipes. core/python/tests/test_v3_compare.py
# holds the same cases for the Python engine.

v3_spec <- function(recipe, ...) {
  spec_from_opt(c(list(recipe = recipe, data = "toy.csv"), list(...)))
}

# F43: a constant sample has one status, and a test that divides by a variance
# of zero stops with a message that names the samples.
test_that("a constant sample is not normal and has no Shapiro P", {
  s <- shapiro_safe(rep(5, 6))
  expect_true(is.na(s$p))
  expect_false(s$normal)
})

test_that("a forced t-test on two constant groups stops and names them", {
  df <- data.frame(group = rep(c("a", "b"), each = 4), value = rep(c(5, 6), each = 4))
  expect_error(recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "value",
                                                    test = "welch_t")),
               "the groups 'a' and 'b' have zero variance")
})

test_that("Welch's ANOVA on a constant group stops and names it", {
  df <- data.frame(group = rep(c("het", "ko", "wt"), each = 4),
                   value = c(5, 5, 5, 5, 1, 2, 3, 4, 2, 4, 6, 9))
  expect_error(recipe_multi_group_compare(df, v3_spec("multi_group_compare", x = "group",
                                                      y = "value", test = "welch")),
               "the group 'het' has zero variance")
})

test_that("a paired t-test on constant differences stops", {
  df <- data.frame(subject = rep(1:5, 2), condition = rep(c("pre", "post"), each = 5),
                   value = c(1:5, 1:5 + 3))
  expect_error(recipe_paired_compare(df, v3_spec("paired_compare", x = "condition", y = "value",
                                                 id = "subject", test = "t")),
               "within-subject differences have zero variance")
})

test_that("correlation stops on a constant column and names it", {
  df <- data.frame(a = c(1, 2, 3, 4, 5), b = rep(2, 5))
  expect_error(recipe_correlation(df, v3_spec("correlation", x = "a", y = "b")),
               "the column 'b' is constant")
})

test_that("feature_panel reports a constant feature as not tested", {
  df <- data.frame(group = rep(c("a", "b"), each = 6),
                   f1 = c(1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5),
                   floor = rep(0.5, 12))
  res <- recipe_feature_panel(df, v3_spec("feature_panel", x = "group"))
  row <- res$stats[res$stats$feature == "floor", ]
  expect_equal(row$note, "not tested: constant")
  expect_true(is.na(row$p_value))
  expect_equal(res$test_meta$n_constant, 1)
  expect_equal(res$test_meta$assumptions$normality$n_checked, 2)
  expect_match(res$methods, "1 feature\\(s\\) were constant")
})

# F61: --test auto keeps the parametric test when the rank test cannot reach
# P < 0.05 at the group sizes, and the methods text gives the attainable P.
test_that("rank_test_min_p gives the smallest attainable P of each rank test", {
  expect_equal(rank_test_min_p(c(3, 3), "mann_whitney"), 0.1)
  expect_equal(rank_test_min_p(5, "signed_rank"), 0.0625)
  expect_lt(rank_test_min_p(c(3, 3, 3), "kruskal"), 0.05)
  expect_true(is.na(rank_test_min_p(c(3, 0), "mann_whitney")))
})

test_that("auto keeps Welch's t-test when the rank test cannot reach P < 0.05", {
  df <- data.frame(condition = rep(c("control", "treated"), each = 3),
                   ct = c(5, 5, 6.2, 7.1, 7.4, 7.6))
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "condition", y = "ct"))
  expect_equal(res$stats$test, "Welch two-sample t-test")
  expect_equal(res$test_meta$rank_test_attainable_p, 0.1)
  expect_match(res$methods, "smallest attainable P = 0.10), so the parametric test was kept",
               fixed = TRUE)
  forced <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "condition", y = "ct",
                                                 test = "wilcoxon"))
  expect_equal(forced$stats$test, "Mann-Whitney U test")
  expect_match(forced$methods, "smallest attainable P = 0.10).", fixed = TRUE)
})

test_that("paired_compare keeps the paired t-test at 5 pairs", {
  df <- data.frame(subject = rep(1:5, 2), condition = rep(c("pre", "post"), each = 5),
                   value = c(1:5, 1:5 + c(1, 1, 1, 1.2, 5)))
  res <- recipe_paired_compare(df, v3_spec("paired_compare", x = "condition", y = "value",
                                           id = "subject"))
  expect_equal(res$stats$test, "paired t-test")
  expect_equal(res$test_meta$rank_test_attainable_p, 0.0625)
  expect_match(res$methods, "so the parametric test was kept")
})

test_that("multi_group_compare keeps the ANOVA at 2 values per group", {
  df <- data.frame(group = rep(c("a", "b", "c"), each = 2), value = c(1, 2, 5, 6, 9, 10.5))
  res <- recipe_multi_group_compare(df, v3_spec("multi_group_compare", x = "group", y = "value"))
  expect_true(res$test_meta$name %in% c("one-way ANOVA", "Welch's one-way ANOVA"))
  expect_gt(res$test_meta$rank_test_attainable_p, 0.05)
  expect_match(res$methods, "so the parametric test was kept")
})

test_that("feature_panel keeps the parametric family at 3 against 3", {
  df <- data.frame(group = rep(c("a", "b"), each = 3),
                   f1 = c(5, 5, 6.2, 7.1, 7.4, 7.6), f2 = c(1, 2, 3, 2, 3, 4.5))
  res <- recipe_feature_panel(df, v3_spec("feature_panel", x = "group"))
  expect_equal(res$test_meta$test, "welch_t")
  expect_true(res$test_meta$kept_parametric)
  expect_match(res$methods, "so the parametric test was kept")
})

# F66: the paired path of two_group_compare names the effect Cohen's dz and
# writes no variance sentence.
test_that("the paired t-test of two_group_compare reports Cohen's dz", {
  df <- data.frame(group = rep(c("a", "b"), each = 6),
                   value = c(1, 2, 3, 4, 5, 6, 2.1, 2.9, 4.4, 4.8, 6.3, 7.2))
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "value",
                                              paired = TRUE, test = "t"))
  expect_equal(res$stats$effect_size_name, "Cohen's dz")
  d <- df$value[1:6] - df$value[7:12]
  expect_equal(res$stats$effect_size, mean(d) / sd(d), tolerance = 1e-9)
  expect_false(grepl("variances", res$methods))
  expect_match(res$methods, "Cohen's dz = -4.82 (a minus b)", fixed = TRUE)
})

# F36: under --paired, --test auto checks the normality of the within-pair
# differences, as paired_compare does.
test_that("two_group_compare --paired chooses the test from the differences", {
  sleep_df <- data.frame(group = rep(c("1", "2"), each = 10),
                         extra = c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0,
                                   1.9, 0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 3.4))
  res <- recipe_two_group_compare(sleep_df, v3_spec("two_group_compare", x = "group",
                                                    y = "extra", paired = TRUE))
  expect_equal(res$stats$test, "Wilcoxon signed-rank test")
  expect_equal(res$stats$shapiro_p_diff, 0.03334161, tolerance = 1e-6)
  expect_equal(res$qc$quantity, "within-pair differences")
  expect_match(res$methods, "normality of the within-pair differences (1 minus 2)", fixed = TRUE)
  expect_error(recipe_two_group_compare(sleep_df[-1, ], v3_spec("two_group_compare", x = "group",
                                                                y = "extra", paired = TRUE)),
               "need the same size")
})

# F65: Cohen's d of two_group_compare uses the pooled SD in both engines.
test_that("two_group_compare reports Cohen's d over the pooled SD", {
  df <- data.frame(group = c(rep("F", 4), rep("M", 7)),
                   hwt = c(7, 7.4, 9.5, 7.2, 9.6, 11.2, 10.1, 12.5, 8.8, 13, 11.4))
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "hwt",
                                              test = "student_t"))
  a <- df$hwt[1:4]; b <- df$hwt[5:11]
  sp <- sqrt((3 * var(a) + 6 * var(b)) / 9)
  expect_equal(res$stats$effect_size, (mean(a) - mean(b)) / sp, tolerance = 1e-12)
  expect_match(res$methods, "Cohen's d is the mean difference over the pooled standard deviation.",
               fixed = TRUE)
})

# F54: a level whose every outcome is missing leaves no empty slot.
test_that("factorial_anova drops a level whose every outcome is missing", {
  df <- data.frame(genotype = rep(c("wt", "het", "ko"), each = 8),
                   treatment = rep(rep(c("vehicle", "drug"), each = 4), 3),
                   response = c(26.3, 19.5, 24.1, 23.5, 25.1, 28.1, 28.1, 23.5, rep(NA, 8),
                                29, 28.3, 34.4, 30.3, 40.2, 42.5, 42.1, 42.5))
  res <- recipe_factorial_anova(df, v3_spec("factorial_anova", x = "genotype", fill = "treatment",
                                            y = "response"))
  expect_equal(levels(res$df_used$genotype), c("ko", "wt"))
})

# F38, F64: one level rule in both engines, with --levels in each comparison
# recipe.
test_that("numeric groups sort as numbers in the comparison recipes", {
  df <- data.frame(conc = rep(c(160, 80), each = 5),
                   total = c(20, 22, 19, 25, 23, 30, 31, 28, 33, 29))
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "conc", y = "total",
                                              test = "welch_t"))
  expect_equal(c(res$stats$group1, res$stats$group2), c("80", "160"))
  expect_gt(res$stats$statistic, 0)
  mg <- data.frame(dose = rep(c(0, 160, 235, 310, 80), each = 4),
                   y = c(1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6))
  res <- recipe_multi_group_compare(mg, v3_spec("multi_group_compare", x = "dose", y = "y",
                                                test = "anova"))
  expect_equal(res$stats$group1[1:4], c("0", "0", "0", "0"))
  expect_equal(res$stats$group2[1:4], c("80", "160", "235", "310"))
})

test_that("text groups sort by code point and --levels sets the order", {
  df <- data.frame(group = rep(c("control", "Drug"), each = 5),
                   v = c(1, 2, 3, 4, 5, 3, 4, 5, 6, 7.5))
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "v",
                                              test = "wilcoxon"))
  expect_equal(res$stats$group1, "Drug")
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "v",
                                              test = "wilcoxon", levels = "control,Drug"))
  expect_equal(res$stats$group1, "control")
})

test_that("paired_compare takes the level rule, not the row order", {
  df <- data.frame(subject = rep(1:6, 2), condition = rep(c("pre", "post"), each = 6),
                   value = c(1:6, c(2.1, 2.9, 4.4, 4.8, 6.3, 7.2)))
  res <- recipe_paired_compare(df, v3_spec("paired_compare", x = "condition", y = "value",
                                           id = "subject", test = "t"))
  expect_equal(res$stats$condition1, "post")
  res <- recipe_paired_compare(df, v3_spec("paired_compare", x = "condition", y = "value",
                                           id = "subject", test = "t", levels = "pre,post"))
  expect_equal(res$stats$condition1, "pre")
  expect_equal(res$stats$contrast, "post - pre")
})

# F44: column names with a space run in the R recipes that build a formula.
test_that("two_group, multi_group and factorial_anova accept a name with a space", {
  df <- data.frame(`treatment group` = rep(c("a", "b"), each = 6),
                   `tumor volume` = c(1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5), check.names = FALSE)
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "treatment group",
                                              y = "tumor volume"))
  expect_equal(res$stats$outcome, "tumor volume")
  mg <- data.frame(`mouse genotype` = rep(c("het", "ko", "wt"), each = 4),
                   `expression level` = c(1, 2, 3, 4, 5, 6, 7, 8.5, 2, 3, 4, 5.5), check.names = FALSE)
  res <- recipe_multi_group_compare(mg, v3_spec("multi_group_compare", x = "mouse genotype",
                                                y = "expression level", test = "anova"))
  expect_equal(nrow(res$stats), 3)
  fa <- data.frame(`mouse genotype` = rep(c("wt", "ko"), each = 8),
                   `drug treatment` = rep(rep(c("vehicle", "drug"), each = 4), 2),
                   `response value` = c(5, 6, 7, 8, 9, 10, 11, 12.5, 4, 5, 6, 7.5, 12, 13, 15, 16),
                   check.names = FALSE)
  res <- recipe_factorial_anova(fa, v3_spec("factorial_anova", x = "mouse genotype",
                                            fill = "drug treatment", y = "response value",
                                            test = "art"))
  expect_equal(res$stats$effect, c("mouse genotype", "drug treatment",
                                   "mouse genotype:drug treatment"))
})

# Write the input and the script into a new temporary folder, run the script
# with Rscript there, and return the exit status and the output.
v3_run_script <- function(script, df) {
  dir <- tempfile("v3_script_")
  dir.create(dir)
  utils::write.csv(df, file.path(dir, "input.csv"), row.names = FALSE)
  writeLines(script, file.path(dir, "script.R"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "script.R",
                                  stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status"), out = paste(out, collapse = "\n"),
       pdf = file.exists(file.path(dir, "figure.pdf")))
}

test_that("the emitted scripts run on names with a space (F44)", {
  df <- data.frame(`treatment group` = rep(c("a", "b"), each = 6),
                   `tumor volume` = c(1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5), check.names = FALSE)
  res <- recipe_two_group_compare(df, v3_spec("two_group_compare", x = "treatment group",
                                              y = "tumor volume"))
  run <- v3_run_script(call_build_script(res$build_script, "input.csv", "figure"), df)
  expect_null(run$status, info = run$out)
  expect_true(run$pdf)
  fa <- data.frame(`mouse genotype` = rep(c("wt", "ko"), each = 8),
                   `drug treatment` = rep(rep(c("vehicle", "drug"), each = 4), 2),
                   `response value` = c(5, 6, 7, 8, 9, 10, 11, 12.5, 4, 5, 6, 7.5, 12, 13, 15, 16),
                   check.names = FALSE)
  res <- recipe_factorial_anova(fa, v3_spec("factorial_anova", x = "mouse genotype",
                                            fill = "drug treatment", y = "response value",
                                            test = "art"))
  run <- v3_run_script(call_build_script(res$build_script, "input.csv", "figure"), fa)
  expect_null(run$status, info = run$out)
  expect_true(run$pdf)
})

# F63: each recipe rejects a --test value that it does not know and lists the
# values that it knows.
test_that("an unknown --test value stops the comparison recipes", {
  pr <- data.frame(arm = rep(c("a", "b"), each = 10), response = rep(c("yes", "no"), 10))
  expect_error(recipe_proportions(pr, v3_spec("proportions", x = "arm", y = "response",
                                              test = "bogus")),
               "proportions does not know --test 'bogus'; use one of: auto, chisq")
  fa <- data.frame(g = rep(c("wt", "ko"), each = 8), t = rep(c("v", "d"), 8), y = 1:16)
  expect_error(recipe_factorial_anova(fa, v3_spec("factorial_anova", x = "g", fill = "t", y = "y",
                                                  test = "bogus")),
               "factorial_anova does not know --test 'bogus'; use one of: auto, anova, aov, art")
  expect_error(recipe_multi_group_compare(fa, v3_spec("multi_group_compare", x = "g", y = "y",
                                                      test = "welch_t")),
               "multi_group_compare does not know --test 'welch_t'")
})
