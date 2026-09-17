# Integration check on the two-group recipe: the stats path must produce the
# forced test, the right group sizes, and a significant p on well-separated data.
test_that("recipe_two_group_compare runs the forced Welch t-test end to end", {
  df <- data.frame(
    group = rep(c("ctrl", "treat"), each = 8),
    value = c(10, 11, 12, 13, 14, 10.5, 11.5, 12.5,
              20, 21, 22, 23, 24, 20.5, 21.5, 22.5)
  )
  opt <- list(recipe = "two_group_compare", data = "toy.csv",
              x = "group", y = "value", test = "welch_t", paired = FALSE,
              p_adjust = "none")
  spec <- spec_from_opt(opt)

  res <- recipe_two_group_compare(df, spec)

  expect_named(res, c("plot", "stats", "test_meta", "resolved", "methods",
                      "build_script", "df_used", "clean_steps", "qc"),
               ignore.order = TRUE)
  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$stats$test, "Welch two-sample t-test")
  expect_equal(res$stats$n1, 8)
  expect_equal(res$stats$n2, 8)
  expect_true(res$stats$p_value < 0.01)
  expect_setequal(c(res$stats$group1, res$stats$group2), c("ctrl", "treat"))
})

test_that("recipe_two_group_compare rejects a grouping column that is not binary", {
  df <- data.frame(group = rep(c("a", "b", "c"), each = 3), value = 1:9)
  opt <- list(recipe = "two_group_compare", data = "toy.csv",
              x = "group", y = "value", test = "auto")
  spec <- spec_from_opt(opt)
  expect_error(recipe_two_group_compare(df, spec), "exactly 2 groups")
})

test_that("the Mann-Whitney branch reports a signed effect size and Levene's p", {
  set.seed(7)
  df <- data.frame(
    group = rep(c("control", "treated"), each = 12),
    value = c(rlnorm(12, 1, 0.4), rlnorm(12, 2, 0.4)),
    stringsAsFactors = FALSE
  )
  opt <- list(recipe = "two_group_compare", data = "toy.csv",
              x = "group", y = "value", test = "wilcoxon", paired = FALSE,
              p_adjust = "none")
  spec <- spec_from_opt(opt)

  out <- recipe_two_group_compare(df, spec)

  expect_equal(out$stats$test, "Mann-Whitney U test")
  expect_equal(out$stats$effect_size_name, "rank-biserial r")
  expect_false(is.na(out$stats$effect_size))
  expect_lt(out$stats$effect_size, 0)   # control ranks below treated
  expect_gte(out$stats$effect_size, -1)
  expect_equal(out$test_meta$assumptions$equal_variance$test, "Levene")
  expect_false(is.na(out$stats$var_equal_p))
})
