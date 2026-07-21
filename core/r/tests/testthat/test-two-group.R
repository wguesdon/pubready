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
