# The correlation heatmap recipe: a base-graphics draw closure plus a pairwise
# stats table with one row per variable pair.
test_that("recipe_correlation_heatmap builds a pairwise stats table", {
  set.seed(2)
  f <- rnorm(60)
  df <- data.frame(
    a = f + rnorm(60, sd = 0.2),
    b = f + rnorm(60, sd = 0.2),
    c = -f + rnorm(60, sd = 0.2),
    d = rnorm(60)
  )
  opt <- list(recipe = "correlation_heatmap", data = "vars.csv", test = "pearson")
  spec <- spec_from_opt(opt)

  res <- recipe_correlation_heatmap(df, spec)

  expect_true(is.function(res$draw))
  expect_null(res$plot)
  expect_equal(nrow(res$stats), 6)   # 4 variables -> 6 pairs
  expect_true(max(res$stats$r) > 0.8)
  expect_true(min(res$stats$r) < -0.8)
  expect_equal(res$test_meta$n_variables, 4)
})
