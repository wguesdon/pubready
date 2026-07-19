# The proportions recipe: chi-squared test on a 2x2 table of arm x response.
test_that("recipe_proportions runs a chi-squared test on a 2x2 table", {
  df <- data.frame(
    arm = c(rep("control", 20), rep("treated", 20)),
    response = c(rep("no", 15), rep("yes", 5), rep("no", 5), rep("yes", 15)),
    stringsAsFactors = FALSE)
  opt <- list(recipe = "proportions", data = "p.csv", x = "arm", y = "response", test = "auto")
  spec <- spec_from_opt(opt)

  res <- recipe_proportions(df, spec)

  expect_equal(res$test_meta$name, "Pearson's chi-squared test")
  expect_true(res$test_meta$p_value < 0.01)
  expect_equal(nrow(res$stats), 4)
})
