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

test_that("the first outcome sits on top of the stack, as in the legend", {
  df <- data.frame(arm = rep(c("control", "treated"), each = 20),
                   response = rep(c("no", "yes", "no", "yes"), c(15, 5, 5, 15)),
                   stringsAsFactors = FALSE)
  spec <- spec_from_opt(list(recipe = "proportions", data = "p.csv", x = "arm",
                             y = "response", test = "auto"))
  res <- recipe_proportions(df, spec)
  bars <- ggplot2::ggplot_build(res$plot)$data[[1]]
  first <- bars[bars$x == 1, ]
  # Outcome "no" (the first level) is the top segment of the first bar.
  expect_equal(first$ymax[first$fill == pubready_palette[1]], 1)
  expect_equal(first$ymin[first$fill == pubready_palette[2]], 0)
})
