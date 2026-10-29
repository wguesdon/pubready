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

test_that("proportions writes the omnibus test and names the odds ratio", {
  df <- data.frame(arm = c(rep("control", 20), rep("treated", 20)),
                   response = c(rep("no", 15), rep("yes", 5), rep("no", 5), rep("yes", 15)),
                   stringsAsFactors = FALSE)
  res <- recipe_proportions(df, spec_from_opt(list(recipe = "proportions", data = "p.csv",
                                                   x = "arm", y = "response", test = "chisq")))
  expect_equal(unique(res$stats$omnibus_statistic), res$test_meta$statistic)
  expect_equal(unique(res$stats$omnibus_p), res$test_meta$p_value)
  expect_match(res$methods, "(chi-square(1) = ", fixed = TRUE)
  expect_match(res$methods, "The odds ratio of no for control versus treated was 9.00.", fixed = TRUE)
})

test_that("multi_group_compare writes the omnibus test on each row", {
  spec <- spec_from_opt(list(recipe = "multi_group_compare", data = "d.csv", x = "group",
                             y = "value", test = "anova"))
  df <- data.frame(group = rep(c("a", "b", "c"), each = 6),
                   value = c(1:6, (1:6) + 2, (1:6) + 5))
  res <- recipe_multi_group_compare(df, spec)
  expect_equal(unique(res$stats$omnibus_p), res$test_meta$omnibus$p_value)
  expect_equal(unique(res$stats$omnibus_df2), 15)
})

test_that("proportions draws the groups in the order of --levels", {
  df <- data.frame(arm = rep(c("b", "c", "a"), each = 10),
                   response = rep(c("no", "yes"), 15), stringsAsFactors = FALSE)
  base <- list(recipe = "proportions", data = "p.csv", x = "arm", y = "response", test = "chisq")
  res <- recipe_proportions(df, spec_from_opt(c(base, levels = "a,b,c")))
  expect_equal(unique(as.character(res$stats$group)), c("a", "b", "c"))
  expect_equal(unique(as.character(recipe_proportions(df, spec_from_opt(base))$stats$group)),
               c("b", "c", "a"))
  expect_error(recipe_proportions(df, spec_from_opt(c(base, levels = "a,b"))), "'c'")
})
