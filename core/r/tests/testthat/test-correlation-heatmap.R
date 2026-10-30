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

# The variables for the area D tests: four columns on 12 rows, with fixed
# values so that both suites see the same numbers.
ch_area_d <- function() {
  data.frame(sample = sprintf("S%02d", 1:12),
             a = c(1.2, 2.3, 3.1, 4.8, 5.2, 6.9, 7.1, 8.4, 9.0, 10.2, 11.5, 12.1),
             b = c(2.0, 2.9, 3.8, 5.1, 5.0, 7.2, 6.8, 8.9, 9.4, 9.9, 12.0, 12.6),
             c = c(5.5, 3.2, 6.1, 2.2, 4.4, 1.9, 6.6, 3.3, 5.0, 2.7, 4.1, 3.9),
             d = c(9.1, 8.2, 7.7, 6.5, 6.6, 5.1, 4.9, 3.8, 3.1, 2.2, 1.7, 0.9),
             stringsAsFactors = FALSE)
}

# A column with one text cell stays a variable; the log names the cell and the
# text column left out (F42).
test_that("correlation_heatmap keeps a mostly numeric column and logs the text cell", {
  df <- ch_area_d()
  df$c <- as.character(df$c)
  df$c[4] <- "<0.1"
  res <- recipe_correlation_heatmap(df, spec_from_opt(list(recipe = "correlation_heatmap",
                                                           data = "v.csv")))
  expect_equal(res$test_meta$n_variables, 4)
  steps <- paste(res$clean_steps, collapse = " | ")
  expect_match(steps, "Coerced column 'c' to numeric; 1 value(s) became NA, for example '<0.1'.",
               fixed = TRUE)
  expect_match(steps, "Left out the text column(s) 'sample'.", fixed = TRUE)
})
