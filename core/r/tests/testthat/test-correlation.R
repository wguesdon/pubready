# The correlation recipe: the method resolver and an end-to-end run that must
# report the forced method, the right n, and a strong, significant coefficient.
test_that("resolve_correlation_method picks the method from normality", {
  expect_equal(resolve_correlation_method("auto", TRUE, TRUE)$method, "pearson")
  expect_equal(resolve_correlation_method("auto", TRUE, FALSE)$method, "spearman")
  expect_equal(resolve_correlation_method("auto", FALSE, TRUE)$method, "spearman")
  expect_equal(resolve_correlation_method("spearman", TRUE, TRUE)$method, "spearman")
  expect_error(resolve_correlation_method("nonsense", TRUE, TRUE),
               "unknown correlation method")
})

test_that("recipe_correlation runs end to end", {
  set.seed(1)
  x <- seq(0, 10, length.out = 40)
  y <- 2 * x + rnorm(40)
  df <- data.frame(x = x, y = y)
  opt <- list(recipe = "correlation", data = "toy.csv", x = "x", y = "y",
              test = "pearson")
  spec <- spec_from_opt(opt)

  res <- recipe_correlation(df, spec)

  expect_s3_class(res$plot, "ggplot")
  expect_equal(res$stats$method, "Pearson correlation")
  expect_equal(res$stats$n, 40)
  expect_true(res$stats$estimate > 0.9)
  expect_true(res$stats$p_value < 1e-6)
})
