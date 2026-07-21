# The normality QC helper: moment stats on the checked quantity and a QQ plot
# with one panel per named vector.
test_that("qc_moment_stats reports skew and kurtosis and guards small n", {
  small <- qc_moment_stats(c(1, 2))
  expect_equal(small$n, 2)
  expect_true(is.na(small$skewness))

  set.seed(1)
  norm <- qc_moment_stats(qnorm(ppoints(200)))
  expect_lt(abs(norm$skewness), 0.2)          # symmetric -> near-zero skew
  expect_lt(abs(norm$kurtosis), 0.3)          # excess kurtosis near zero
  expect_true(norm$shapiro_p > 0.05)
})

test_that("qc_normality builds one QQ panel per named vector", {
  qc <- list(quantity = "group values",
             panels = list(ctrl = rnorm(30), treat = rnorm(30)))
  out <- qc_normality(qc)
  expect_s3_class(out$plot, "ggplot")
  expect_equal(nrow(out$panels), 2)
  expect_setequal(out$panels$panel, c("ctrl", "treat"))
  expect_true(all(c("shapiro_p", "skewness", "kurtosis") %in% names(out$panels)))
})

test_that("two_group_compare exposes the checked group values for QC", {
  df <- data.frame(group = rep(c("a", "b"), each = 10), value = c(rnorm(10), rnorm(10) + 3))
  spec <- spec_from_opt(list(recipe = "two_group_compare", data = "toy.csv",
                             x = "group", y = "value", test = "auto"))
  res <- recipe_two_group_compare(df, spec)
  expect_equal(res$qc$quantity, "group values")
  expect_setequal(names(res$qc$panels), c("a", "b"))
})
