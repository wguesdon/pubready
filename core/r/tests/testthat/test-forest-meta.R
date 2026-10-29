# forest_meta: the DerSimonian-Laird truncation and the log-axis ticks. The
# Python suite (test_forest_meta.py) checks the same numbers.

fm_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "forest_meta", data = "d.csv", test = "auto"),
                                  list(...)))
}

# Five studies that agree closely, so Cochran's Q is below its df.
homogeneous_studies <- function() {
  z <- stats::qnorm(0.975)
  lr <- log(c(0.80, 0.78, 0.82, 0.79, 0.81)); se <- c(0.20, 0.25, 0.30, 0.22, 0.28)
  data.frame(study = paste("s", 1:5), hr = exp(lr), lo = exp(lr - z * se), hi = exp(lr + z * se))
}

test_that("forest_meta truncates tau2 and I2 at zero when Q is below its df", {
  res <- recipe_forest_meta(homogeneous_studies(),
                            fm_spec(y = "hr", ci_low = "lo", ci_high = "hi", label = "study"))
  het <- res$test_meta$heterogeneity
  expect_lt(het$q, het$df)
  expect_equal(het$tau2, 0)
  expect_equal(het$i2, 0)
  expect_equal(res$test_meta$random$estimate, res$test_meta$fixed$estimate, tolerance = 1e-10)
  st <- res$stats[res$stats$row_type == "study", ]
  expect_true(all(st$weight_random_pct > 0))
  expect_equal(st$weight_random_pct, st$weight_fixed_pct, tolerance = 1e-10)
})

test_that("log_ticks_ref puts a tick on each side of the null line", {
  expect_equal(log_ticks_ref(0.4, 1.57, 1), c(0.5, 1, 2))
  expect_equal(log_ticks_ref(0.7, 0.95, 1), c(0.5, 1))
  expect_equal(log_ticks_ref(0.04, 18, 1), c(0.1, 1, 10))
  expect_equal(log_ticks_ref(0.15, 6, 1), c(0.2, 0.5, 1, 2, 5))
})

test_that("forest_meta names the tau2 truncation only when it occurs", {
  spec <- fm_spec(y = "hr", ci_low = "lo", ci_high = "hi", label = "study")
  homo <- recipe_forest_meta(homogeneous_studies(), spec)
  expect_match(homo$methods, "tau2 was truncated at zero", fixed = TRUE)
  z <- stats::qnorm(0.975)
  lr <- log(c(0.3, 0.9, 1.5, 0.5, 2.0)); se <- rep(0.1, 5)
  hetero <- data.frame(study = paste("s", 1:5), hr = exp(lr),
                       lo = exp(lr - z * se), hi = exp(lr + z * se))
  res <- recipe_forest_meta(hetero, spec)
  expect_gt(res$test_meta$heterogeneity$q, res$test_meta$heterogeneity$df)
  expect_false(grepl("truncated", res$methods))
})
