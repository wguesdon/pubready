# dose_response: the "not determined" rule. The Python suite
# (test_dose_response.py) runs the same data and expects the same result.

drc_spec <- function(...) {
  spec_from_opt(utils::modifyList(list(recipe = "dose_response", data = "d.csv", test = "auto"),
                                  list(...)))
}

drc_value <- function(dose, bottom, top, ec50, hill) {
  bottom + (top - bottom) / (1 + 10^((log10(ec50) - log10(dose)) * hill))
}

# Two falling curves and one flat curve, with the same fixed noise in both
# engines.
drc_three_curves <- function() {
  dose <- rep(c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30), each = 3)
  noise <- 2 * sin(seq_along(dose) * 1.7)
  data.frame(
    compound = rep(c("active A", "active B", "inactive"), each = length(dose)),
    dose = rep(dose, 3),
    response = c(drc_value(dose, 5, 100, 0.3, -1) + noise,
                 drc_value(dose, 5, 100, 3, -1) + noise,
                 100 + noise))
}

test_that("dose_response reports a flat curve as not determined and keeps the figure", {
  res <- recipe_dose_response(drc_three_curves(),
                              drc_spec(x = "dose", y = "response", group = "compound"))
  st <- res$stats
  expect_equal(st$fit_status, c("fitted", "fitted", "not determined"))
  expect_true(is.na(st$ec50[3]))
  expect_true(all(is.na(st[3, c("ec50_ci_low", "ec50_ci_high", "hill")])))
  expect_false(is.na(st$status_reason[3]))
  expect_equal(st$ec50[1:2], c(0.3, 3), tolerance = 0.1)
  expect_equal(res$test_meta$n_determined, 2L)
  expect_equal(res$test_meta$shared_ec50_test$df1, 1)
  expect_match(res$methods, "not determined")
  expect_match(res$plot$labels$subtitle, "Not determined, left out of the test: inactive")
})

test_that("the not-determined rule names the same reason as the Python engine", {
  lx <- c(-2, -1, 0, 1)
  expect_equal(.drc_status_reason(FALSE, 0, 100, 0, lx, 0.1, 10, 100), "the fit did not converge")
  expect_equal(.drc_status_reason(TRUE, 0, 5, 0, lx, 0.1, 10, 100),
               "the response span is below 10% of the data range")
  expect_equal(.drc_status_reason(TRUE, 0, 100, 1.5, lx, 0.1, 10, 100),
               "the EC50 lies outside the tested doses")
  expect_equal(.drc_status_reason(TRUE, 0, 100, 0, lx, 0, Inf, 100),
               "the 95% CI of the EC50 is not finite")
  expect_true(is.na(.drc_status_reason(TRUE, 0, 100, 0, lx, 0.1, 10, 100)))
})
