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
  expect_true(is.na(st$midpoint[3]))
  expect_true(all(is.na(st[3, c("midpoint_ci_low", "midpoint_ci_high", "hill")])))
  expect_equal(st$status_reason[3], "the mean responses vary by less than 10% of the data range")
  expect_equal(st$midpoint[1:2], c(0.3, 3), tolerance = 0.1)
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

test_that("the check before the fit finds a flat curve as the Python engine does", {
  lx <- rep(c(-2, -1, 0, 1), each = 2)
  expect_true(.drc_is_flat(lx, c(100, 102, 99, 101, 103, 101, 100, 100), 100))
  expect_false(.drc_is_flat(lx, c(100, 102, 90, 92, 50, 52, 10, 12), 100))
  # A single curve is its own reference, so it is never flat.
  y <- c(100, 101, 100, 101, 101, 102, 100, 101)
  expect_false(.drc_is_flat(lx, y, 1.5))
})

test_that("a single flat curve is not determined by the dose ANOVA, as in Python", {
  d <- data.frame(dose = rep(c(0.01, 0.1, 1, 10, 100), each = 3),
                  response = 100 + rep(c(1, -1, 0), 5))
  res <- recipe_dose_response(d, drc_spec(x = "dose", y = "response"))
  expect_equal(res$stats$fit_status, "not determined")
  expect_equal(res$stats$status_reason,
               "the responses do not differ between the doses (one-way ANOVA, P >= 0.05)")
})

test_that("the dose ANOVA matches stats::anova and needs a replicate", {
  lx <- rep(c(-2, -1, 0, 1), each = 3)
  y <- c(1, 2, 3, 2, 3, 4, 5, 6, 7, 3, 3, 4)
  ref <- stats::anova(stats::lm(y ~ factor(lx)))[["Pr(>F)"]][1]
  expect_equal(.drc_dose_anova_p(lx, y), ref, tolerance = 1e-12)
  expect_true(is.na(.drc_dose_anova_p(c(-2, -1, 0, 1), c(1, 2, 3, 4))))
})

test_that("the dose and kinetics methods say what each point shows", {
  expect_match(.drc_points_text(TRUE), "mean +/- SEM of the replicates", fixed = TRUE)
  expect_match(.drc_points_text(FALSE), "single response", fixed = TRUE)
  expect_match(.enz_points_text(TRUE), "small points show each replicate", fixed = TRUE)
  expect_match(.enz_points_text(FALSE), "single rate", fixed = TRUE)
})

# The midpoint columns say midpoint, so an IC50 is not stored as "ec50" (F123).
test_that("the dose_response stats table names the midpoint columns", {
  res <- recipe_dose_response(drc_three_curves(),
                              drc_spec(x = "dose", y = "response", group = "compound"))
  expect_true(all(c("midpoint_type", "midpoint", "log_midpoint", "se_log_midpoint",
                    "midpoint_ci_low", "midpoint_ci_high") %in% names(res$stats)))
  expect_false(any(grepl("ec50", names(res$stats))))
  expect_equal(res$stats$midpoint_type[1], "IC50")
})

# The methods text counts the rows with a dose of zero or below (F76).
test_that("the dose_response methods count the rows with a zero dose", {
  df <- drc_three_curves()
  df <- rbind(df, data.frame(compound = c("active A", "active B"), dose = 0, response = 100))
  res <- recipe_dose_response(df, drc_spec(x = "dose", y = "response", group = "compound"))
  expect_match(res$methods, "2 row(s) with a dose of zero or below were left out of the fit",
               fixed = TRUE)
})

# A response on the scale of 1e12 is divided by a power of ten for the fit, so
# nls converges, and the plateaus come back on the original scale (F80).
test_that("dose_response fits a response on the scale of 1e12", {
  df <- drc_three_curves()
  df <- df[df$compound == "active A", ]
  df$response <- df$response * 1e10
  res <- recipe_dose_response(df, drc_spec(x = "dose", y = "response"))
  expect_equal(res$stats$fit_status, "fitted")
  expect_equal(res$stats$midpoint, 0.3, tolerance = 0.1)
  expect_equal(res$stats$top, 1e12, tolerance = 0.1)
  expect_match(res$methods, "The responses were divided by 1e+11 for the fit", fixed = TRUE)
})

# A long curve name wraps in the legend and the figure grows with it (F81).
test_that("a long dose_response curve name wraps and widens the figure", {
  df <- drc_three_curves()
  df$compound[df$compound == "active A"] <- "a very long compound name that goes on and on beyond sixty characters"
  res <- recipe_dose_response(df, drc_spec(x = "dose", y = "response", group = "compound"))
  expect_gt(res$width, 6.2)
  expect_equal(.drc_width(c("cmpd_A (IC50 = 0.045)", "compound")), 6.2)
  expect_equal(.drc_width(c("a very long compound name", "that goes on")), 4.6 + 0.2 + 0.065 * 25)
})
