# Behaviour of the six recipes added in v1.0.0 on small synthetic inputs. The
# numbers themselves are checked against public benchmarks in test-benchmarks.R.

new_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

fpl_value <- function(dose, bottom, top, ec50, hill) {
  bottom + (top - bottom) / (1 + 10^((log10(ec50) - log10(dose)) * hill))
}

test_that("dose_response reports an IC50 for a falling curve and drops the zero dose", {
  set.seed(1)
  dose <- rep(c(0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10), each = 3)
  df <- data.frame(dose = dose,
                   response = fpl_value(pmax(dose, 1e-6), 5, 100, 0.3, -1.2) + stats::rnorm(24, 0, 1))
  res <- recipe_dose_response(df, new_spec("dose_response", x = "dose", y = "response"))
  expect_equal(res$stats$midpoint, "IC50")
  expect_lt(res$stats$hill, 0)
  expect_equal(res$stats$ec50, 0.3, tolerance = 0.1)
  expect_true(any(grepl("zero or below", res$clean_steps)))
  expect_equal(res$stats$n, 21)
})

test_that("dose_response needs four doses per curve", {
  df <- data.frame(dose = c(1, 2, 3, 1, 2, 3), response = c(1, 2, 3, 1, 2, 3))
  expect_error(recipe_dose_response(df, new_spec("dose_response", x = "dose", y = "response")),
               "at least 4 distinct doses")
})

test_that("roc_curve reads a marker whose low values mark the case", {
  df <- data.frame(status = rep(c("case", "control"), each = 20),
                   marker = c(1:20, 11:30))
  res <- recipe_roc_curve(df, new_spec("roc_curve", y = "status", x = "marker", positive = "case"))
  expect_equal(res$stats$direction, ">")
  expect_gt(res$stats$auc, 0.5)
  expect_equal(res$test_meta$case_level, "case")
})

test_that("roc_curve takes the second sorted level as the case by default", {
  df <- data.frame(y = rep(c(0, 1), each = 15), m = c(1:15, 6:20))
  res <- recipe_roc_curve(df, new_spec("roc_curve", y = "y", x = "m"))
  expect_equal(res$test_meta$case_level, "1")
  expect_equal(res$stats$direction, "<")
})

test_that("roc_curve needs exactly two outcome levels", {
  df <- data.frame(y = rep(c("a", "b", "c"), each = 5), m = 1:15)
  expect_error(recipe_roc_curve(df, new_spec("roc_curve", y = "y", x = "m")),
               "exactly 2 outcome levels")
})

test_that("bland_altman recovers a constant offset as the bias", {
  set.seed(2)
  a <- stats::runif(30, 10, 100)
  df <- data.frame(a = a + 2 + stats::rnorm(30, 0, 0.5), b = a)
  res <- recipe_bland_altman(df, new_spec("bland_altman", x = "a", y = "b"))
  expect_equal(res$stats$bias, 2, tolerance = 0.1)
  expect_lt(res$stats$loa_low, res$stats$bias)
  expect_gt(res$stats$loa_high, res$stats$bias)
  expect_named(res$qc$panels, "a - b")
})

test_that("bland_altman gives the y axis enough breaks", {
  set.seed(4)
  a <- stats::runif(40, 40, 220)
  df <- data.frame(a = a + 6 + stats::rnorm(40, 0, 4.7), b = a)
  res <- recipe_bland_altman(df, new_spec("bland_altman", x = "a", y = "b"))
  built <- ggplot2::ggplot_build(res$plot)
  brk <- built$layout$panel_params[[1]]$y$breaks
  expect_gte(sum(is.finite(brk)), 4)
})

test_that("forest_meta pools a difference measure from standard errors by inverse variance", {
  df <- data.frame(study = c("s1", "s2", "s3"), md = c(1.0, 2.0, 1.5), se = c(0.5, 0.25, 1.0))
  res <- recipe_forest_meta(df, new_spec("forest_meta", y = "md", se = "se", label = "study"))
  w <- 1 / df$se^2
  fe <- res$stats[res$stats$label == "Fixed effect", ]
  expect_equal(res$test_meta$measure, "difference")
  expect_equal(fe$estimate, sum(w * df$md) / sum(w), tolerance = 1e-10)
  expect_equal(sum(res$stats$weight_fixed_pct[res$stats$row_type == "study"]), 100)
})

test_that("forest_meta picks the ratio scale for CIs that are symmetric on the log scale", {
  z <- stats::qnorm(0.975)
  lr <- log(c(0.6, 0.8, 0.7, 0.5)); se <- c(0.2, 0.3, 0.25, 0.4)
  df <- data.frame(or = exp(lr), lo = exp(lr - z * se), hi = exp(lr + z * se))
  res <- recipe_forest_meta(df, new_spec("forest_meta", y = "or", ci_low = "lo", ci_high = "hi"))
  expect_equal(res$test_meta$measure, "ratio")
  expect_equal(res$stats$se[res$stats$row_type == "study"], se, tolerance = 1e-10)
})

test_that("timecourse without a group runs a one-way repeated-measures ANOVA", {
  set.seed(3)
  df <- expand.grid(id = paste0("s", 1:8), time = c(0, 1, 2, 3))
  df$value <- 10 + 2 * df$time + stats::rnorm(nrow(df))
  res <- recipe_timecourse(df, new_spec("timecourse", x = "time", y = "value", id = "id"))
  expect_equal(res$test_meta$name, "one-way repeated-measures ANOVA")
  expect_equal(nrow(res$stats), 1)
  expect_lt(res$stats$p_value, 0.001)
})

test_that("timecourse drops an incomplete subject and honours --test gg", {
  set.seed(4)
  df <- expand.grid(id = paste0("s", 1:10), time = c(1, 2, 3))
  df$arm <- ifelse(df$id %in% paste0("s", 1:5), "a", "b")
  df$value <- stats::rnorm(nrow(df))
  df <- df[!(df$id == "s1" & df$time == 3), ]
  res <- recipe_timecourse(df, new_spec("timecourse", x = "time", y = "value", id = "id",
                                        group = "arm", test = "gg"))
  expect_equal(res$test_meta$n_subjects, 9)
  expect_true(res$test_meta$correction_applied)
  expect_equal(res$stats$correction[res$stats$effect == "time"], "Greenhouse-Geisser")
})

test_that("timecourse refuses two values for one subject at one time", {
  df <- data.frame(id = c("a", "a", "a", "b", "b"), time = c(1, 1, 2, 1, 2), value = 1:5)
  expect_error(recipe_timecourse(df, new_spec("timecourse", x = "time", y = "value", id = "id")),
               "one value per subject and time")
})

test_that("timecourse refuses a subject in two groups", {
  df <- data.frame(id = c("a", "a", "b", "b"), time = c(1, 2, 1, 2), value = 1:4,
                   arm = c("x", "y", "x", "x"))
  expect_error(recipe_timecourse(df, new_spec("timecourse", x = "time", y = "value", id = "id",
                                              group = "arm")),
               "one group only")
})

test_that("waterfall reads one precomputed change per patient and sorts the bars", {
  df <- data.frame(patient = c("p1", "p2", "p3", "p4"), best = c(-45, 25, -10, 5))
  res <- recipe_waterfall(df, new_spec("waterfall", y = "best", id = "patient"))
  expect_equal(res$stats$patient, c("p2", "p4", "p3", "p1"))
  expect_equal(res$stats$category_from_change,
               c("progressive disease", "stable disease", "stable disease", "partial response"))
})

test_that("waterfall without --x refuses two rows for one patient", {
  df <- data.frame(patient = c("p1", "p1"), best = c(-45, 25))
  expect_error(recipe_waterfall(df, new_spec("waterfall", y = "best", id = "patient")),
               "one row per patient")
})

test_that("swimmer sorts the courses, counts the ongoing ones and the markers", {
  df <- data.frame(pt = c("a", "b", "c"), months = c(3, 10, 6), ended = c(1, 0, 1),
                   response = c(NA, 2, 4))
  res <- recipe_swimmer(df, new_spec("swimmer", id = "pt", time = "months", event = "ended",
                                     markers = "response"))
  expect_equal(res$stats$patient, c("b", "c", "a"))
  expect_equal(res$test_meta$n_ongoing, 1)
  expect_equal(res$test_meta$markers$response, 2)
})

test_that("swimmer refuses a status that is not coded 1 or 0", {
  df <- data.frame(pt = c("a", "b"), months = c(3, 10), ended = c(1, 2))
  expect_error(recipe_swimmer(df, new_spec("swimmer", id = "pt", time = "months", event = "ended")),
               "1 = the course ended")
})
