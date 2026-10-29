# Behaviour of the cumulative_incidence, calibration and sankey recipes on small
# synthetic inputs. The numbers are also checked against public benchmarks in
# test-benchmarks.R. The Python suite checks the same inputs against the values
# that these R recipes report.

clin_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

cif_tiny <- function() {
  data.frame(
    time = c(2, 3, 3, 5, 6, 6, 8, 9, 11, 12, 12, 14, 1, 3, 4, 4, 6, 7, 9, 9, 10, 13, 15, 15,
             2, 2, 5, 7, 8, 8, 10, 11, 12, 14, 14, 16),
    status = c(1, 0, 2, 1, 1, 0, 2, 1, 0, 1, 2, 0, 2, 1, 1, 2, 0, 1, 1, 0, 2, 1, 0, 1,
               1, 1, 0, 2, 1, 2, 0, 1, 1, 2, 0, 0),
    group = rep(c("A", "B", "C"), each = 12))
}

test_that("cumulative_incidence matches cmprsk::cuminc on tied times with three groups", {
  d <- cif_tiny()
  res <- recipe_cumulative_incidence(d, clin_spec("cumulative_incidence", time = "time",
                                                  event = "status", x = "group"))
  fit <- cmprsk::cuminc(d$time, d$status, d$group)
  r1 <- res$stats[res$stats$cause == "1", ]
  expect_equal(r1$gray_statistic[1], unname(fit$Tests["1", "stat"]), tolerance = 1e-10)
  expect_equal(res$stats$gray_df[1], 2)
  tp <- cmprsk::timepoints(fit, c(5, 7.5, 12.5))$est
  expect_equal(r1$cif_7.5, unname(tp[c("A 1", "B 1", "C 1"), "7.5"]), tolerance = 1e-12)
  expect_equal(names(res$stats)[5:7], c("cif_5", "cif_7.5", "cif_12.5"))
  expect_equal(res$test_meta$cause_of_interest, "1")
})

test_that("cumulative_incidence puts the cause named by --cause first", {
  res <- recipe_cumulative_incidence(cif_tiny(), clin_spec("cumulative_incidence", time = "time",
                                                           event = "status", x = "group",
                                                           cause = "2"))
  expect_equal(res$stats$cause[1], "2")
  expect_equal(unlist(res$test_meta$causes), c("2", "1"))
  expect_error(recipe_cumulative_incidence(cif_tiny(), clin_spec("cumulative_incidence",
                                                                 time = "time", event = "status",
                                                                 cause = "7")),
               "is not an event code")
})

test_that("cumulative_incidence with one cause and one group equals 1 - Kaplan-Meier", {
  d <- data.frame(time = c(1, 2, 2, 3, 5, 6, 7, 9), status = c(1, 0, 1, 1, 0, 1, 0, 1))
  res <- recipe_cumulative_incidence(d, clin_spec("cumulative_incidence", time = "time",
                                                  event = "status"))
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)
  cv <- res$curves
  expect_true(is.na(res$stats$gray_p))
  expect_equal(.cif_value_at(6, cv$time, cv$est), 1 - summary(km, times = 6)$surv, tolerance = 1e-12)
  expect_true(is.na(.cif_value_at(10, cv$time, cv$est)))
})

test_that("calibration groups tied predictions as ResourceSelection::hoslem.test does", {
  set.seed(11)
  p <- round(stats::runif(200, 0.05, 0.95), 2)
  y <- stats::rbinom(200, 1, p)
  res <- recipe_calibration(data.frame(y = y, p = p), clin_spec("calibration", y = "y", x = "p"))
  hl <- ResourceSelection::hoslem.test(y, p, g = 10)
  expect_equal(res$stats$hl_statistic, unname(hl$statistic), tolerance = 1e-12)
  expect_equal(res$stats$hl_df, unname(hl$parameter))
  expect_equal(res$stats$hl_p, hl$p.value, tolerance = 1e-10)
  expect_equal(res$stats$brier, mean((p - y)^2))
  expect_equal(res$stats$case_level, "1")
})

test_that("calibration leaves out an empty decile group, where hoslem.test returns NaN", {
  # The type-7 deciles of these tied values give the break 0.63 between 0.6 and
  # 0.7, so the group (0.6, 0.63] is empty. hoslem.test divides 0 by 0 there;
  # the recipe keeps the 8 groups that hold subjects, so df = 6. The Python
  # suite checks the same value.
  p <- rep(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), c(3, 4, 3, 4, 3, 4, 2, 3, 4))
  y <- c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1)
  res <- recipe_calibration(data.frame(y = y, p = p), clin_spec("calibration", y = "y", x = "p"))
  hl <- suppressWarnings(ResourceSelection::hoslem.test(y, p, g = 10))
  expect_true(is.nan(unname(hl$statistic)))
  expect_equal(res$stats$hl_groups, 8)
  expect_equal(res$stats$hl_df, 6)
  expect_equal(res$stats$hl_statistic, 3.812047479810638, tolerance = 1e-12)
})

test_that("calibration recovers a slope of 1 and an intercept of 0 for true risks", {
  set.seed(7)
  p <- stats::plogis(stats::rnorm(4000, -0.5, 1.2))
  y <- stats::rbinom(4000, 1, p)
  d <- data.frame(outcome = ifelse(y == 1, "event", "none"), risk = p)
  res <- recipe_calibration(d, clin_spec("calibration", y = "outcome", x = "risk", positive = "event"))
  expect_lt(res$stats$slope_ci_low, 1)
  expect_gt(res$stats$slope_ci_high, 1)
  expect_lt(res$stats$intercept_ci_low, 0)
  expect_gt(res$stats$intercept_ci_high, 0)
})

test_that("calibration rejects values outside [0, 1]", {
  d <- data.frame(y = c(0, 1, 0, 1), p = c(0.2, 1.3, 0.1, 0.8))
  expect_error(recipe_calibration(d, clin_spec("calibration", y = "y", x = "p")),
               "outside \\[0, 1\\]")
})

test_that("sankey counts the transitions of the patients present at both visits", {
  d <- data.frame(pid = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4),
                  visit = c(8, 16, 24, 8, 24, 8, 16, 24, 8, 16),
                  state = c("SD", "PR", "PR", "SD", "PD", "PR", "PR", "CR", "SD", "SD"))
  res <- recipe_sankey(d, clin_spec("sankey", id = "pid", x = "visit", y = "state"))
  expect_equal(unlist(res$test_meta$visits), c("8", "16", "24"))
  expect_equal(unlist(res$test_meta$states), c("SD", "PR", "PD", "CR"))
  expect_equal(sum(res$stats$n[res$stats$from_visit == "8"]), 3)
  expect_equal(res$test_meta$name, "Sankey diagram (descriptive)")
  expect_equal(res$test_meta$no_test_reason, "more than two visits")
})

test_that("sankey runs Bowker's test and drops the empty pairs of cells", {
  from <- c(rep("a", 10), rep("b", 8), rep("c", 6))
  to <- c(rep("a", 5), rep("b", 4), "c", rep("b", 6), rep("a", 2), rep("c", 6))
  d <- data.frame(id = rep(seq_along(from), 2), visit = rep(c("pre", "post"), each = length(from)),
                  state = c(from, to))
  res <- recipe_sankey(d, clin_spec("sankey", id = "id", x = "visit", y = "state"))
  expect_equal(unlist(res$test_meta$visits), c("post", "pre"))
  expect_equal(res$test_meta$name, "Bowker's test of symmetry")
  # The pairs a-b and a-c hold patients who changed state; no patient moved
  # between b and c in either direction, so that pair is empty.
  expect_equal(res$test_meta$empty_pairs, 1)
  expect_equal(res$test_meta$df, 2)
})

test_that("sankey uses McNemar's test for two states and matches stats::mcnemar.test", {
  from <- c(rep("no", 12), rep("yes", 8))
  to <- c(rep("no", 7), rep("yes", 5), rep("no", 2), rep("yes", 6))
  d <- data.frame(id = rep(1:20, 2), visit = rep(c(0, 1), each = 20), state = c(from, to))
  res <- recipe_sankey(d, clin_spec("sankey", id = "id", x = "visit", y = "state"))
  mt <- stats::mcnemar.test(table(from, to), correct = FALSE)
  expect_equal(res$test_meta$name, "McNemar's test")
  expect_equal(res$test_meta$statistic, unname(mt$statistic), tolerance = 1e-12)
  expect_equal(res$test_meta$p_value, mt$p.value, tolerance = 1e-12)
})

sk_visits <- function() {
  data.frame(pid = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4),
             visit = c(8, 16, 24, 8, 24, 8, 16, 24, 8, 16),
             state = c("SD", "PR", "PR", "SD", "PD", "PR", "PR", "CR", "SD", "SD"))
}

test_that("sankey --levels sets the state order and leaves out an unused state", {
  res <- recipe_sankey(sk_visits(), clin_spec("sankey", id = "pid", x = "visit", y = "state",
                                              levels = "CR, PR,SD,PD,NE"))
  expect_equal(unlist(res$test_meta$states), c("CR", "PR", "SD", "PD"))
  expect_match(res$methods, "states in the order set by the analyst (CR, PR, SD, PD)", fixed = TRUE)
})

test_that("sankey --levels stops when a state of the data is not listed", {
  expect_error(recipe_sankey(sk_visits(), clin_spec("sankey", id = "pid", x = "visit",
                                                    y = "state", levels = "CR,PR,SD")),
               "--levels does not list the state(s) 'PD' of the data", fixed = TRUE)
})

test_that("sankey needs one row per patient and visit", {
  d <- data.frame(id = c(1, 1, 2), visit = c(1, 1, 2), state = c("a", "b", "a"))
  expect_error(recipe_sankey(d, clin_spec("sankey", id = "id", x = "visit", y = "state")),
               "one row per patient and visit")
})
