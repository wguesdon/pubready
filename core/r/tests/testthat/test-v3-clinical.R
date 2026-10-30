# The v3 fixes of the clinical recipes. The Python suite (test_v3_clinical.py)
# runs the same cases and expects the same results.

v3c_spec <- function(recipe, ...) {
  spec_from_opt(utils::modifyList(list(recipe = recipe, data = "d.csv", test = "auto"),
                                  list(...)))
}

# Write the input and the script into a new temporary folder, run the script
# with Rscript there, and return the exit status and the output.
v3c_run_script <- function(script, df) {
  dir <- tempfile("script_")
  dir.create(dir)
  utils::write.csv(df, file.path(dir, "input.csv"), row.names = FALSE)
  writeLines(script, file.path(dir, "script.R"))
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "script.R",
                                  stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status"), out = paste(out, collapse = "\n"))
}

v3c_km_space <- function() {
  data.frame(`follow up months` = c(2, 5, 5, 8, 10, 10, 14, 18, 21.8, 3, 6, 10, 12, 15, 17, 19, 20, 20),
             `death event` = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1),
             `treatment arm` = rep(c("A", "B"), each = 9),
             check.names = FALSE)
}

test_that("survival_km reads columns whose names hold a space (F44)", {
  res <- recipe_survival_km(v3c_km_space(), v3c_spec("survival_km", time = "follow up months",
                                                     event = "death event", x = "treatment arm"))
  expect_equal(res$stats$group, c("A", "B"))
  expect_equal(res$stats$events, c(7L, 6L))
  d <- stats::setNames(v3c_km_space(), c("t", "e", "g"))
  ref <- survival::survdiff(survival::Surv(t, e) ~ g, data = d)
  expect_equal(res$test_meta$log_rank$chisq, unname(ref$chisq), tolerance = 1e-12)
  expect_match(res$methods, "groups of treatment arm")
  script <- call_build_script(res$build_script, "input.csv", "figure")
  run <- v3c_run_script(script, v3c_km_space())
  expect_true(is.null(run$status) || run$status == 0, info = run$out)
})

test_that("an event column holds 0, 1 or a logical, else the message lists its values (F56)", {
  expect_equal(event_codes(c(TRUE, FALSE, NA), "e"), c(1L, 0L, NA))
  expect_equal(event_codes(c("1", " 0", NA, ""), "e"), c(1L, 0L, NA, NA))
  expect_equal(event_codes(c(1, 0, 1, NA), "e"), c(1L, 0L, 1L, NA))
  expect_error(event_codes(c("dead", "alive", "dead"), "status"),
               "event column 'status' must be coded 1 = event, 0 = censored, but it holds 'alive', 'dead'",
               fixed = TRUE)
  expect_error(event_codes(c(1, 1.5, 0, 2, 10), "e"), "but it holds '1.5', '2', '10'", fixed = TRUE)
  expect_error(event_codes(c(0, 1, 2:9), "e"), "'6', '7' and 2 more", fixed = TRUE)
  d <- data.frame(time = c(3, 5, 8, 9), event = c("dead", "alive", "dead", "alive"), arm = c("a", "a", "b", "b"))
  expect_error(recipe_survival_km(d, v3c_spec("survival_km", time = "time", event = "event", x = "arm")),
               "holds 'alive', 'dead'")
  expect_error(recipe_cox_forest(d, v3c_spec("cox_forest", time = "time", event = "event", covariates = "arm")),
               "holds 'alive', 'dead'")
  s <- data.frame(patient = c("P1", "P2"), months = c(4, 7), ended = c("yes", "no"))
  expect_error(recipe_swimmer(s, v3c_spec("swimmer", id = "patient", time = "months", event = "ended")),
               "must be coded 1 = the course ended, 0 = ongoing, but it holds 'no', 'yes'", fixed = TRUE)
})

test_that("the median survival takes the midpoint of a step at exactly 0.5 (F67)", {
  d <- data.frame(time = c(1:10, seq(2, 20, 2), 1, 2, 2.5, 3),
                  event = c(rep(1, 20), 1, 1, 0, 1),
                  arm = rep(c("a", "b", "c"), c(10, 10, 4)))
  res <- recipe_survival_km(d, v3c_spec("survival_km", time = "time", event = "event", x = "arm"))
  expect_equal(res$stats$median_survival, c(5.5, 11, 2.5))
  ref <- summary(survival::survfit(survival::Surv(time, event) ~ arm, data = d))$table[, "median"]
  expect_equal(unname(res$stats$median_survival), unname(ref))
})

test_that("survival_km names one group and draws follow-up times above 5e6 (F92)", {
  d <- data.frame(time = c(3, 5, 8, 9), event = c(1, 0, 1, 1), arm = "control")
  expect_error(recipe_survival_km(d, v3c_spec("survival_km", time = "time", event = "event", x = "arm")),
               "column 'arm' holds 1 group ('control', n = 4)", fixed = TRUE)
  expect_equal(km_break_step(1.038e8), 2e7)
  expect_equal(km_break_step(5.5e12), 1e12)
  d2 <- data.frame(time = c(1, 3, 5, 8, 9, 10) * 1e7, event = c(1, 0, 1, 1, 1, 0),
                   arm = rep(c("a", "b"), each = 3))
  res <- recipe_survival_km(d2, v3c_spec("survival_km", time = "time", event = "event", x = "arm"))
  expect_equal(res$test_meta$at_risk$times, c(0, 2, 4, 6, 8, 10) * 1e7)
})

test_that("survival_km, cox_forest and cumulative_incidence share one level order (F64)", {
  d <- data.frame(time = c(5, 8, 12, 3, 9, 14, 7, 10, 15, 4, 6, 11),
                  event = c(1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0),
                  dose = rep(c(80, 160, 235), each = 4),
                  inherit = rep(c("X-linked", "autosomal"), 6),
                  status = c(1, 2, 0, 1, 1, 0, 2, 1, 0, 1, 2, 1))
  km <- recipe_survival_km(d, v3c_spec("survival_km", time = "time", event = "event", x = "dose"))
  expect_equal(km$stats$group, c("80", "160", "235"))
  cx <- recipe_cox_forest(d, v3c_spec("cox_forest", time = "time", event = "event",
                                      covariates = "inherit"))
  expect_equal(cx$stats$term, "inheritautosomal")
  ci <- recipe_cumulative_incidence(d, v3c_spec("cumulative_incidence", time = "time",
                                                event = "status", x = "dose"))
  expect_equal(unlist(ci$test_meta$groups), c("80", "160", "235"))
  expect_equal(unique(ci$stats$group), c("80", "160", "235"))
})

v3c_cox_data <- function() {
  data.frame(time = c(5, 8, 12, 3, 9, 14, 7, 10, 15, 4, 6, 11, 13, 2, 16, 9),
             event = c(1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1),
             arm = rep(c("control", "treated"), 8),
             age = c(61, 55, 70, 48, 66, 59, 52, 73, 45, 68, 57, 63, 50, 71, 58, 64))
}

test_that("a text cell in a numeric Cox covariate becomes NA and is logged (F42)", {
  d <- v3c_cox_data()
  d$age <- as.character(d$age)
  d$age[3] <- "45 y"
  res <- recipe_cox_forest(d, v3c_spec("cox_forest", time = "time", event = "event",
                                       covariates = "arm,age"))
  expect_equal(res$stats$term, c("armtreated", "age"))
  expect_equal(res$test_meta$n, 15L)
  expect_equal(res$clean_steps[1],
               "Coerced column 'age' to numeric; 1 value(s) became NA, for example '45 y'.")
  expect_match(res$clean_steps[2], "Dropped 1 row")
})

v3c_cox_spec <- function(covs) {
  v3c_spec("cox_forest", time = "time", event = "event", covariates = covs)
}

test_that("cox_forest names a constant or collinear covariate and a data set with no event (F91)", {
  d <- v3c_cox_data()
  d$site <- "A"
  d$dose <- 10
  d$age_copy <- 2 * d$age + 1
  d$arm2 <- d$arm
  expect_error(recipe_cox_forest(d, v3c_cox_spec("arm,age,site")),
               "covariate 'site' holds one value ('A') in the 16 subject(s) kept", fixed = TRUE)
  expect_error(recipe_cox_forest(d, v3c_cox_spec("arm,dose")), "covariate 'dose' holds one value ('10')",
               fixed = TRUE)
  expect_error(recipe_cox_forest(d, v3c_cox_spec("arm,age,age_copy")),
               "covariate 'age_copy' is collinear with the covariates before it", fixed = TRUE)
  expect_error(recipe_cox_forest(d, v3c_cox_spec("arm,age,arm2")), "covariate 'arm2' is collinear", fixed = TRUE)
  d$event <- 0
  expect_error(recipe_cox_forest(d, v3c_cox_spec("arm,age")),
               "column 'event' holds no 1 among the 16 subject(s) kept", fixed = TRUE)
  x <- cbind(a = c(1, 2, 3, 5), b = c(2, 1, 4, 4))
  expect_equal(cox_collinear_columns(cbind(x, x[, 1] + 2 * x[, 2])), 3L)
  expect_equal(cox_collinear_columns(cbind(x[, 1], 3 - x[, 1], x[, 2])), 2L)
  expect_equal(cox_collinear_columns(x), integer(0))
})

test_that("a covariate that separates the events is marked not estimable (F57)", {
  d <- v3c_cox_data()
  d$marker <- ifelse(d$event == 1, "pos", "neg")
  # coxph warns that the coefficient may be infinite; the recipe marks the term.
  res <- suppressWarnings(recipe_cox_forest(d, v3c_cox_spec("arm,age,marker")))
  expect_equal(res$stats$note, c("estimated", "estimated", "not estimable, separation"))
  expect_true(all(is.na(unlist(res$stats[3, c("HR", "ci_low", "ci_high", "z", "p_value")]))))
  expect_true(is.finite(res$stats$HR[1]))
  expect_match(res$methods, "The term markerpos was not estimable because it separates the events",
               fixed = TRUE)
  expect_match(res$methods, "held the highest value of markerpos", fixed = TRUE)
  time <- c(1, 2, 3, 4, 5, 6)
  event <- c(1, 1, 0, 1, 0, 0)
  x <- cbind(up = c(1, 1, 0, 1, 0, 0), down = c(0, 0, 1, 0, 1, 1), mixed = c(1, 0, 1, 0, 1, 0))
  expect_equal(cox_separation(time, event, x), c(1L, -1L, 0L))
  ok <- recipe_cox_forest(v3c_cox_data(), v3c_cox_spec("arm,age"))
  expect_equal(ok$stats$note, c("estimated", "estimated"))
})

test_that("cumulative_incidence stops on a text or a negative status code (F58)", {
  d <- data.frame(time = c(2, 3, 5, 6, 8, 9), arm = rep(c("a", "b"), 3),
                  status = c("relapse", "censored", "death", "relapse", "censored", "death"))
  sp <- v3c_spec("cumulative_incidence", time = "time", event = "status", x = "arm")
  expect_error(recipe_cumulative_incidence(d, sp),
               "column 'status' must hold 0 for a censored subject and a positive number for each cause, but it holds 'censored', 'death', 'relapse'",
               fixed = TRUE)
  d$status <- c(1, -1, 2, 1, -1, 2)
  expect_error(recipe_cumulative_incidence(d, sp), "but it holds '-1'", fixed = TRUE)
  d$status <- c(1, 0, 2, 1, 0, 2)
  expect_equal(unlist(recipe_cumulative_incidence(d, sp)$test_meta$causes), c("1", "2"))
})

test_that("sankey visits and timecourse time points share one time order (F39)", {
  expect_equal(time_level_order(c("week 8", "week 24", "week 16"))$levels,
               c("week 8", "week 16", "week 24"))
  expect_equal(time_level_order(c("M0", "M12", "M3"))$rule,
               "the order of the number that each visit holds")
  expect_equal(time_level_order(c("Day 14", "Day -7", "Day 0"))$levels, c("Day -7", "Day 0", "Day 14"))
  expect_equal(time_level_order(c("week-8", "week-2"))$levels, c("week-2", "week-8"))
  o <- time_level_order(c("10", "2", "1"))
  expect_equal(o$levels, c("1", "2", "10"))
  expect_equal(o$rule, "numeric order")
  o <- time_level_order(c("week 4", "baseline", "C1D8"), what = "time point")
  expect_equal(o$levels, c("week 4", "baseline", "C1D8"))
  expect_equal(o$rule, "the order of first appearance, because the time points do not each hold one number")
  expect_equal(time_level_order(c("b1", "a1", "c0"))$levels, c("c0", "a1", "b1"))
  o <- time_level_order(c("M0", "M3"), c("M3", "M0", "M6"))
  expect_equal(o$levels, c("M3", "M0"))
  expect_equal(o$rule, "the order that --visit_levels gives")
  expect_error(time_level_order(c("M0", "M3"), "M0"), "--visit_levels does not list the visit(s) 'M3'",
               fixed = TRUE)
  d <- data.frame(patient = rep(sprintf("P%d", 1:4), each = 3),
                  visit = rep(c("week 8", "week 24", "week 16"), 4),
                  response = c("PR", "CR", "CR", "SD", "PR", "PR", "PD", "PD", "SD", "PR", "PR", "CR"))
  res <- recipe_sankey(d, v3c_spec("sankey", id = "patient", x = "visit", y = "response"))
  expect_equal(unlist(res$test_meta$visits), c("week 8", "week 16", "week 24"))
  expect_true("The visits of 'visit' are in the order of the number that each visit holds: week 8, week 16, week 24." %in%
                res$clean_steps)
  tc <- data.frame(id = rep(sprintf("m%d", 1:6), each = 3), day = rep(c("day 10", "day 2", "day 1"), 6),
                   y = c(5, 3, 1, 6, 3, 2, 7, 4, 1, 6, 2, 2, 8, 4, 2, 5, 3, 1))
  res <- recipe_timecourse(tc, v3c_spec("timecourse", x = "day", y = "y", id = "id"))
  expect_true("The time points of 'day' are in the order of the number that each time point holds: day 1, day 2, day 10." %in%
                res$clean_steps)
  res <- recipe_timecourse(tc, v3c_spec("timecourse", x = "day", y = "y", id = "id",
                                        levels = "day 10,day 2,day 1"))
  expect_true("The time points of 'day' are in the order that --levels gives: day 10, day 2, day 1." %in%
                res$clean_steps)
})

test_that("timecourse stops on an outcome with zero variance (F90)", {
  tc <- data.frame(id = rep(sprintf("m%d", 1:6), each = 3), day = rep(c(0, 7, 14), 6),
                   arm = rep(c("a", "b"), each = 9), y = 100)
  expect_error(recipe_timecourse(tc, v3c_spec("timecourse", x = "day", y = "y", id = "id", group = "arm")),
               "the outcome 'y' has zero variance: every value of the 6 complete subject(s) is 100",
               fixed = TRUE)
})

test_that("Mauchly's test is not available with too few subjects, and GG applies (F95)", {
  tc <- data.frame(id = rep(sprintf("m%d", 1:4), each = 4), day = rep(c(0, 7, 14, 21), 4),
                   arm = rep(c("a", "b"), each = 8),
                   y = c(10, 14, 19, 25, 11, 13, 20, 27, 9, 12, 14, 15, 10, 11, 15, 17))
  res <- recipe_timecourse(tc, v3c_spec("timecourse", x = "day", y = "y", id = "id", group = "arm"))
  expect_false(res$test_meta$sphericity$available)
  expect_true(is.na(res$test_meta$sphericity$p_value))
  expect_true(res$test_meta$correction_applied)
  expect_equal(res$stats$correction, c("none", "Greenhouse-Geisser", "Greenhouse-Geisser"))
  expect_match(res$methods, "Mauchly's test was not available, because it needs at least 5 complete subjects (4 subject(s) in 2 group(s) at 4 time points)",
               fixed = TRUE)
})

test_that("waterfall names a text time and timecourse names a set with no complete subject (F92)", {
  wf <- data.frame(patient = rep(c("P1", "P2"), each = 2), week = rep(c("week 0", "week 8"), 2),
                   lesion = c(50, 40, 30, 36))
  expect_error(recipe_waterfall(wf, v3c_spec("waterfall", x = "week", y = "lesion", id = "patient")),
               "column 'week' holds no number (4 of 4 values are text, for example 'week 0')", fixed = TRUE)
  expect_error(recipe_spider_response(wf, v3c_spec("spider_response", x = "week", y = "lesion", id = "patient")),
               "column 'week' holds no number (4 of 4 values are text", fixed = TRUE)
  tc <- data.frame(id = c("a", "a", "b", "b", "c", "c"), day = c(0, 7, 7, 14, 0, 14), y = 1:6)
  expect_error(recipe_timecourse(tc, v3c_spec("timecourse", x = "day", y = "y", id = "id")),
               "none of the 3 subject(s) of 'id' has a value at every one of the 3 time points of 'day' (0, 7, 14)",
               fixed = TRUE)
})

test_that("waterfall and spider_response name a patient with two rows at one time (F59)", {
  d <- data.frame(patient = c("P1", "P1", "P1", "P2", "P2", "P3", "P3", "P3"),
                  week = c(0, 0, 8, 0, 8, 0, 8, 8),
                  lesion = c(50, 52, 40, 30, 36, 20, 21, 22))
  msg <- "needs one row per patient and time point, and 2 patient(s) have two rows or more at one time: 'P1' at week = 0, 'P3' at week = 8"
  expect_error(recipe_waterfall(d, v3c_spec("waterfall", x = "week", y = "lesion", id = "patient")),
               paste("waterfall", msg), fixed = TRUE)
  expect_error(recipe_spider_response(d, v3c_spec("spider_response", x = "week", y = "lesion", id = "patient")),
               paste("spider_response", msg), fixed = TRUE)
})

test_that("the spider and swimmer data logs count each drop reason on its own (F124)", {
  d <- data.frame(patient = c("P1", "P1", "P2", "P2", "P3", "P4", "P4", "P5"),
                  week = c(0, 8, 0, 8, 0, 0, 8, NA),
                  lesion = c(50, 40, 0, 10, 30, 20, 25, 33))
  res <- recipe_spider_response(d, v3c_spec("spider_response", x = "week", y = "lesion", id = "patient"))
  expect_equal(res$clean_steps, c(
    "Kept 2 patient(s) across 4 measurement(s).",
    "Dropped 1 row(s) with a missing patient, time or value.",
    "Dropped 1 patient(s) (2 row(s)) whose baseline was zero, because a percentage change needs a non-zero baseline.",
    "Dropped 1 patient(s) (1 row(s)) with no visit after the baseline, because a response needs a second measurement."))
  s <- data.frame(patient = c("P1", "P2", "P3", "P4"), months = c(4, -2, NA, 7))
  res <- recipe_swimmer(s, v3c_spec("swimmer", id = "patient", time = "months"))
  expect_equal(res$clean_steps, c("Dropped 1 row(s) with a missing patient, duration, status or arm.",
                                  "Dropped 1 row(s) with a negative duration.",
                                  "Kept 2 of 4 patient(s)."))
})
