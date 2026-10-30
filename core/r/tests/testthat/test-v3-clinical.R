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
