# The cox_forest figure: the hazard-ratio axis ticks and the display rows. The
# Python engine has the same cases in test_cox_forest.py.
test_that(".cox_log_ticks covers the CIs on both sides of 1", {
  expect_equal(.cox_log_ticks(0.34, 1.94), c(0.5, 1, 2))
  expect_equal(.cox_log_ticks(0.44, 3.54), c(0.5, 1, 2))
  expect_equal(.cox_log_ticks(0.8, 1.3), c(0.7, 1, 1.5))
  expect_equal(.cox_log_ticks(0.99, 1.02), c(0.8, 1, 1.25))
  expect_equal(.cox_log_ticks(0.001, 500), c(0.001, 0.01, 0.1, 1, 10, 100))
})

test_that(".cox_forest_stars stops at three stars", {
  expect_equal(.cox_forest_stars(c(1e-9, 5e-4, 5e-3, 0.02, 0.05, 0.3, NA)),
               c("***", "***", "**", "*", "", "", ""))
})

test_that("recipe_cox_forest shows a reference row per factor", {
  set.seed(1)
  n <- 120
  df <- data.frame(time = stats::rexp(n, 0.1), event = stats::rbinom(n, 1, 0.7),
                   arm = sample(c("control", "treated"), n, TRUE),
                   age = stats::rnorm(n, 60, 8),
                   stage = sample(c("I", "II", "III"), n, TRUE),
                   stringsAsFactors = FALSE)
  spec <- spec_from_opt(list(recipe = "cox_forest", data = "trial.csv", time = "time",
                             event = "event", covariates = "arm,age,stage"))
  res <- recipe_cox_forest(df, spec)
  expect_equal(res$stats$term, c("armtreated", "age", "stageII", "stageIII"))
  cox <- survival::coxph(survival::Surv(time, event) ~ arm + age + stage,
                         data = transform(df, arm = factor(arm), stage = factor(stage)))
  rows <- .cox_forest_rows(cox, res$stats, c("arm", "age", "stage"))
  expect_equal(rows$label, c("arm: control", "arm: treated", "age",
                             "stage: I", "stage: II", "stage: III"))
  expect_equal(rows$hr_txt[rows$ref], c("reference", "reference"))
  expect_equal(rows$HR[!rows$ref], res$stats$HR)
})
