# survival_km: the number at risk counts time >= t at breaks that both engines
# compute in the same way. The Python suite (test_survival_km.py) expects the
# same numbers.

km_data <- function() {
  data.frame(
    arm = rep(c("A", "B"), each = 9),
    months = c(2, 5, 5, 8, 10, 10, 14, 18, 21.8, 3, 6, 10, 12, 15, 17, 19, 20, 20),
    event = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1))
}

test_that("the number at risk counts the subjects whose time is at or after each break", {
  spec <- spec_from_opt(list(recipe = "survival_km", data = "d.csv", test = "auto",
                             time = "months", event = "event", x = "arm"))
  res <- recipe_survival_km(km_data(), spec)
  ar <- res$test_meta$at_risk
  expect_equal(ar$times, c(0, 5, 10, 15, 20))
  expect_equal(ar$counts$A, c(9L, 8L, 5L, 2L, 1L))
  expect_equal(ar$counts$B, c(9L, 8L, 7L, 5L, 2L))
  fit <- survival::survfit(survival::Surv(months, event) ~ arm, data = km_data())
  ref <- summary(fit, times = ar$times, extend = TRUE)$n.risk
  expect_equal(c(ar$counts$A, ar$counts$B), ref)
})

test_that("the x breaks follow the 1-2-5 step that gives at most six intervals", {
  expect_equal(km_break_step(21.8), 5)
  expect_equal(km_break_step(60), 10)
  expect_equal(km_break_step(1), 0.2)
  expect_equal(km_breaks(21.8, 5), c(0, 5, 10, 15, 20))
})
