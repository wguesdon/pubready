# The paired recipe: paired t-test end to end with the conditions kept in order.
test_that("recipe_paired_compare runs the paired t-test", {
  pre  <- c(10, 11, 12, 13, 14, 15, 16, 17)
  post <- pre - 5 + 0.3 * (-1)^(seq_along(pre))
  df <- data.frame(
    subject   = rep(paste0("s", 1:8), 2),
    condition = rep(c("pre", "post"), each = 8),
    value     = c(pre, post), stringsAsFactors = FALSE)
  opt <- list(recipe = "paired_compare", data = "p.csv", x = "condition", y = "value",
              id = "subject", test = "t")
  spec <- spec_from_opt(opt)

  res <- recipe_paired_compare(df, spec)

  expect_equal(res$stats$n_pairs, 8)
  expect_equal(res$stats$test, "paired t-test")
  expect_equal(res$stats$condition1, "pre")
  expect_true(res$stats$p_value < 0.001)
})
