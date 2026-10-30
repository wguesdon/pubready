# The v3 fixes of the comparison recipes. core/python/tests/test_v3_compare.py
# holds the same cases for the Python engine.

v3_spec <- function(recipe, ...) {
  spec_from_opt(c(list(recipe = recipe, data = "toy.csv"), list(...)))
}

# F43: a constant sample has one status, and a test that divides by a variance
# of zero stops with a message that names the samples.
test_that("a constant sample is not normal and has no Shapiro P", {
  s <- shapiro_safe(rep(5, 6))
  expect_true(is.na(s$p))
  expect_false(s$normal)
})

test_that("a forced t-test on two constant groups stops and names them", {
  df <- data.frame(group = rep(c("a", "b"), each = 4), value = rep(c(5, 6), each = 4))
  expect_error(recipe_two_group_compare(df, v3_spec("two_group_compare", x = "group", y = "value",
                                                    test = "welch_t")),
               "the groups 'a' and 'b' have zero variance")
})

test_that("Welch's ANOVA on a constant group stops and names it", {
  df <- data.frame(group = rep(c("het", "ko", "wt"), each = 4),
                   value = c(5, 5, 5, 5, 1, 2, 3, 4, 2, 4, 6, 9))
  expect_error(recipe_multi_group_compare(df, v3_spec("multi_group_compare", x = "group",
                                                      y = "value", test = "welch")),
               "the group 'het' has zero variance")
})

test_that("a paired t-test on constant differences stops", {
  df <- data.frame(subject = rep(1:5, 2), condition = rep(c("pre", "post"), each = 5),
                   value = c(1:5, 1:5 + 3))
  expect_error(recipe_paired_compare(df, v3_spec("paired_compare", x = "condition", y = "value",
                                                 id = "subject", test = "t")),
               "within-subject differences have zero variance")
})

test_that("correlation stops on a constant column and names it", {
  df <- data.frame(a = c(1, 2, 3, 4, 5), b = rep(2, 5))
  expect_error(recipe_correlation(df, v3_spec("correlation", x = "a", y = "b")),
               "the column 'b' is constant")
})

test_that("feature_panel reports a constant feature as not tested", {
  df <- data.frame(group = rep(c("a", "b"), each = 6),
                   f1 = c(1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5),
                   floor = rep(0.5, 12))
  res <- recipe_feature_panel(df, v3_spec("feature_panel", x = "group"))
  row <- res$stats[res$stats$feature == "floor", ]
  expect_equal(row$note, "not tested: constant")
  expect_true(is.na(row$p_value))
  expect_equal(res$test_meta$n_constant, 1)
  expect_equal(res$test_meta$assumptions$normality$n_checked, 2)
  expect_match(res$methods, "1 feature\\(s\\) were constant")
})
