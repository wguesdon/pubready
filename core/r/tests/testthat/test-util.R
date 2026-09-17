test_that("%||% returns the fallback only for NULL or empty", {
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% 2, 2)
  expect_equal(character(0) %||% "x", "x")
  expect_equal(0 %||% 9, 0)
  expect_equal(FALSE %||% 9, FALSE)
})

test_that("cap_first capitalizes the first letter and passes edge cases", {
  expect_equal(cap_first("hello"), "Hello")
  expect_equal(cap_first("Hello"), "Hello")
  expect_equal(cap_first(""), "")
})

test_that("shapiro_safe guards small and large samples", {
  small <- shapiro_safe(c(1, 2))
  expect_true(is.na(small$normal))

  big <- shapiro_safe(seq_len(6000))
  expect_equal(big$n, 6000L)
  expect_true(isTRUE(big$normal))

  norm <- shapiro_safe(qnorm(ppoints(40)))
  expect_true(isTRUE(norm$normal))
})

test_that("resolve_two_group_method picks the test from assumptions", {
  student <- resolve_two_group_method("auto", all_normal = TRUE, equal_var = TRUE, paired = FALSE)
  expect_equal(student$family, "t")
  expect_equal(student$label, "Student's two-sample t-test")

  welch <- resolve_two_group_method("auto", all_normal = TRUE, equal_var = FALSE, paired = FALSE)
  expect_equal(welch$label, "Welch two-sample t-test")

  mwu <- resolve_two_group_method("auto", all_normal = FALSE, equal_var = TRUE, paired = FALSE)
  expect_equal(mwu$family, "wilcoxon")
  expect_equal(mwu$label, "Mann-Whitney U test")
})

test_that("resolve_two_group_method honors forced tests and paired labels", {
  forced <- resolve_two_group_method("welch_t", all_normal = FALSE, equal_var = TRUE, paired = FALSE)
  expect_equal(forced$family, "t")
  expect_false(isTRUE(forced$var_equal))

  paired_t <- resolve_two_group_method("t", all_normal = TRUE, equal_var = TRUE, paired = TRUE)
  expect_equal(paired_t$label, "paired t-test")

  signed_rank <- resolve_two_group_method("wilcoxon", all_normal = FALSE, equal_var = TRUE, paired = TRUE)
  expect_equal(signed_rank$label, "Wilcoxon signed-rank test")

  expect_error(resolve_two_group_method("nonsense", TRUE, TRUE, FALSE), "unknown test method")
})

test_that("rank_biserial_unpaired maps U onto the full range", {
  expect_equal(rank_biserial_unpaired(0, 6, 6), -1)
  expect_equal(rank_biserial_unpaired(36, 6, 6), 1)
  expect_equal(rank_biserial_unpaired(18, 6, 6), 0)
})

test_that("rank_biserial_paired follows Kerby's formula and drops zeros", {
  expect_equal(rank_biserial_paired(c(1, 2, 3, -4)), 0.2)
  expect_equal(rank_biserial_paired(c(0, 1, 2, 3, -4)), 0.2)
  expect_true(is.na(rank_biserial_paired(c(0, 0))))
})
