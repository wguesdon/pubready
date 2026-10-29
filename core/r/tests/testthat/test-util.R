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

test_that("p_stars mirrors util.p_stars of the Python engine", {
  expect_equal(p_stars(c(1e-5, 5e-4, 5e-3, 2e-2, 0.2, NA)),
               c("****", "***", "**", "*", "ns", "ns"))
  expect_equal(p_stars(1e-4), "***")
  expect_equal(p_stars(0.05), "ns")
})

test_that("p_phrase writes P < 0.001 without an equals sign", {
  expect_equal(p_phrase(1e-5), "P < 0.001")
  expect_equal(p_phrase(0.0091), "P = 0.0091")
  expect_equal(p_phrase(0.22), "P = 0.22")
  expect_equal(p_phrase(NA_real_), "P = NA")
})

test_that("sig_text keeps the trailing zeros and drops a trailing point", {
  expect_equal(sig_text(3.7), "3.70")
  expect_equal(sig_text(3), "3.00")
  expect_equal(sig_text(213), "213")
  expect_equal(sig_text(c(4.58, -22.04)), c("4.58", "-22.0"))
  expect_equal(sig_text(0.0104, 2), "0.010")
  expect_equal(sig_text(7.1e-5, 2), "7.1e-05")
})

test_that("df_text writes a whole df as an integer and a Welch df with two decimals", {
  expect_equal(df_text(5), "5")
  expect_equal(df_text(29.9520363861042), "29.95")
})
