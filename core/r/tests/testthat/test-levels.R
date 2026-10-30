test_that("level_labels writes whole numbers with no decimal point", {
  expect_equal(level_labels(c(80, 160, NA, 0.5, 2.5e-9, 1e5)),
               c("80", "160", NA, "0.5", "2.5e-09", "100000"))
})

test_that("level_order sorts numbers as numbers", {
  expect_equal(level_order(c("0", "160", "235", "310", "80")),
               c("0", "80", "160", "235", "310"))
  expect_equal(level_order(c(80, 0, NA, 160)), c("0", "80", "160"))
})

test_that("level_order sorts text by code point, not by the locale", {
  expect_equal(level_order(c("control", "Drug", "non-USA", "USA")),
               c("Drug", "USA", "control", "non-USA"))
})

test_that("level_order follows --levels and leaves out unused levels", {
  expect_equal(level_order(c("wt", "ko", "het"), c("wt", "het", "ko", "rescue")),
               c("wt", "het", "ko"))
})

test_that("level_order stops on a value that --levels does not list", {
  expect_error(level_order(c("wt", "ko"), "wt"), "does not list the group\\(s\\) 'ko'")
})
