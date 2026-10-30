make_spec <- function(x = "group", y = "value") {
  list(data = list(x = x, y = y),
       test = list(method = "auto", paired = FALSE, p_adjust = "none"),
       appearance = default_appearance())
}

test_that("clean_tidy passes clean data through and logs that nothing changed", {
  df <- data.frame(group = c("A", "A", "B", "B"), value = c(1, 2, 3, 4))
  cl <- clean_tidy(df, make_spec())
  expect_equal(nrow(cl$df), 4)
  expect_match(cl$steps[1], "No cleaning applied")
})

test_that("clean_tidy drops missing rows and records the row counts", {
  df <- data.frame(group = c("A", "A", NA, "B"), value = c(1, NA, 3, 4))
  cl <- clean_tidy(df, make_spec())
  expect_equal(nrow(cl$df), 2)
  expect_true(any(grepl("Dropped 2 row", cl$steps)))
})

test_that("clean_tidy errors on a missing column", {
  df <- data.frame(group = c("A", "B"), value = c(1, 2))
  expect_error(clean_tidy(df, make_spec(y = "missing")), "not found")
})

test_that("coerce_numeric_columns keeps a mostly numeric column", {
  df <- data.frame(tnf = c("1.5", "2.0", "<0.1", "3.2", "4.1"),
                   group = c("a", "a", "b", "b", "b"),
                   il6 = c(1, 2, 3, 4, 5), stringsAsFactors = FALSE)
  out <- coerce_numeric_columns(df)
  expect_equal(out$numeric, c("tnf", "il6"))
  expect_equal(out$text, "group")
  expect_equal(sum(is.na(out$df$tnf)), 1)
  expect_equal(out$steps,
               "Coerced column 'tnf' to numeric; 1 value(s) became NA, for example '<0.1'.")
})

test_that("coerce_numeric_columns reads a logical column as text", {
  df <- data.frame(flag = c(TRUE, FALSE, TRUE), obj = c(TRUE, NA, FALSE), v = c(1, 2, 3))
  out <- coerce_numeric_columns(df)
  expect_equal(out$numeric, "v")
  expect_equal(out$text, c("flag", "obj"))
})
