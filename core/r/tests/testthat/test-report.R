test_that("report_markdown_table escapes table values", {
  frame <- data.frame(name = "a|b", value = 2)

  table <- report_markdown_table(frame)

  expect_match(table, "a\\\\|b", fixed = TRUE)
  expect_match(table, "| name | value |", fixed = TRUE)
})
