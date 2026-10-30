test_that("report_markdown_table escapes table values", {
  frame <- data.frame(name = "a|b", value = 2)

  table <- report_markdown_table(frame)

  # One backslash escapes the pipe; two would render a backslash and split the cell.
  expect_match(table, "| a\\|b | 2 |", fixed = TRUE)
  expect_match(table, "| name | value |", fixed = TRUE)
})
