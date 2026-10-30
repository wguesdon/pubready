# Tests of the area B fixes of the v3 findings: the read rule, the bundle, the
# spec and the inspect report. The Python suite (tests/test_v3_core.py) holds
# the same cases with the same expected values.

write_tmp_csv <- function(lines) {
  path <- tempfile(fileext = ".csv")
  writeLines(lines, path)
  path
}

test_that("a decimal comma stops the read and names the column (F40)", {
  path <- write_tmp_csv(c("group,value", 'control,"46,2"', "treated,51.0"))
  expect_error(read_tidy(path),
               "column 'value' holds '46,2', a number with a comma", fixed = TRUE)
})

test_that("a thousands mark stops the read and names the column (F40)", {
  path <- write_tmp_csv(c("group,value", 'control,"4,380"', "treated,4570"))
  expect_error(read_tidy(path),
               "column 'value' holds '4,380', a number with a comma", fixed = TRUE)
})

test_that("a text value with a comma stays text (F40)", {
  path <- write_tmp_csv(c("id,state", '1,"edema, no diuretics"', "2,no edema"))
  df <- read_tidy(path)
  expect_equal(df$state, c("edema, no diuretics", "no edema"))
})

test_that("only an empty cell or NA reads as missing (F41)", {
  path <- write_tmp_csv(c("arm,response,score", "a,None,1", "a,N/A,", "b,NULL,NA",
                          "b,nan,4", "c,,5"))
  rd <- read_tidy_logged(path)
  expect_equal(rd$df$response, c("None", "N/A", "NULL", "nan", NA))
  expect_equal(rd$df$score, c(1, NA, NA, 4, 5))
  expect_equal(rd$steps,
               "Read 3 cell(s) that hold \"NA\" or nothing as missing: 'response' (1), 'score' (2).")
})

test_that("the spaces around a value and a column name are removed (F41)", {
  path <- write_tmp_csv(c("group , value", "control ,52.2", " control, 49.9", "treated,48.6"))
  rd <- read_tidy_logged(path)
  expect_equal(names(rd$df), c("group", "value"))
  expect_equal(rd$df$group, c("control", "control", "treated"))
  expect_equal(rd$df$value, c(52.2, 49.9, 48.6))
  expect_equal(rd$steps, c(
    "Removed the spaces around 2 column name(s).",
    "Removed the spaces around 3 value(s): 'group' (2), 'value' (1)."))
})

test_that("a column with one text cell stays text, and logical tokens read as logical (F41)", {
  path <- write_tmp_csv(c("a,b,c,d", "1,TRUE,T,1", "nan,false,F,2", "3,,T,3"))
  df <- read_tidy(path)
  expect_equal(df$a, c("1", "nan", "3"))
  expect_equal(df$b, c(TRUE, FALSE, NA))
  expect_equal(df$c, c("T", "F", "T"))
  expect_equal(df$d, c(1, 2, 3))
})

test_that("a subnormal p value keeps its value (F128)", {
  path <- write_tmp_csv(c("snp,p", "rs1,1e-320", "rs2,0.5", "rs3,2.5E-10"))
  df <- read_tidy(path)
  expect_identical(df$p[1], 1e-320)
  expect_true(df$p[1] > 0 && df$p[1] < .Machine$double.xmin)
  expect_identical(df$p[3], 2.5e-10)
})
