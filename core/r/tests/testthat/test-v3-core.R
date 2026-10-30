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

test_that("the manifest keeps the full precision of each number (F68)", {
  txt <- manifest_json(list(intercept = 0.00508612345678912, p = 1.5e-12))
  expect_match(txt, "0.00508612345678912", fixed = TRUE)
  back <- jsonlite::fromJSON(txt)
  expect_identical(back$intercept, 0.00508612345678912)
  expect_identical(back$p, 1.5e-12)
})

test_that("the manifest writes NaN, infinity and NA as null (F69)", {
  txt <- manifest_json(list(a = NaN, b = Inf, c = NA_real_, d = c(1.5, -Inf)))
  expect_equal(gsub("\\s", "", as.character(txt)), '{"a":null,"b":null,"c":null,"d":[1.5,null]}')
})

# The arguments block of the manifest, the same list in both engines (F101).
manifest_arg_names <- c(
  "x", "y", "fill", "facet", "time", "event", "covariates", "annotation", "id",
  "group", "label", "positive", "se", "ci_low", "ci_high", "measure", "markers",
  "ref", "cause", "levels", "visit_levels", "features", "test", "paired",
  "p_adjust", "scale", "cluster", "fc_cutoff", "p_cutoff", "top_n",
  "pd_threshold", "pr_threshold", "geom", "theme", "xlab", "ylab", "title",
  "palette")

test_that("the manifest records every plot option (F101)", {
  spec <- spec_from_opt(list(recipe = "sankey", data = "d.csv", x = "week", y = "response",
                             id = "patient", levels = "CR|PR", covariates = "age",
                             xlab = "Week", test = "auto"))
  args <- manifest_arguments(spec)
  expect_equal(names(args), manifest_arg_names)
  txt <- gsub("\\s", "", as.character(manifest_json(list(arguments = args))))
  expect_match(txt, '"id":"patient"', fixed = TRUE)
  expect_match(txt, '"levels":["CR","PR"]', fixed = TRUE)
  expect_match(txt, '"covariates":["age"]', fixed = TRUE)
  expect_match(txt, '"fill":null', fixed = TRUE)
  expect_match(txt, '"xlab":"Week"', fixed = TRUE)
  expect_match(txt, '"paired":false', fixed = TRUE)
  expect_match(txt, '"top_n":15', fixed = TRUE)
})

# One CSV form for a stats table in both engines (F126).
stats_csv_expected <- c(
  "name,n,p,ok,d,flag",
  "a,3,0.5,TRUE,-Inf,TRUE",
  "b,NA,NA,FALSE,2,FALSE",
  "c,5,Inf,NA,NA,TRUE")

test_that("the stats CSV writes NA, booleans, infinity and integers in one form (F126)", {
  df <- data.frame(name = c("a", "b", "c"), n = c(3, NA, 5), p = c(0.5, NaN, Inf),
                   ok = c(TRUE, FALSE, NA), d = c(-Inf, 2, NA), flag = c(TRUE, FALSE, TRUE),
                   stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".csv")
  write_stats_csv(df, path)
  expect_equal(readLines(path), stats_csv_expected)
})

# The wrap rule of the QC header, the same lines in both engines (F119).
test_that("a one-panel QC figure wraps its title and subtitle (F119)", {
  lay <- qc_layout("differences between the methods (new_assay - reference_assay)", 1)
  expect_equal(lay$subtitle, c("Points near the line support normality;",
                               "mild, judged deviation is usually",
                               "acceptable."))
  expect_equal(lay$title, c("Normality QC:", "differences between the", "methods (new_assay -",
                            "reference_assay)"))
  expect_true(all(nchar(lay$title) <= 23))
  expect_equal(lay$width, 3.3)
  expect_equal(lay$height, 3.9 + 0.25 * 3 + 0.18 * 2)
  three <- qc_layout("group values", 3)
  expect_equal(length(three$subtitle), 1)
  expect_equal(three$title, "Normality QC: group values")
  res <- qc_normality(list(quantity = "x", panels = list(a = rnorm(10))))
  expect_equal(c(res$width, res$height), c(3.3, 3.9 + 0.18 * 2))
})
