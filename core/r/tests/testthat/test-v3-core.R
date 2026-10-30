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
  "ref", "cause", "levels", "visit_levels", "features", "test", "paired", "scale",
  "cluster", "fc_cutoff", "p_cutoff", "top_n",
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

test_that("the report table writes a factor with NA and no warning (F120)", {
  df <- data.frame(group = factor(c("a", NA, "b|c")), n = c(1, NA, 3))
  expect_silent(txt <- report_markdown_table(df))
  expect_equal(strsplit(txt, "\n")[[1]],
               c("| group | n |", "| --- | --- |", "| a | 1 |", "|  |  |", "| b\\|c | 3 |"))
})

test_that("the report table writes whole numbers as Python does (F126)", {
  df <- data.frame(n = c(1, NA, 3), p = c(0.5, Inf, NaN))
  expect_equal(strsplit(report_markdown_table(df), "\n")[[1]],
               c("| n | p |", "| --- | --- |", "| 1 | 0.5 |", "|  | Inf |", "| 3 |  |"))
})

# The inspect report, the same words in both engines (F72).
inspect_csv <- c("extra,group,ID,sex", "0.7,1,1,f", "-1.6,1,2,m", "-0.2,1,3,f", "-1.2,1,4,m",
                 "-0.1,1,5,f", "1.9,2,1,f", "0.8,2,2,m", "1.1,2,3,f", "0.1,2,4,m", "-0.1,2,5,f")
inspect_expected <- c(
  "Rows: 10   Columns: 4",
  "",
  "  extra                numeric   missing=0  unique=9  [-1.6 .. 1.9, mean 0.14]",
  paste0("  group                numeric   missing=0  unique=2  levels: \"1\" (5), \"2\" (5)",
         "  [1 .. 2, mean 1.5]"),
  paste0("  ID                   numeric   missing=0  unique=5  levels: \"1\" (2), \"2\" (2), ",
         "\"3\" (2), \"4\" (2), \"5\" (2)  [1 .. 5, mean 3]"),
  "  sex                  text      missing=0  unique=2  levels: \"f\" (6), \"m\" (4)",
  "",
  "Candidate grouping columns (x):  group (numeric code), ID (numeric code), sex",
  "Candidate outcome columns (y):   extra",
  "Possible subject id columns:     ID (repeats across the levels of group)")

test_that("inspect gives level counts, numeric codes and a subject id (F72)", {
  info <- inspect_data(read_tidy(write_tmp_csv(inspect_csv)))
  expect_equal(capture.output(print_inspect(info)), inspect_expected)
  expect_equal(info$id_candidates[[1]],
               list(column = "ID", across = "group", n_subjects = 5L, complete = TRUE))
  expect_equal(unlist(info$group_candidates), c("group", "ID", "sex"))
  expect_equal(unlist(info$outcome_candidates), "extra")
})

# The options of figkit plot and figkit diagnose, the same set in both engines (F73).
plot_option_names <- c(
  "annotation", "ci_high", "ci_low", "cause", "cluster", "covariates", "data", "engine",
  "event", "facet", "fc_cutoff", "features", "fill", "geom", "group", "id", "label",
  "levels", "markers", "measure", "out", "p_cutoff", "paired", "palette",
  "pd_threshold", "positive", "pr_threshold", "recipe", "ref", "scale", "se", "sheet",
  "stamp", "test", "theme", "time", "title", "top_n", "visit_levels", "x", "xlab", "y",
  "ylab")

test_that("diagnose and plot accept the same full option list (F73)", {
  opts <- plot_option_list()
  flags <- vapply(opts, methods::slot, character(1), "long_flag")
  expect_equal(sort(sub("^--", "", flags), method = "radix"), sort(plot_option_names, method = "radix"))
  opt <- optparse::parse_args(optparse::OptionParser(option_list = opts),
                    args = c("--recipe", "survival_km", "--data", "d.csv", "--time", "t",
                             "--ref", "vehicle"))
  expect_equal(opt$time, "t")
  expect_equal(opt$ref, "vehicle")
  expect_null(opt$test)
  expect_equal(spec_from_opt(opt)$test$method, "auto")
})

test_that("render finds the bundle input in the folder of the config (F74)", {
  dir <- file.path(tempdir(), "bundle_f74")
  dir.create(dir, showWarnings = FALSE)
  writeLines("a,b", file.path(dir, "input_x.csv"))
  cfg <- file.path(dir, "plot_config.yaml")
  expect_equal(resolve_config_path("input_x.csv", cfg), file.path(dir, "input_x.csv"))
  expect_equal(resolve_config_path("example/other.csv", cfg), "example/other.csv")
  expect_equal(resolve_config_path("/abs/path.csv", cfg), "/abs/path.csv")
  expect_null(resolve_config_path(NULL, cfg))
})

test_that("--p_adjust is gone and an old config drops the key (F75)", {
  path <- tempfile(fileext = ".yaml")
  writeLines(c("engine: r", "recipe: two_group_compare", "data:", "  x: group",
               "test:", "  method: auto", "  p_adjust: bonferroni"), path)
  spec <- spec_from_config(path)
  expect_null(spec$test$p_adjust)
  expect_equal(spec$test$method, "auto")
  expect_false("p_adjust" %in% names(spec_from_opt(list(recipe = "x", data = "d.csv"))$test))
})

# The option table of the recipes, one file for both engines (F71).
test_that("the option table names each recipe and only parser options (F71)", {
  table <- read_recipe_options()
  recipes <- setdiff(names(table), "always")
  expect_true(all(vapply(paste0("recipe_", recipes), exists, logical(1), mode = "function")))
  ns_recipes <- sub("^recipe_", "", ls(environment(recipe_heatmap), pattern = "^recipe_"))
  expect_setequal(recipes, ns_recipes)
  flags <- sub("^--", "", vapply(plot_option_list(), methods::slot, character(1), "long_flag"))
  expect_length(setdiff(unlist(table), flags), 0)
})

test_that("an option that a recipe does not read stops the run (F71)", {
  table <- list(always = c("recipe", "data"), a = c("x", "title"), b = "x", c = "title")
  expect_null(check_recipe_options("a", c("recipe", "x", "title"), table))
  expect_null(check_recipe_options("zzz", c("title"), table))
  expect_error(check_recipe_options("b", c("recipe", "x", "title"), table),
               "--title is not an option of b; it is read by: a, c", fixed = TRUE)
  expect_error(check_recipe_options("b", "p_adjust", table),
               "--p_adjust is not an option of b; no recipe reads it", fixed = TRUE)
  expect_equal(given_options(c("--recipe", "heatmap", "--ylab=-log10 p", "-x", "--title", "t")),
               c("recipe", "ylab", "title"))
  expect_error(check_recipe_options("heatmap", "title"),
               "--title is not an option of heatmap; it is read by: two_group_compare, ",
               fixed = TRUE)
  expect_null(check_recipe_options("sankey", c("x", "y", "id", "levels", "visit_levels")))
})
