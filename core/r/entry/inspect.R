#!/usr/bin/env Rscript
# Entry: figkit inspect. Report the shape of a tidy data file so the agent and
# the scientist can pick x and y and reason about the design.
suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--data",   type = "character"),
  make_option("--sheet",  type = "character", default = NULL),
  make_option("--format", type = "character", default = "text")
)
opt <- parse_args(OptionParser(option_list = opts))
if (is.null(opt$data)) stop("missing required --data")

source(file.path(Sys.getenv("PUBPLOT_CORE", "/opt/pubplot/core"), "r", "bootstrap.R"))

df   <- read_tidy(opt$data, opt$sheet)
info <- inspect_data(df)

if (identical(opt$format, "json")) {
  cat(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, na = "null"))
  cat("\n")
} else {
  print_inspect(info)
}
