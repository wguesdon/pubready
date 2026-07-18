#!/usr/bin/env Rscript
# Entry: figkit render. Rebuild a figure from an edited plot_config.yaml.
suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--config", type = "character"),
  make_option("--data",   type = "character", default = NULL),
  make_option("--sheet",  type = "character", default = NULL),
  make_option("--out",    type = "character", default = "pubplot_output")
)
opt <- parse_args(OptionParser(option_list = opts))
if (is.null(opt$config)) stop("missing required --config")

source(file.path(Sys.getenv("PUBPLOT_CORE", "/opt/pubplot/core"), "r", "bootstrap.R"))

spec <- spec_from_config(opt$config)
raw_input <- opt$data %||% spec$data$file
if (is.null(raw_input)) stop("no input data: pass --data or set data.file in the config")

run_recipe(spec, raw_input = raw_input, out_root = opt$out, sheet = opt$sheet)
