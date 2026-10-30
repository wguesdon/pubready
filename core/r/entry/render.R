#!/usr/bin/env Rscript
# Entry: figkit render. Rebuild a figure from an edited plot_config.yaml.
suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--config", type = "character"),
  make_option("--data",   type = "character", default = NULL),
  make_option("--sheet",  type = "character", default = NULL),
  make_option("--stamp",  type = "character", default = NULL),
  make_option("--engine", type = "character", default = NULL),  # read by the figkit dispatcher for routing; ignored here
  make_option("--out",    type = "character", default = "pubready_output")
)
opt <- parse_args(OptionParser(option_list = opts))
if (is.null(opt$config)) stop("missing required --config")

source(file.path(Sys.getenv("PUBREADY_CORE", "/opt/pubready/core"), "r", "bootstrap.R"))

spec <- spec_from_config(opt$config)
# The config of a bundle names its input and its annotation by the bare file
# name, so each is found in the folder of the config first.
spec$data$annotation <- resolve_config_path(spec$data$annotation, opt$config)
raw_input <- opt$data %||% resolve_config_path(spec$data$file, opt$config)
if (is.null(raw_input)) stop("no input data: pass --data or set data.file in the config")

run_recipe(spec, raw_input = raw_input, out_root = opt$out, sheet = opt$sheet,
           stamp = opt$stamp)
