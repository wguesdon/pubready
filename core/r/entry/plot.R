#!/usr/bin/env Rscript
# Entry: figkit plot. Build a spec from CLI args and run the recipe.
suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--recipe",   type = "character"),
  make_option("--data",     type = "character"),
  make_option("--x",        type = "character"),
  make_option("--y",        type = "character"),
  make_option("--test",     type = "character", default = "auto"),
  make_option("--paired",   action = "store_true", default = FALSE),
  make_option("--p_adjust", type = "character", default = "none"),
  make_option("--geom",     type = "character", default = "box"),
  make_option("--xlab",     type = "character", default = NULL),
  make_option("--ylab",     type = "character", default = NULL),
  make_option("--title",    type = "character", default = NULL),
  make_option("--palette",  type = "character", default = NULL),
  make_option("--sheet",    type = "character", default = NULL),
  make_option("--out",      type = "character", default = "pubplot_output")
)
opt <- parse_args(OptionParser(option_list = opts))

for (req in c("recipe", "data", "x", "y")) {
  if (is.null(opt[[req]])) stop(sprintf("missing required --%s", req))
}

source(file.path(Sys.getenv("PUBPLOT_CORE", "/opt/pubplot/core"), "r", "bootstrap.R"))

spec <- spec_from_opt(opt)
run_recipe(spec, raw_input = opt$data, out_root = opt$out, sheet = opt$sheet)
