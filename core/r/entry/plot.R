#!/usr/bin/env Rscript
# Entry: figkit plot. Build a spec from CLI args and run the recipe.
suppressPackageStartupMessages(library(optparse))

source(file.path(Sys.getenv("PUBREADY_CORE", "/opt/pubready/core"), "r", "bootstrap.R"))

opt <- parse_args(OptionParser(option_list = plot_option_list()))

for (req in c("recipe", "data")) {
  if (is.null(opt[[req]])) stop(sprintf("missing required --%s", req))
}
check_recipe_options(opt$recipe, given_options(commandArgs(trailingOnly = TRUE)))

spec <- spec_from_opt(opt)
run_recipe(spec, raw_input = opt$data, out_root = opt$out, sheet = opt$sheet,
           stamp = opt$stamp)
