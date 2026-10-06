#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(metafor))

raw <- read.csv("input_hazard_ratio.csv", check.names = FALSE)
yi <- log(raw[["hazard_ratio"]])
sei <- (log(raw[["ci_high"]]) - log(raw[["ci_low"]])) / (2 * qnorm(0.975))
fe <- rma(yi = yi, sei = sei, method = "FE"); print(fe)
re <- rma(yi = yi, sei = sei, method = "DL"); print(re)
pdf("figure_hazard_ratio_meta.pdf", width = 8.2, height = 1.6 + 0.3 * (length(yi) + 3))
forest(re, slab = raw[["trial"]], atransf = exp)
invisible(dev.off())
cat("Reproduced figure.\n")
