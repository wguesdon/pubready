#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ComplexHeatmap))

raw <- read.csv("input_set_membership.csv", check.names = FALSE, row.names = 1)
bin <- as.data.frame(lapply(raw, function(v) as.numeric(v) == 1))
bin <- bin[rowSums(as.matrix(bin)) > 0, , drop = FALSE]
m <- make_comb_mat(bin)

cairo_pdf("figure_set_membership_upset.pdf", width = 7, height = 4.5)
draw(UpSet(m)); dev.off()
cat("Reproduced figure.\n")
