#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })

raw <- read.csv("input_correlation_vars.csv", check.names = FALSE)
mat <- as.matrix(raw[, sapply(raw, is.numeric), drop = FALSE])
M <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
col_fun <- colorRamp2(c(-1, 0, 1), c("#3B6DB3", "white", "#C1432B"))
cell_fun <- function(j, i, x, y, w, h, fill)
  grid::grid.text(sprintf("%.2f", M[i, j]), x, y, gp = grid::gpar(fontsize = 7))
ht <- Heatmap(M, name = "correlation", col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE, cell_fun = cell_fun, rect_gp = grid::gpar(col = "white", lwd = 1))

cairo_pdf("figure_correlation_vars_corrhm.pdf", width = 6, height = 6)
draw(ht, merge_legends = TRUE); dev.off()
cat("Reproduced figure.\n")
