#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })

raw <- read.csv("input_expression_matrix.csv", check.names = FALSE, row.names = 1)
mat <- as.matrix(raw); storage.mode(mat) <- "double"
keep <- apply(mat, 1, function(r) all(is.finite(r)) && sd(r) > 0); mat <- mat[keep, ]
disp <- t(scale(t(mat)))
disp[!is.finite(disp)] <- 0
ann <- read.csv("expression_annotation.csv", check.names = FALSE, row.names = 1)
ann <- ann[colnames(mat), , drop = FALSE]
top <- HeatmapAnnotation(df = ann)
col_fun <- colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))
ht <- Heatmap(disp, name = "z-score", col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE, top_annotation = top)

cairo_pdf("figure_expression_matrix_demo.pdf", width = 7, height = 7)
draw(ht, merge_legends = TRUE); dev.off()
cat("Reproduced figure.\n")
