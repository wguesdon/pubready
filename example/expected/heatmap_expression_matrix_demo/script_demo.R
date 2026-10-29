#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })

# The leaf order of the figure: complete linkage on the Euclidean distance of
# the drawn values, and at each merge the branch with the higher mean of the
# unscaled values first. On a tie, the branch with the earlier input index.
.hm_branch <- function (k, w, sum_w, size, first) 
{
    if (k < 0) 
        c(w[-k], 1, -k)
    else c(sum_w[k], size[k], first[k])
}

.hm_order_tree <- function (hc, w) 
{
    m <- hc$merge
    nm <- nrow(m)
    sum_w <- numeric(nm)
    size <- numeric(nm)
    first <- numeric(nm)
    for (i in seq_len(nm)) {
        a <- .hm_branch(m[i, 1], w, sum_w, size, first)
        b <- .hm_branch(m[i, 2], w, sum_w, size, first)
        ma <- a[1]/a[2]
        mb <- b[1]/b[2]
        if (mb > ma || (mb == ma && b[3] < a[3])) 
            m[i, ] <- m[i, 2:1]
        sum_w[i] <- a[1] + b[1]
        size[i] <- a[2] + b[2]
        first[i] <- min(a[3], b[3])
    }
    hc$merge <- m
    hc$order <- .hm_leaf_order(m)
    hc
}

.hm_leaf_order <- function (m) 
{
    out <- integer(0)
    stack <- nrow(m)
    while (length(stack) > 0) {
        k <- stack[1]
        stack <- stack[-1]
        if (k < 0) 
            out <- c(out, -k)
        else stack <- c(m[k, 1], m[k, 2], stack)
    }
    out
}

raw <- read.csv("input_expression_matrix.csv", check.names = FALSE, row.names = 1)
mat <- as.matrix(raw); storage.mode(mat) <- "double"
keep <- apply(mat, 1, function(r) all(is.finite(r)) && sd(r) > 0); mat <- mat[keep, ]
disp <- t(scale(t(mat)))
disp[!is.finite(disp)] <- 0
ann <- read.csv("expression_annotation.csv", check.names = FALSE, row.names = 1)
ann <- ann[colnames(mat), , drop = FALSE]
top <- HeatmapAnnotation(df = ann)
col_fun <- colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))
rows <- if (nrow(disp) > 1) as.dendrogram(.hm_order_tree(hclust(dist(disp), "complete"), rowMeans(mat))) else FALSE
cols <- if (ncol(disp) > 1) as.dendrogram(.hm_order_tree(hclust(dist(t(disp)), "complete"), colMeans(mat))) else FALSE
ht <- Heatmap(disp, name = "z-score", col = col_fun, cluster_rows = rows, cluster_columns = cols, top_annotation = top)

cairo_pdf("figure_expression_matrix_demo.pdf", width = 7, height = 7)
draw(ht, merge_legends = TRUE); dev.off()
cat("Reproduced figure.\n")
