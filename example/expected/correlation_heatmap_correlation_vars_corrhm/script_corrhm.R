#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })

raw <- read.csv("input_correlation_vars.csv", check.names = FALSE)
# The variables of the figure; a cell that does not read as a number is NA.
vars <- c("il6", "tnf", "crp", "il1b", "il10", "tgfb", "il4", "ifng")
num <- raw[vars]
for (v in vars) num[[v]] <- suppressWarnings(as.numeric(as.character(num[[v]])))
mat <- as.matrix(num)
M <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
# The leaf order of the figure: complete linkage on the Euclidean distance of
# the matrix rows, and at each merge the branch with the higher mean correlation
# first. On a tie, the branch with the earlier input variable.
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

rows <- as.dendrogram(.hm_order_tree(hclust(dist(M), "complete"), rowMeans(M)))
cols <- as.dendrogram(.hm_order_tree(hclust(dist(t(M)), "complete"), colMeans(M)))
col_fun <- colorRamp2(c(-1, 0, 1), c("#3B6DB3", "white", "#C1432B"))
cell_fun <- function(j, i, x, y, w, h, fill)
  grid::grid.text(sprintf("%.2f", M[i, j]), x, y, gp = grid::gpar(fontsize = 7))
ht <- Heatmap(M, name = "correlation", col = col_fun, cluster_rows = rows, cluster_columns = cols, cell_fun = cell_fun, rect_gp = grid::gpar(col = "white", lwd = 1))

cairo_pdf("figure_correlation_vars_corrhm.pdf", width = 6, height = 6)
draw(ht, merge_legends = TRUE); dev.off()
cat("Reproduced figure.\n")
