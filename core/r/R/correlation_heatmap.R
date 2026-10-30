# Recipe: clustered correlation heatmap among the numeric columns of a table.
# Input is a tidy CSV whose numeric columns are the variables. The correlation
# matrix among those columns is drawn with ComplexHeatmap, rows and columns
# hierarchically clustered, with each cell showing the coefficient and BH-
# adjusted significance stars. Pearson by default; --test spearman|kendall.

# A column is a variable when it has 3 finite values or more and some variance.
.ch_usable <- function(v) sum(is.finite(v)) >= 3 && isTRUE(stats::sd(v, na.rm = TRUE) > 0)

# Stop when the values of a column differ but their variance is not a finite
# positive number at double precision: values near 1e300 overflow it to Inf, and
# values near 1e-300 underflow it to 0 (F83). distribution uses the same check.
# The Python engine has the same rule and message.
.check_finite_variance <- function(v, col, recipe) {
  v <- v[is.finite(v)]
  if (length(v) < 2 || length(unique(v)) < 2) return(invisible(TRUE))
  s2 <- stats::var(v)
  if (!is.finite(s2) || s2 <= 0) {
    stop(sprintf("%s: the variance of column '%s' is not finite at double precision (largest absolute value %s); rescale the column, for example to other units or to log10",
                 recipe, col, sprintf("%.2g", max(abs(v)))))
  }
  invisible(TRUE)
}

recipe_correlation_heatmap <- function(df, spec) {
  # coerce_numeric_columns() decides what is numeric: a column with one text
  # cell such as "<0.1" stays a variable, the cell becomes NA, and the data log
  # names the cell and each text column left out (F42).
  co <- coerce_numeric_columns(df)
  if (length(co$numeric) < 2) {
    stop("correlation_heatmap needs at least 2 numeric columns")
  }
  mat <- as.matrix(co$df[, co$numeric, drop = FALSE])
  for (cn in colnames(mat)) .check_finite_variance(mat[, cn], cn, "correlation_heatmap")

  keepc <- apply(mat, 2, .ch_usable)
  mat <- mat[, keepc, drop = FALSE]
  if (ncol(mat) < 2) stop("correlation_heatmap needs at least 2 numeric columns with variance")
  clean_steps <- c(
    co$steps,
    if (length(co$text) > 0) sprintf("Left out the text column(s) %s.", .pca_name_list(co$text)),
    if (sum(!keepc) > 0) sprintf("Dropped %d column(s) with fewer than 3 values or no variance: %s.",
                                 sum(!keepc), .pca_name_list(co$numeric[!keepc])))
  if (length(clean_steps) == 0) clean_steps <- "No cleaning applied; numeric columns used as-is."
  df_used <- as.data.frame(mat)

  method <- tolower(spec$test$method %||% "auto")
  if (method %in% c("auto", "")) method <- "pearson"
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop(sprintf("unknown correlation method '%s' (use pearson, spearman, or kendall)", method))
  }

  M    <- stats::cor(mat, method = method, use = "pairwise.complete.obs")
  vars <- colnames(M); k <- length(vars)

  # Pairwise p-values, BH-adjusted over the upper triangle then mirrored.
  P <- matrix(NA_real_, k, k, dimnames = list(vars, vars))
  for (i in seq_len(k)) for (j in seq_len(k)) {
    if (i == j) { P[i, j] <- 0; next }
    P[i, j] <- tryCatch(stats::cor.test(mat[, i], mat[, j], method = method)$p.value,
                        error = function(e) NA_real_)
  }
  ut      <- upper.tri(P)
  P_ut    <- P[ut]
  padj_ut <- stats::p.adjust(P_ut, "BH")
  padj <- matrix(0, k, k, dimnames = list(vars, vars))
  padj[ut] <- padj_ut
  padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]

  cluster <- spec$appearance$cluster %||% "both"
  cluster_rows <- cluster %in% c("both", "rows")
  cluster_cols <- cluster %in% c("both", "columns")
  # The leaf rule of heatmap (.hm_tree): complete linkage on the Euclidean
  # distance of the matrix rows, and at each merge the branch with the higher mean
  # correlation first; on a tie, the branch with the earlier input variable.
  row_tree <- if (cluster_rows && k > 1) .hm_tree(M, rowMeans(M)) else NULL
  col_tree <- if (cluster_cols && k > 1) .hm_tree(t(M), colMeans(M)) else NULL
  row_order <- if (is.null(row_tree)) seq_len(k) else row_tree$order
  col_order <- if (is.null(col_tree)) seq_len(k) else col_tree$order

  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3B6DB3", "white", "#C1432B"))
  star_of <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
  }
  cell_fun <- function(j, i, x, y, width, height, fill) {
    lab <- if (i == j) "1" else paste0(sprintf("%.2f", M[i, j]), star_of(padj[i, j]))
    grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 7))
  }

  ht <- ComplexHeatmap::Heatmap(
    M, name = "correlation", col = col_fun,
    cluster_rows = if (is.null(row_tree)) FALSE else stats::as.dendrogram(row_tree),
    cluster_columns = if (is.null(col_tree)) FALSE else stats::as.dendrogram(col_tree),
    rect_gp = grid::gpar(col = "white", lwd = 1),
    cell_fun = cell_fun,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    heatmap_legend_param = list(
      title = "correlation", at = c(-1, -0.5, 0, 0.5, 1),
      title_gp = grid::gpar(fontsize = 9, fontface = "bold")))
  draw_fn <- function() ComplexHeatmap::draw(ht, merge_legends = TRUE)

  ii <- row(M)[ut]; jj <- col(M)[ut]
  stats_df <- data.frame(
    var1 = vars[ii], var2 = vars[jj],
    r = M[ut], p_value = P_ut, p_adj = padj_ut,
    significant = padj_ut < 0.05, stringsAsFactors = FALSE)

  n_sig <- sum(padj_ut < 0.05, na.rm = TRUE)
  test_meta <- list(
    name = sprintf("Clustered %s correlation heatmap (ComplexHeatmap)", method),
    method = method, clustering = cluster,
    n_variables = k, n_pairs = length(P_ut), n_significant = n_sig,
    row_order = as.list(vars[row_order]), column_order = as.list(vars[col_order]))

  resolved <- list(family = "correlation_heatmap",
                   label = sprintf("%s correlation heatmap", method))
  methods  <- .ch_methods(method, cluster, k, n_sig)
  build_script <- function(in_name, fig_stub) {
    .ch_emit_script(in_name, fig_stub, method, cluster_rows, cluster_cols, vars)
  }

  list(plot = NULL, draw = draw_fn, stats = stats_df, test_meta = test_meta,
       resolved = resolved, methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       width = max(4.8, 1.4 + 0.55 * k), height = max(4.5, 1.2 + 0.55 * k))
}

.ch_methods <- function(method, cluster, k, n_sig) {
  v    <- pkg_versions(c("ComplexHeatmap", "circlize"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  clust_txt <- switch(cluster,
    both = "Variables were hierarchically clustered on both axes.",
    rows = "Rows were hierarchically clustered; column order was preserved.",
    columns = "Columns were hierarchically clustered; row order was preserved.",
    "No clustering was applied.")
  paste0(
    sprintf("Pairwise %s correlations among %d variables were computed and drawn as a heatmap with ComplexHeatmap. ",
            cap_first(method), k),
    clust_txt,
    " Each cell shows the coefficient with Benjamini-Hochberg adjusted significance stars",
    sprintf(" (%d of the %d pairs were significant at adjusted P < 0.05).",
            n_sig, choose(k, 2)),
    sprintf(" Rendered in R %s with ComplexHeatmap %s and circlize %s.",
            rver, v$ComplexHeatmap, v$circlize))
}

.ch_emit_script <- function(in_name, fig_stub, method, cluster_rows, cluster_cols, vars) {
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    "# The variables of the figure; a cell that does not read as a number is NA.",
    sprintf("vars <- %s", paste(deparse(vars, width.cutoff = 500L), collapse = "")),
    "num <- raw[vars]",
    "for (v in vars) num[[v]] <- suppressWarnings(as.numeric(as.character(num[[v]])))",
    "mat <- as.matrix(num)",
    sprintf('M <- cor(mat, method = "%s", use = "pairwise.complete.obs")', method),
    "# The leaf order of the figure: complete linkage on the Euclidean distance of",
    "# the matrix rows, and at each merge the branch with the higher mean correlation",
    "# first. On a tie, the branch with the earlier input variable.",
    .hm_fn_lines(".hm_branch", .hm_branch),
    .hm_fn_lines(".hm_order_tree", .hm_order_tree),
    .hm_fn_lines(".hm_leaf_order", .hm_leaf_order),
    if (cluster_rows) 'rows <- as.dendrogram(.hm_order_tree(hclust(dist(M), "complete"), rowMeans(M)))' else "rows <- FALSE",
    if (cluster_cols) 'cols <- as.dendrogram(.hm_order_tree(hclust(dist(t(M)), "complete"), colMeans(M)))' else "cols <- FALSE",
    'col_fun <- colorRamp2(c(-1, 0, 1), c("#3B6DB3", "white", "#C1432B"))',
    "cell_fun <- function(j, i, x, y, w, h, fill)",
    '  grid::grid.text(sprintf("%.2f", M[i, j]), x, y, gp = grid::gpar(fontsize = 7))',
    'ht <- Heatmap(M, name = "correlation", col = col_fun, cluster_rows = rows, cluster_columns = cols, cell_fun = cell_fun, rect_gp = grid::gpar(col = "white", lwd = 1))',
    "",
    sprintf('cairo_pdf("%s.pdf", width = 6, height = 6)', fig_stub),
    "draw(ht, merge_legends = TRUE); dev.off()",
    'cat("Reproduced figure.\\n")'
  )
}
