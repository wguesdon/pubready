# Recipe: clustered correlation heatmap among the numeric columns of a table.
# Input is a tidy CSV whose numeric columns are the variables. The correlation
# matrix among those columns is drawn with ComplexHeatmap, rows and columns
# hierarchically clustered, with each cell showing the coefficient and BH-
# adjusted significance stars. Pearson by default; --test spearman|kendall.

recipe_correlation_heatmap <- function(df, spec) {
  num <- vapply(df, is.numeric, logical(1))
  if (sum(num) < 2) {
    stop("correlation_heatmap needs at least 2 numeric columns")
  }
  mat <- as.matrix(df[, num, drop = FALSE])
  n0 <- ncol(mat)

  keepc <- apply(mat, 2, function(v) sum(is.finite(v)) >= 3 && stats::sd(v, na.rm = TRUE) > 0)
  mat <- mat[, keepc, drop = FALSE]
  if (ncol(mat) < 2) stop("correlation_heatmap needs at least 2 numeric columns with variance")
  clean_steps <- if (sum(!keepc) > 0) {
    sprintf("Dropped %d non-numeric or zero-variance column(s) (%d -> %d variables).",
            sum(!keepc), n0, ncol(mat))
  } else "No cleaning applied; numeric columns used as-is."
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
    cluster_rows = cluster_rows, cluster_columns = cluster_cols,
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
    n_variables = k, n_pairs = length(P_ut), n_significant = n_sig)

  resolved <- list(family = "correlation_heatmap",
                   label = sprintf("%s correlation heatmap", method))
  methods  <- .ch_methods(method, cluster, k, n_sig)
  build_script <- function(in_name, fig_stub) {
    .ch_emit_script(in_name, fig_stub, method, cluster_rows, cluster_cols)
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
            method, k),
    clust_txt,
    " Each cell shows the coefficient with Benjamini-Hochberg adjusted significance stars",
    sprintf(" (%d of the %d pairs were significant at adjusted P < 0.05).",
            n_sig, choose(k, 2)),
    sprintf(" Rendered in R %s with ComplexHeatmap %s and circlize %s.",
            rver, v$ComplexHeatmap, v$circlize))
}

.ch_emit_script <- function(in_name, fig_stub, method, cluster_rows, cluster_cols) {
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    "mat <- as.matrix(raw[, sapply(raw, is.numeric), drop = FALSE])",
    sprintf('M <- cor(mat, method = "%s", use = "pairwise.complete.obs")', method),
    'col_fun <- colorRamp2(c(-1, 0, 1), c("#3B6DB3", "white", "#C1432B"))',
    "cell_fun <- function(j, i, x, y, w, h, fill)",
    '  grid::grid.text(sprintf("%.2f", M[i, j]), x, y, gp = grid::gpar(fontsize = 7))',
    sprintf('ht <- Heatmap(M, name = "correlation", col = col_fun, cluster_rows = %s, cluster_columns = %s, cell_fun = cell_fun, rect_gp = grid::gpar(col = "white", lwd = 1))',
            if (cluster_rows) "TRUE" else "FALSE", if (cluster_cols) "TRUE" else "FALSE"),
    "",
    sprintf('cairo_pdf("%s.pdf", width = 6, height = 6)', fig_stub),
    "draw(ht, merge_legends = TRUE); dev.off()",
    'cat("Reproduced figure.\\n")'
  )
}
