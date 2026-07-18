# Recipe: clustered heatmap of a features-by-samples matrix (ComplexHeatmap).
# Input is a matrix CSV: first column = feature id, remaining columns = samples.
# An optional annotation CSV (--annotation) labels the columns; its first column
# is the sample id. By default rows are z-scored and both axes are clustered.
# When the primary annotation column has exactly two groups, each feature gets a
# Welch t-test (BH-adjusted) and significance stars are drawn beside it.

recipe_heatmap <- function(df, spec) {
  feat <- as.character(df[[1]])
  mat  <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- feat

  n0 <- nrow(mat)
  keep <- apply(mat, 1, function(r) all(is.finite(r)) && stats::sd(r) > 0)
  mat <- mat[keep, , drop = FALSE]
  clean_steps <- if (sum(!keep) > 0) {
    sprintf("Dropped %d feature(s) with missing values or zero variance (%d -> %d).",
            sum(!keep), n0, nrow(mat))
  } else "No cleaning applied; input used as-is."
  df_used <- df[df[[1]] %in% rownames(mat), , drop = FALSE]

  scale_mode <- spec$appearance$scale %||% "row"
  cluster    <- spec$appearance$cluster %||% "both"
  disp <- switch(scale_mode,
    row    = t(scale(t(mat))),
    column = scale(mat),
    none   = mat,
    t(scale(t(mat))))
  disp[!is.finite(disp)] <- 0
  legend_name <- if (scale_mode == "none") "value" else "z-score"

  cluster_rows <- cluster %in% c("both", "rows")
  cluster_cols <- cluster %in% c("both", "columns")

  # Optional column annotation.
  top_anno <- NULL; ann <- NULL; group_col <- NULL
  if (!is.null(spec$data$annotation) && nzchar(spec$data$annotation)) {
    araw <- read_tidy(spec$data$annotation)
    rn <- as.character(araw[[1]])
    ann <- araw[, -1, drop = FALSE]
    rownames(ann) <- rn
    ann <- ann[colnames(mat), , drop = FALSE]
    for (cn in names(ann)) ann[[cn]] <- as.character(ann[[cn]])
    ann_cols <- lapply(names(ann), function(cn) {
      lv <- sort(unique(ann[[cn]]))
      stats::setNames(rep(pubplot_palette, length.out = length(lv)), lv)
    })
    names(ann_cols) <- names(ann)
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = ann, col = ann_cols,
      annotation_name_gp = grid::gpar(fontsize = 9))
    group_col <- names(ann)[1]
  }

  # Per-feature differential test when the primary annotation has two groups.
  right_anno <- NULL
  if (!is.null(group_col) && length(unique(ann[[group_col]])) == 2) {
    g  <- factor(ann[[group_col]]); lv <- levels(g)
    p  <- apply(mat, 1, function(r)
      tryCatch(stats::t.test(r[g == lv[1]], r[g == lv[2]])$p.value,
               error = function(e) NA_real_))
    padj  <- stats::p.adjust(p, "BH")
    stars <- ifelse(is.na(padj), "",
             ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**",
             ifelse(padj < 0.05, "*", ""))))
    right_anno <- ComplexHeatmap::rowAnnotation(
      ` ` = ComplexHeatmap::anno_text(stars, gp = grid::gpar(fontsize = 10)))
    stats_df <- data.frame(
      feature = rownames(mat),
      mean_g1 = rowMeans(mat[, g == lv[1], drop = FALSE]),
      mean_g2 = rowMeans(mat[, g == lv[2], drop = FALSE]),
      p_value = p, p_adj = padj, significant = padj < 0.05,
      stringsAsFactors = FALSE)
    names(stats_df)[2:3] <- paste0("mean_", lv)
    test_label <- sprintf("per-feature Welch t-test (%s vs %s), BH-adjusted", lv[1], lv[2])
    n_sig <- sum(padj < 0.05, na.rm = TRUE)
  } else {
    stats_df <- data.frame(feature = rownames(mat), mean = rowMeans(mat),
                           sd = apply(mat, 1, stats::sd), stringsAsFactors = FALSE)
    test_label <- NULL; n_sig <- NA_integer_
  }

  col_fun <- if (scale_mode == "none") {
    q <- stats::quantile(disp, c(0.02, 0.5, 0.98), na.rm = TRUE)
    circlize::colorRamp2(q, c("#3B6DB3", "white", "#C1432B"))
  } else {
    circlize::colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))
  }

  ht <- ComplexHeatmap::Heatmap(
    disp, name = legend_name, col = col_fun,
    cluster_rows = cluster_rows, cluster_columns = cluster_cols,
    top_annotation = top_anno, right_annotation = right_anno,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 9, fontface = "bold")))
  draw_fn <- function() ComplexHeatmap::draw(ht, merge_legends = TRUE)

  test_meta <- list(
    name = "Clustered heatmap (ComplexHeatmap)",
    scale = scale_mode, clustering = cluster,
    n_features = nrow(mat), n_samples = ncol(mat),
    per_feature_test = test_label, n_significant = n_sig)

  resolved <- list(family = "heatmap", label = "heatmap")
  methods  <- .hm_methods(scale_mode, cluster, nrow(mat), ncol(mat), test_label, n_sig)
  build_script <- function(in_name, fig_stub) {
    .hm_emit_script(spec, in_name, fig_stub, scale_mode, cluster_rows, cluster_cols, legend_name)
  }

  list(plot = NULL, draw = draw_fn, stats = stats_df, test_meta = test_meta,
       resolved = resolved, methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       extra_inputs = spec$data$annotation,
       width = max(5.5, 2.2 + 0.24 * ncol(mat)),
       height = max(5.0, 1.5 + 0.16 * nrow(mat)))
}

.hm_methods <- function(scale_mode, cluster, n_feat, n_samp, test_label, n_sig) {
  v    <- pkg_versions(c("ComplexHeatmap", "circlize"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  scale_txt <- switch(scale_mode,
    row = "Values were z-scored across samples within each feature (row scaling).",
    column = "Values were z-scored within each sample (column scaling).",
    "Raw values were shown without scaling.")
  clust_txt <- switch(cluster,
    both = "Rows and columns were hierarchically clustered.",
    rows = "Rows were hierarchically clustered; column order was preserved.",
    columns = "Columns were hierarchically clustered; row order was preserved.",
    "No clustering was applied.")
  test_txt <- if (!is.null(test_label)) {
    sprintf(" A %s was applied per feature; %d feature(s) were significant (adjusted P < 0.05) and are marked with stars.",
            test_label, n_sig)
  } else ""
  paste0(
    sprintf("A heatmap of %d features across %d samples was drawn with ComplexHeatmap. ",
            n_feat, n_samp),
    scale_txt, " ", clust_txt, test_txt,
    sprintf(" Rendered in R %s with ComplexHeatmap %s and circlize %s.",
            rver, v$ComplexHeatmap, v$circlize))
}

.hm_emit_script <- function(spec, in_name, fig_stub, scale_mode, cluster_rows, cluster_cols, legend_name) {
  ann_name <- if (!is.null(spec$data$annotation) && nzchar(spec$data$annotation)) {
    basename(spec$data$annotation)
  } else NULL
  disp_line <- switch(scale_mode,
    row    = "disp <- t(scale(t(mat)))",
    column = "disp <- scale(mat)",
    "disp <- mat")
  ann_lines <- if (!is.null(ann_name)) c(
    sprintf('ann <- read.csv("%s", check.names = FALSE, row.names = 1)', ann_name),
    "ann <- ann[colnames(mat), , drop = FALSE]",
    "top <- HeatmapAnnotation(df = ann)"
  ) else "top <- NULL"

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE, row.names = 1)', in_name),
    "mat <- as.matrix(raw); storage.mode(mat) <- \"double\"",
    "keep <- apply(mat, 1, function(r) all(is.finite(r)) && sd(r) > 0); mat <- mat[keep, ]",
    disp_line,
    "disp[!is.finite(disp)] <- 0",
    ann_lines,
    sprintf('col_fun <- colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))'),
    sprintf('ht <- Heatmap(disp, name = "%s", col = col_fun, cluster_rows = %s, cluster_columns = %s, top_annotation = top)',
            legend_name, if (cluster_rows) "TRUE" else "FALSE", if (cluster_cols) "TRUE" else "FALSE"),
    "",
    sprintf('cairo_pdf("%s.pdf", width = 7, height = 7)', fig_stub),
    "draw(ht, merge_legends = TRUE); dev.off()",
    'cat("Reproduced figure.\\n")'
  )
}
