# Recipe: clustered heatmap of a features-by-samples matrix (ComplexHeatmap).
# Input is a matrix CSV: first column = feature id, remaining columns = samples.
# An optional annotation CSV (--annotation) labels the columns; its first column
# is the sample id. By default rows are z-scored and both axes are clustered.
# With --test welch_t and an annotation whose first column holds two groups,
# each feature gets a Welch t-test (BH-adjusted) and significance stars are
# drawn beside it. The default, --test auto, runs no test (F97).
#
# Both engines draw the same leaf order (.hm_tree): a complete-linkage tree on
# the Euclidean distance of the drawn values, and at each merge the branch with
# the higher mean of the unscaled values first. On a tie, the branch that holds
# the earlier input row or column comes first. ComplexHeatmap would reorder by
# the means of the drawn values, which are rounding noise after a row z-score.

# The tree of one heatmap axis. `disp` has one row per leaf; `w` is the weight of
# each leaf, the mean of its unscaled values.
.hm_tree <- function(disp, w) {
  .hm_order_tree(stats::hclust(stats::dist(disp), method = "complete"), w)
}

# The summed weight, the size and the first input index of a merge entry: a
# negative entry is a leaf, a positive one an earlier merge.
.hm_branch <- function(k, w, sum_w, size, first) {
  if (k < 0) c(w[-k], 1, -k) else c(sum_w[k], size[k], first[k])
}

# Swap the two branches of each merge so that the branch with the higher mean
# weight comes first, or on a tie the branch with the earlier input index.
# as.dendrogram() draws the first entry of a merge row on the left.
.hm_order_tree <- function(hc, w) {
  m <- hc$merge
  nm <- nrow(m)
  sum_w <- numeric(nm); size <- numeric(nm); first <- numeric(nm)
  for (i in seq_len(nm)) {
    a <- .hm_branch(m[i, 1], w, sum_w, size, first)
    b <- .hm_branch(m[i, 2], w, sum_w, size, first)
    ma <- a[1] / a[2]; mb <- b[1] / b[2]
    if (mb > ma || (mb == ma && b[3] < a[3])) m[i, ] <- m[i, 2:1]
    sum_w[i] <- a[1] + b[1]; size[i] <- a[2] + b[2]; first[i] <- min(a[3], b[3])
  }
  hc$merge <- m
  hc$order <- .hm_leaf_order(m)
  hc
}

# The leaves of a merge matrix from left to right.
.hm_leaf_order <- function(m) {
  out <- integer(0)
  stack <- nrow(m)
  while (length(stack) > 0) {
    k <- stack[1]; stack <- stack[-1]
    if (k < 0) out <- c(out, -k) else stack <- c(m[k, 1], m[k, 2], stack)
  }
  out
}

# Whether to run the per-feature Welch test. --test auto, the default, draws no
# test; --test welch_t runs it and needs an annotation whose first column holds
# two groups. Any other value stops the run.
.hm_welch_requested <- function(method, group_col, ann) {
  method <- method %||% "auto"
  if (!method %in% c("auto", "welch_t")) {
    stop(sprintf("heatmap: unknown --test '%s'; use auto (no test) or welch_t", method))
  }
  if (method == "auto") return(FALSE)
  if (is.null(group_col)) {
    stop("heatmap: --test welch_t needs --annotation, a table whose first column after the sample id holds two groups")
  }
  lv <- unique(ann[[group_col]][!is.na(ann[[group_col]])])
  if (length(lv) != 2) {
    stop(sprintf("heatmap: --test welch_t compares two groups, but the annotation column '%s' holds %d",
                 group_col, length(lv)))
  }
  TRUE
}

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
  row_tree <- if (cluster_rows && nrow(disp) > 1) .hm_tree(disp, rowMeans(mat)) else NULL
  col_tree <- if (cluster_cols && ncol(disp) > 1) .hm_tree(t(disp), colMeans(mat)) else NULL
  row_order <- if (is.null(row_tree)) seq_len(nrow(disp)) else row_tree$order
  col_order <- if (is.null(col_tree)) seq_len(ncol(disp)) else col_tree$order

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
      stats::setNames(rep(pubready_palette, length.out = length(lv)), lv)
    })
    names(ann_cols) <- names(ann)
    top_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = ann, col = ann_cols,
      annotation_name_gp = grid::gpar(fontsize = 9))
    group_col <- names(ann)[1]
  }

  # Per-feature differential test, only with --test welch_t. On genes that a DE
  # test chose, a second test on the same samples is circular, so the default
  # draws no test (F97).
  run_test <- .hm_welch_requested(spec$test$method, group_col, ann)
  right_anno <- NULL
  if (run_test) {
    g  <- factor(ann[[group_col]]); lv <- levels(g)
    p  <- apply(mat, 1, function(r)
      tryCatch(stats::t.test(r[g == lv[1]], r[g == lv[2]])$p.value,
               error = function(e) NA_real_))
    padj  <- stats::p.adjust(p, "BH")
    stars <- ifelse(is.na(padj), "",
             ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**",
             ifelse(padj < 0.05, "*", ""))))
    # anno_text cannot size a column of empty strings, so a null result, with
    # no star to draw, gets no star column (F46).
    if (any(nzchar(stars))) {
      right_anno <- ComplexHeatmap::rowAnnotation(
        ` ` = ComplexHeatmap::anno_text(stars, gp = grid::gpar(fontsize = 10)))
    }
    stats_df <- data.frame(
      feature = rownames(mat),
      mean_g1 = rowMeans(mat[, g == lv[1], drop = FALSE]),
      mean_g2 = rowMeans(mat[, g == lv[2], drop = FALSE]),
      p_value = p, p_adj = padj, significant = padj < 0.05,
      stringsAsFactors = FALSE)
    names(stats_df)[2:3] <- paste0("mean_", lv)
    test_label <- sprintf("per-feature Welch t-test (%s vs %s), BH-adjusted", lv[1], lv[2])
    test_levels <- lv[1:2]
    n_sig <- sum(padj < 0.05, na.rm = TRUE)
  } else {
    stats_df <- data.frame(feature = rownames(mat), mean = rowMeans(mat),
                           sd = apply(mat, 1, stats::sd), stringsAsFactors = FALSE)
    test_label <- NULL; test_levels <- NULL; n_sig <- NA_integer_
  }

  col_fun <- if (scale_mode == "none") {
    q <- stats::quantile(disp, c(0.02, 0.5, 0.98), na.rm = TRUE)
    circlize::colorRamp2(q, c("#3B6DB3", "white", "#C1432B"))
  } else {
    circlize::colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))
  }

  ht <- ComplexHeatmap::Heatmap(
    disp, name = legend_name, col = col_fun,
    cluster_rows = if (is.null(row_tree)) FALSE else stats::as.dendrogram(row_tree),
    cluster_columns = if (is.null(col_tree)) FALSE else stats::as.dendrogram(col_tree),
    top_annotation = top_anno, right_annotation = right_anno,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 9, fontface = "bold")))
  draw_fn <- function() ComplexHeatmap::draw(ht, merge_legends = TRUE)

  test_meta <- list(
    name = "Clustered heatmap (ComplexHeatmap)",
    scale = scale_mode, clustering = cluster,
    n_features = nrow(mat), n_samples = ncol(mat),
    per_feature_test = test_label, n_significant = n_sig,
    row_order = as.list(rownames(disp)[row_order]),
    column_order = as.list(colnames(disp)[col_order]))

  resolved <- list(family = "heatmap", label = "heatmap")
  methods  <- .hm_methods(scale_mode, cluster, nrow(mat), ncol(mat), test_levels, n_sig)
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

.hm_methods <- function(scale_mode, cluster, n_feat, n_samp, test_levels, n_sig) {
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
  if (cluster %in% c("both", "rows", "columns")) {
    clust_txt <- paste(clust_txt, "Clustering used complete linkage on the Euclidean distance of the displayed values, and at each merge the branch with the higher mean of the unscaled values was drawn first.")
  }
  sig_txt <- if (isTRUE(n_sig > 0)) {
    sprintf("%d of the %d features had an adjusted P below 0.05 and are marked with stars.", n_sig, n_feat)
  } else "no feature had an adjusted P below 0.05, so no star is drawn."
  test_txt <- if (!is.null(test_levels)) {
    sprintf(" Each feature was compared between %s and %s with a Welch two-sample t-test, which treats the samples as independent, and the P values were adjusted across the %d features with the Benjamini-Hochberg method; %s",
            test_levels[1], test_levels[2], n_feat, sig_txt)
  } else " No statistical test was applied."
  paste0(
    sprintf("A heatmap of %d features across %d samples was drawn with ComplexHeatmap. ",
            n_feat, n_samp),
    scale_txt, " ", clust_txt, test_txt,
    sprintf(" Rendered in R %s with ComplexHeatmap %s and circlize %s.",
            rver, v$ComplexHeatmap, v$circlize))
}

# The source of a function as script lines that define it under its name.
.hm_fn_lines <- function(name, fn) {
  src <- deparse(fn)
  c(paste(name, "<-", src[1]), src[-1], "")
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
  row_line <- if (cluster_rows) {
    'rows <- if (nrow(disp) > 1) as.dendrogram(.hm_order_tree(hclust(dist(disp), "complete"), rowMeans(mat))) else FALSE'
  } else "rows <- FALSE"
  col_line <- if (cluster_cols) {
    'cols <- if (ncol(disp) > 1) as.dendrogram(.hm_order_tree(hclust(dist(t(disp)), "complete"), colMeans(mat))) else FALSE'
  } else "cols <- FALSE"

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })",
    "",
    "# The leaf order of the figure: complete linkage on the Euclidean distance of",
    "# the drawn values, and at each merge the branch with the higher mean of the",
    "# unscaled values first. On a tie, the branch with the earlier input index.",
    .hm_fn_lines(".hm_branch", .hm_branch),
    .hm_fn_lines(".hm_order_tree", .hm_order_tree),
    .hm_fn_lines(".hm_leaf_order", .hm_leaf_order),
    sprintf('raw <- read.csv("%s", check.names = FALSE, row.names = 1)', in_name),
    "mat <- as.matrix(raw); storage.mode(mat) <- \"double\"",
    "keep <- apply(mat, 1, function(r) all(is.finite(r)) && sd(r) > 0); mat <- mat[keep, ]",
    disp_line,
    "disp[!is.finite(disp)] <- 0",
    ann_lines,
    sprintf('col_fun <- colorRamp2(c(-2, 0, 2), c("#3B6DB3", "white", "#C1432B"))'),
    row_line,
    col_line,
    sprintf('ht <- Heatmap(disp, name = "%s", col = col_fun, cluster_rows = rows, cluster_columns = cols, top_annotation = top)',
            legend_name),
    "",
    sprintf('cairo_pdf("%s.pdf", width = 7, height = 7)', fig_stub),
    "draw(ht, merge_legends = TRUE); dev.off()",
    'cat("Reproduced figure.\\n")'
  )
}
