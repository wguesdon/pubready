# Recipe: MA plot of differential-expression results. M is the log2 fold change
# and A the mean expression of each gene. Input is a results table: DESeq2
# (log2FoldChange, padj, baseMean), limma (logFC, adj.P.Val, AveExpr) or edgeR
# (logFC, FDR, logCPM) names are detected; --x sets the mean column and --y the
# p column. DESeq2's baseMean is drawn on the log10 scale, the other means are
# already logarithmic. Genes that pass the p cutoff and the fold-change cutoff are
# red (up) or blue (down); the rest are grey and drawn first. No test is run
# here: the differential test is upstream, so the recipe describes its result.

.ma_fc_cols    <- c("log2FoldChange", "logFC", "log2FC", "log2fc", "avg_log2FC")
.ma_padj_cols  <- c("padj", "adj.P.Val", "FDR", "p_val_adj", "qvalue")
.ma_rawp_cols  <- c("pvalue", "P.Value", "PValue", "p_value", "pval", "p.value")
.ma_mean_cols  <- c("baseMean", "AveExpr", "logCPM")
.ma_mean_titles <- c(baseMean = "log10 mean of normalized counts",
                     AveExpr = "Average log2 expression", logCPM = "log2 CPM")

# Above this many grey points, the grey points that fall in the same cell of a
# fine grid are drawn once, so the vector outputs stay small.
.ma_thin_above <- 20000
.ma_grid_bins  <- 600

# A gene label longer than this is cut and ends in an ellipsis on the figure;
# the stats table keeps the full name.
.ma_label_max <- 18

.ma_short_label <- function(s, max_chars = .ma_label_max) {
  ifelse(nchar(s) > max_chars, paste0(substr(s, 1, max_chars - 1), "\u2026"), s)
}

.ma_first_hit <- function(cands, nm) {
  hit <- cands[cands %in% nm]
  if (length(hit)) hit[1] else NA_character_
}

.ma_is_text <- function(v) is.character(v) || is.factor(v)

.ma_columns <- function(df, spec) {
  nm <- names(df)
  p_col <- spec$data$y
  if (is.null(p_col) || !nzchar(p_col)) p_col <- .ma_first_hit(c(.ma_padj_cols, .ma_rawp_cols), nm)
  a_col <- spec$data$x
  if (is.null(a_col) || !nzchar(a_col)) a_col <- .ma_first_hit(.ma_mean_cols, nm)
  lab <- spec$data$label
  if (is.null(lab) || !nzchar(lab)) {
    txt <- nm[vapply(df, .ma_is_text, logical(1))]
    lab <- if (length(txt)) txt[1] else NA_character_
  }
  list(fc = .ma_first_hit(.ma_fc_cols, nm), p = p_col, a = a_col, lab = lab)
}

# TRUE for the first point in each cell of an nbin x nbin grid over the range.
.ma_grid_first <- function(x, y, nbin) {
  span_x <- max(diff(range(x)), 1e-12)
  span_y <- max(diff(range(y)), 1e-12)
  bx <- floor((x - min(x)) / span_x * nbin)
  by <- floor((y - min(y)) / span_y * nbin)
  !duplicated(bx * (nbin + 1) + by)
}

recipe_ma_plot <- function(df, spec) {
  cols <- .ma_columns(df, spec)
  if (is.na(cols$fc)) {
    stop(sprintf("ma_plot: no log2 fold-change column found; name it one of %s",
                 paste(.ma_fc_cols, collapse = ", ")))
  }
  if (is.na(cols$p)) stop("ma_plot: no p-value column found; set it with --y")
  if (is.na(cols$a)) {
    stop("ma_plot: no mean-expression column (baseMean, AveExpr, logCPM) found; set it with --x")
  }
  for (cn in c(cols$p, cols$a, if (!is.na(cols$lab)) cols$lab)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  fc <- suppressWarnings(as.numeric(df[[cols$fc]]))
  p  <- suppressWarnings(as.numeric(df[[cols$p]]))
  a_raw <- suppressWarnings(as.numeric(df[[cols$a]]))
  log_a <- identical(cols$a, "baseMean")
  a <- if (log_a) ifelse(is.finite(a_raw) & a_raw > 0, log10(pmax(a_raw, 1e-300)), NA_real_) else a_raw
  lab <- if (!is.na(cols$lab)) as.character(df[[cols$lab]]) else as.character(seq_len(nrow(df)))
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop(sprintf("ma_plot: column '%s' holds values outside 0 to 1", cols$p))

  n0 <- nrow(df)
  keep <- is.finite(fc) & is.finite(a)
  fc <- fc[keep]; p <- p[keep]; a <- a[keep]; lab <- lab[keep]
  if (length(fc) < 1) stop("ma_plot: no usable rows after cleaning")
  n_missing_p <- sum(is.na(p))
  clean_steps <- c(
    if (sum(!keep) > 0) {
      sprintf("Dropped %d row(s) with a missing or non-finite fold change or mean%s (%d -> %d).",
              sum(!keep), if (log_a) ", or a baseMean of zero" else "", n0, length(fc))
    },
    if (n_missing_p > 0) {
      sprintf("Kept %d row(s) with a missing p as not significant.", n_missing_p)
    },
    if (log_a) "Plotted log10(baseMean) as the mean expression A.")
  if (length(clean_steps) == 0) clean_steps <- "No cleaning applied; input used as-is."

  fcc  <- spec$appearance$fc_cutoff %||% 1.0
  pcc  <- spec$appearance$p_cutoff %||% 0.05
  topn <- as.integer(spec$appearance$top_n %||% 15)

  sig  <- !is.na(p) & p <= pcc
  up   <- sig & fc >= fcc
  down <- sig & fc <= -fcc
  cat  <- ifelse(up, "Up", ifelse(down, "Down", "NS"))
  n_up <- sum(up); n_down <- sum(down); n_ns <- length(fc) - n_up - n_down

  hit <- which(up | down)
  hit <- hit[order(p[hit], -abs(fc[hit]), hit)]
  stats_df <- data.frame(gene = lab[hit], A = a[hit], log2FC = fc[hit], p_value = p[hit],
                         direction = ifelse(fc[hit] > 0, "up", "down"), stringsAsFactors = FALSE)
  top <- utils::head(hit, topn)
  lab_df <- data.frame(A = a[top], M = fc[top], gene = .ma_short_label(lab[top]),
                       stringsAsFactors = FALSE)

  ns_idx <- which(cat == "NS")
  thinned <- n_ns > .ma_thin_above
  if (thinned) {
    ns_idx <- ns_idx[.ma_grid_first(a[ns_idx], fc[ns_idx], .ma_grid_bins)]
  }
  pdat <- data.frame(A = a, M = fc, cat = factor(cat, levels = c("NS", "Down", "Up")))
  ns_dat <- pdat[ns_idx, , drop = FALSE]
  # The most significant genes are drawn last, so they sit on top.
  sig_dat <- pdat[rev(hit), , drop = FALSE]

  plab <- if (cols$p %in% .ma_rawp_cols) "p" else cols$p
  a_title <- if (cols$a %in% names(.ma_mean_titles)) .ma_mean_titles[[cols$a]] else cols$a
  # The y range covers the data and both cutoff lines, with a 6% margin.
  ylo <- min(fc, -1.5 * fcc); yhi <- max(fc, 1.5 * fcc)
  ylim <- c(ylo, yhi) + c(-0.06, 0.06) * (yhi - ylo)
  pal <- c(NS = "grey70", Down = "#3B6DB3", Up = "#C1432B")
  sub <- sprintf("%d up / %d down (%s <= %s, |log2FC| >= %s)", n_up, n_down, plab,
                 format(pcc), format(fcc))

  p_plot <- ggplot2::ggplot(mapping = ggplot2::aes(x = A, y = M, colour = cat)) +
    ggplot2::geom_point(data = ns_dat, size = 0.45, stroke = 0, shape = 16) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45, colour = "grey15") +
    ggplot2::geom_hline(yintercept = c(-fcc, fcc), linetype = "dashed",
                        linewidth = 0.35, colour = "grey40") +
    ggplot2::geom_point(data = sig_dat, size = 1.2, alpha = 0.85, stroke = 0, shape = 16) +
    ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(x = A, y = M, label = gene),
                             inherit.aes = FALSE, size = if (nrow(lab_df) > 10) 2.6 else 3,
                             colour = "black",
                             max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey40", segment.size = 0.3,
                             box.padding = 0.35, point.padding = 0.15, force = 1.5,
                             max.time = 2, max.iter = 30000, seed = 1) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL,
                                 labels = c(NS = "Not significant", Down = "Down", Up = "Up")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    ggplot2::scale_y_continuous(limits = ylim, breaks = scales::breaks_extended(7)) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_extended(6)) +
    ggplot2::labs(x = spec$appearance$x_label %||% a_title,
                  y = spec$appearance$y_label %||% "log2 fold change",
                  title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "bottom",
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  test_meta <- list(
    name = "MA plot", fc_cutoff = fcc, p_cutoff = pcc, p_column = cols$p,
    mean_column = cols$a, mean_log10 = log_a,
    n_total = length(fc), n_up = n_up, n_down = n_down, n_not_significant = n_ns,
    n_missing_p = n_missing_p, n_grey_drawn = length(ns_idx))
  resolved <- list(family = "ma_plot", label = "MA plot")
  methods <- .ma_methods(test_meta, plab, a_title, nrow(lab_df), thinned)

  list(plot = p_plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods,
       build_script = list(fn = .ma_emit_script,
                           args = list(cols = cols, log_a = log_a, fcc = fcc, pcc = pcc,
                                       topn = topn, a_title = a_title)),
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = stats::setNames(data.frame(lab, a, fc, p, stringsAsFactors = FALSE),
                                 c(if (is.na(cols$lab)) "label" else cols$lab, "A", cols$fc, cols$p)),
       clean_steps = clean_steps, width = 6.0, height = 5.4)
}

.ma_methods <- function(tm, plab, a_title, n_labels, thinned) {
  v    <- pkg_versions(c("ggplot2", "ggrepel"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  paste0(
    sprintf("Differential-expression results for %d genes were shown as an MA plot of the log2 fold change (M) against the mean expression (A, %s). ",
            tm$n_total, a_title),
    sprintf("Genes with %s <= %s and |log2 fold change| >= %s were called significant: %d were up-regulated, %d down-regulated and %d not significant. ",
            plab, format(tm$p_cutoff), format(tm$fc_cutoff), tm$n_up, tm$n_down, tm$n_not_significant),
    if (tm$n_missing_p > 0) sprintf("%d genes without a p value were counted as not significant. ", tm$n_missing_p) else "",
    if (n_labels > 0) sprintf("The %d significant genes with the smallest p are labelled. ", n_labels) else "",
    if (thinned) sprintf("Grey points that fall in the same cell of a %d by %d grid were drawn once (%d of %d drawn); every count uses all genes. ",
                         .ma_grid_bins, .ma_grid_bins, tm$n_grey_drawn, tm$n_not_significant) else "",
    "The differential test itself was performed upstream; this figure summarises its output. ",
    sprintf("Drawn in R %s with ggplot2 %s and ggrepel %s.", rver, v$ggplot2, v$ggrepel))
}

.ma_emit_script <- function(cols, log_a, fcc, pcc, topn, a_title, in_name, fig_stub) {
  lab_expr <- if (!is.na(cols$lab)) sprintf('as.character(raw[["%s"]])', cols$lab) else
    "as.character(seq_len(nrow(raw)))"
  a_expr <- if (log_a) sprintf('log10(as.numeric(raw[["%s"]]))', cols$a) else
    sprintf('as.numeric(raw[["%s"]])', cols$a)
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('m <- as.numeric(raw[["%s"]]); p <- as.numeric(raw[["%s"]])', cols$fc, cols$p),
    sprintf("a <- %s; lab <- %s", a_expr, lab_expr),
    "keep <- is.finite(m) & is.finite(a); m <- m[keep]; p <- p[keep]; a <- a[keep]; lab <- lab[keep]",
    sprintf("fcc <- %s; pcc <- %s", fcc, pcc),
    "sig <- !is.na(p) & p <= pcc",
    'cat <- factor(ifelse(sig & m >= fcc, "Up", ifelse(sig & m <= -fcc, "Down", "NS")), levels = c("NS", "Down", "Up"))',
    "d <- data.frame(A = a, M = m, cat = cat, lab = lab)",
    'hits <- d[d$cat != "NS", ]; hits <- hits[order(p[d$cat != "NS"], -abs(hits$M)), ]',
    sprintf("top <- head(hits, %d)", topn),
    "lo <- min(m, -1.5 * fcc); hi <- max(m, 1.5 * fcc); lim <- c(lo, hi) + c(-0.06, 0.06) * (hi - lo)",
    "g <- ggplot(d, aes(A, M, colour = cat)) +",
    '  geom_point(data = d[d$cat == "NS", ], size = 0.45, shape = 16) +',
    '  geom_hline(yintercept = 0) + geom_hline(yintercept = c(-fcc, fcc), linetype = "dashed", colour = "grey40") +',
    '  geom_point(data = hits, size = 1.2, alpha = 0.85, shape = 16) +',
    '  geom_text_repel(data = top, aes(label = lab), colour = "black", size = 3, seed = 1, max.overlaps = Inf) +',
    '  scale_colour_manual(values = c(NS = "grey70", Down = "#3B6DB3", Up = "#C1432B"), drop = FALSE, name = NULL) +',
    sprintf('  scale_y_continuous(limits = lim) + labs(x = "%s", y = "log2 fold change") +', a_title),
    '  theme_classic(base_size = 13) + theme(legend.position = "bottom")',
    sprintf('ggsave("%s.pdf", g, width = 6, height = 5.4)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
