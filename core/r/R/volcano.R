# Recipe: volcano plot of differential-expression results, EnhancedVolcano style.
# Input is a results table with a log2 fold-change column and a p-value column;
# common DESeq2 / limma / edgeR names are auto-detected, or set with --x / --y. A
# gene-label column is auto-detected or set with --label. Points are colored by
# whether they pass the fold-change cutoff, the p cutoff, both, or neither, and
# the top hits are labeled with ggrepel. No test is run here; the differential
# test is upstream, so the recipe is descriptive of that result.

.volcano_columns <- function(df, spec) {
  nm <- names(df)
  pick <- function(cands, override) {
    if (!is.null(override) && nzchar(override)) return(override)
    hit <- cands[cands %in% nm]
    if (length(hit)) hit[1] else NA_character_
  }
  list(
    fc  = pick(c("log2FoldChange", "logFC", "log2FC", "log2fc", "avg_log2FC"), spec$data$x),
    p   = pick(c("pvalue", "P.Value", "PValue", "p_value", "pval", "p.value"), spec$data$y),
    lab = {
      lab <- spec$data$label
      if (is.null(lab) || !nzchar(lab)) {
        chr <- nm[vapply(df, function(c) is.character(c) || is.factor(c), logical(1))]
        if (length(chr)) chr[1] else NA_character_
      } else lab
    })
}

recipe_volcano <- function(df, spec) {
  cols <- .volcano_columns(df, spec)
  if (is.na(cols$fc)) stop("volcano: no log2 fold-change column found; set it with --x")
  if (is.na(cols$p))  stop("volcano: no p-value column found; set it with --y")

  fc  <- suppressWarnings(as.numeric(df[[cols$fc]]))
  p   <- suppressWarnings(as.numeric(df[[cols$p]]))
  lab <- if (!is.na(cols$lab)) as.character(df[[cols$lab]]) else as.character(seq_len(nrow(df)))

  n0 <- nrow(df)
  keep <- is.finite(fc) & is.finite(p) & p > 0
  fc <- fc[keep]; p <- p[keep]; lab <- lab[keep]
  if (length(fc) < 1) stop("volcano: no usable rows after cleaning")
  clean_steps <- if (sum(!keep) > 0) {
    sprintf("Dropped %d row(s) with missing or non-finite fold change or p (%d -> %d).",
            sum(!keep), n0, length(fc))
  } else "No cleaning applied; input used as-is."
  df_used <- stats::setNames(data.frame(lab, fc, p, stringsAsFactors = FALSE),
                             c(cols$lab %||% "label", cols$fc, cols$p))

  fcc  <- spec$appearance$fc_cutoff %||% 1.0
  pcc  <- spec$appearance$p_cutoff %||% 0.05
  topn <- spec$appearance$top_n %||% 15

  neglog <- -log10(p)
  passfc <- abs(fc) >= fcc
  passp  <- p <= pcc
  cat <- ifelse(passfc & passp, "FC and p", ifelse(passp, "p", ifelse(passfc, "FC", "NS")))
  cat <- factor(cat, levels = c("NS", "FC", "p", "FC and p"))
  pal <- c(NS = "grey70", FC = "#2E8B57", p = "#3B6DB3", `FC and p` = "#C1432B")

  pdat <- data.frame(fc = fc, neglog = neglog, cat = cat, lab = lab, stringsAsFactors = FALSE)
  sig  <- pdat[pdat$cat == "FC and p", , drop = FALSE]
  sig  <- sig[order(-(abs(sig$fc) * sig$neglog)), , drop = FALSE]
  lab_df <- if (nrow(sig) > 0) utils::head(sig, topn) else
            utils::head(pdat[order(-pdat$neglog), , drop = FALSE], topn)

  n_up   <- sum(fc >= fcc & p <= pcc)
  n_down <- sum(fc <= -fcc & p <= pcc)

  p_plot <- ggplot2::ggplot(pdat, ggplot2::aes(x = fc, y = neglog, colour = cat)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
    ggplot2::geom_vline(xintercept = c(-fcc, fcc), linetype = "dashed",
                        linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(pcc), linetype = "dashed",
                        linewidth = 0.3, colour = "grey40") +
    ggrepel::geom_text_repel(data = lab_df, ggplot2::aes(label = lab),
                             size = 3, colour = "black", max.overlaps = 20,
                             min.segment.length = 0, seed = 1) +
    ggplot2::labs(
      x     = spec$appearance$x_label %||% "log2 fold change",
      y     = spec$appearance$y_label %||% "-log10 p",
      title = spec$appearance$title) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "bottom")

  stats_df <- data.frame(
    gene = lab_df$lab, log2FC = lab_df$fc, p_value = 10^(-lab_df$neglog),
    neg_log10_p = lab_df$neglog,
    direction = ifelse(lab_df$fc > 0, "up", "down"), stringsAsFactors = FALSE)

  test_meta <- list(
    name = "Volcano plot", fc_cutoff = fcc, p_cutoff = pcc,
    n_total = length(fc), n_up = n_up, n_down = n_down, n_significant = n_up + n_down)

  resolved <- list(family = "volcano", label = "volcano plot")
  methods  <- .volcano_methods(length(fc), fcc, pcc, n_up, n_down)
  build_script <- function(in_name, fig_stub) {
    .volcano_emit_script(in_name, fig_stub, cols, fcc, pcc, topn)
  }

  list(plot = p_plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps, width = 6.0, height = 5.6)
}

.volcano_methods <- function(n, fcc, pcc, n_up, n_down) {
  v    <- pkg_versions(c("ggplot2", "ggrepel"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  paste0(
    sprintf("Differential-expression results for %d features were shown as a volcano plot. ", n),
    sprintf("Features with |log2 fold change| >= %s and p <= %s were called significant; ", fcc, pcc),
    sprintf("%d were up-regulated and %d down-regulated. ", n_up, n_down),
    "The differential test itself was performed upstream; this figure summarises its output. ",
    sprintf("Drawn in R %s with ggplot2 %s and ggrepel %s.", rver, v$ggplot2, v$ggrepel))
}

.volcano_emit_script <- function(in_name, fig_stub, cols, fcc, pcc, topn) {
  lab_expr <- if (!is.na(cols$lab)) sprintf('as.character(raw[["%s"]])', cols$lab) else "as.character(seq_len(nrow(raw)))"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('fc <- as.numeric(raw[["%s"]]); p <- as.numeric(raw[["%s"]])', cols$fc, cols$p),
    sprintf('lab <- %s', lab_expr),
    "keep <- is.finite(fc) & is.finite(p) & p > 0; fc <- fc[keep]; p <- p[keep]; lab <- lab[keep]",
    sprintf('fcc <- %s; pcc <- %s', fcc, pcc),
    'cat <- ifelse(abs(fc) >= fcc & p <= pcc, "FC and p", ifelse(p <= pcc, "p", ifelse(abs(fc) >= fcc, "FC", "NS")))',
    'cat <- factor(cat, levels = c("NS", "FC", "p", "FC and p"))',
    'pal <- c(NS = "grey70", FC = "#2E8B57", p = "#3B6DB3", `FC and p` = "#C1432B")',
    "d <- data.frame(fc = fc, neglog = -log10(p), cat = cat, lab = lab)",
    'sig <- d[d$cat == "FC and p", ]; sig <- sig[order(-(abs(sig$fc) * sig$neglog)), ]',
    sprintf('lab_df <- head(sig, %d)', topn),
    "p <- ggplot(d, aes(fc, neglog, colour = cat)) + geom_point(size = 1.4, alpha = 0.7) +",
    "  scale_colour_manual(values = pal, drop = FALSE, name = NULL) +",
    "  geom_vline(xintercept = c(-fcc, fcc), linetype = \"dashed\") +",
    "  geom_hline(yintercept = -log10(pcc), linetype = \"dashed\") +",
    "  geom_text_repel(data = lab_df, aes(label = lab), size = 3, colour = \"black\", seed = 1) +",
    '  labs(x = "log2 fold change", y = "-log10 p") + theme_classic() + theme(legend.position = "bottom")',
    sprintf('ggsave("%s.pdf", p, width = 6, height = 5.6)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
