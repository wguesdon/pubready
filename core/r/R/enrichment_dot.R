# Recipe: dot plot of GSEA / ORA enrichment results, clusterProfiler style.
# Input is an enrichment table; common column names are auto-detected. When an
# NES column is present the plot is a GSEA dot plot (x = NES); otherwise an ORA
# dot plot (x = gene ratio). Dot color encodes the adjusted p, size the gene
# count. Descriptive: no test is run here, the enrichment is computed upstream.

.enrich_columns <- function(df, spec) {
  nm <- names(df)
  pick <- function(cands, override = NULL) {
    if (!is.null(override) && nzchar(override)) return(override)
    hit <- cands[cands %in% nm]; if (length(hit)) hit[1] else NA_character_
  }
  list(
    term  = pick(c("Description", "ID", "pathway", "term", "Term", "NAME"), spec$data$label),
    padj  = pick(c("p.adjust", "padj", "qvalue", "q_value", "FDR", "pvalue", "p_value")),
    count = pick(c("Count", "setSize", "count", "size")),
    ratio = pick(c("GeneRatio", "gene_ratio", "generatio")),
    nes   = pick(c("NES", "nes", "enrichmentScore", "ES")))
}

.parse_ratio <- function(x) {
  vapply(as.character(x), function(s) {
    if (grepl("/", s, fixed = TRUE)) {
      ab <- suppressWarnings(as.numeric(strsplit(s, "/", fixed = TRUE)[[1]]))
      if (length(ab) == 2 && is.finite(ab[2]) && ab[2] != 0) ab[1] / ab[2] else NA_real_
    } else suppressWarnings(as.numeric(s))
  }, numeric(1), USE.NAMES = FALSE)
}

recipe_enrichment_dot <- function(df, spec) {
  cols <- .enrich_columns(df, spec)
  if (is.na(cols$term)) stop("enrichment_dot: no term/Description column found; set it with --label")
  if (is.na(cols$padj)) stop("enrichment_dot: no adjusted-p column found")

  term  <- as.character(df[[cols$term]])
  padj  <- suppressWarnings(as.numeric(df[[cols$padj]]))
  count <- if (!is.na(cols$count)) suppressWarnings(as.numeric(df[[cols$count]])) else rep(NA_real_, length(term))
  gsea_mode <- !is.na(cols$nes)
  if (gsea_mode) {
    xval <- suppressWarnings(as.numeric(df[[cols$nes]])); xlab_default <- "NES"
  } else if (!is.na(cols$ratio)) {
    xval <- .parse_ratio(df[[cols$ratio]]); xlab_default <- "gene ratio"
  } else {
    xval <- -log10(padj); xlab_default <- "-log10 adjusted p"
  }

  keep <- is.finite(padj) & is.finite(xval) & !is.na(term)
  term <- term[keep]; padj <- padj[keep]; count <- count[keep]; xval <- xval[keep]
  if (length(term) < 1) stop("enrichment_dot: no usable rows")
  df_used <- data.frame(term = term, x = xval, p_adjust = padj, count = count, stringsAsFactors = FALSE)

  topn <- spec$appearance$top_n %||% 15
  idx  <- utils::head(order(padj), topn)
  pd <- data.frame(term = term[idx], x = xval[idx], padj = padj[idx],
                   count = count[idx], stringsAsFactors = FALSE)
  pd$term <- factor(pd$term, levels = pd$term[order(pd$x)])
  has_size <- !all(is.na(pd$count))

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = term))
  p <- if (has_size) {
    p + ggplot2::geom_point(ggplot2::aes(size = count, colour = padj))
  } else {
    p + ggplot2::geom_point(ggplot2::aes(colour = padj), size = 3)
  }
  p <- p +
    ggplot2::scale_colour_gradient(low = "#C1432B", high = "#3B6DB3", name = "p.adjust") +
    ggplot2::scale_size_continuous(name = "count") +
    ggplot2::labs(x = spec$appearance$x_label %||% xlab_default, y = NULL,
                  title = spec$appearance$title) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right", axis.title.y = ggplot2::element_blank())
  if (gsea_mode) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 linewidth = 0.3, colour = "grey50")
  }

  stats_df <- data.frame(term = as.character(pd$term), x = pd$x, p_adjust = pd$padj,
                         count = pd$count, stringsAsFactors = FALSE)
  names(stats_df)[2] <- if (gsea_mode) "NES" else if (!is.na(cols$ratio)) "gene_ratio" else "neg_log10_padj"

  test_meta <- list(name = if (gsea_mode) "GSEA dot plot" else "ORA dot plot",
                    n_terms_total = length(term), n_terms_shown = nrow(pd), x_axis = xlab_default)
  resolved <- list(family = "enrichment_dot",
                   label = if (gsea_mode) "GSEA dot plot" else "ORA dot plot")
  methods <- .enrich_methods(gsea_mode, length(term), nrow(pd), xlab_default)
  build_script <- function(in_name, fig_stub) .enrich_emit_script(in_name, fig_stub, cols, gsea_mode, topn)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = "No cleaning applied; enrichment table used as-is.",
       width = 7.2, height = max(4.0, 0.34 * nrow(pd) + 1.6))
}

.enrich_methods <- function(gsea_mode, n_total, n_shown, xlab) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  kind <- if (gsea_mode) "GSEA" else "over-representation (ORA)"
  paste0(
    sprintf("The top %d of %d enriched terms from a %s analysis were shown as a dot plot. ",
            n_shown, n_total, kind),
    sprintf("The x-axis is %s, dot color the adjusted p-value, and dot size the gene count. ", xlab),
    "The enrichment analysis was performed upstream; this figure summarises its output. ",
    sprintf("Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.enrich_emit_script <- function(in_name, fig_stub, cols, gsea_mode, topn) {
  xline <- if (gsea_mode) {
    sprintf('d$x <- as.numeric(raw[["%s"]]); xlab <- "NES"', cols$nes)
  } else if (!is.na(cols$ratio)) {
    sprintf('d$x <- sapply(strsplit(as.character(raw[["%s"]]), "/"), function(a) as.numeric(a[1])/as.numeric(a[2])); xlab <- "gene ratio"', cols$ratio)
  } else {
    sprintf('d$x <- -log10(as.numeric(raw[["%s"]])); xlab <- "-log10 adjusted p"', cols$padj)
  }
  vline <- if (gsea_mode) '  + geom_vline(xintercept = 0, linetype = "dashed")' else ""
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(term = as.character(raw[["%s"]]), padj = as.numeric(raw[["%s"]]))',
            cols$term, cols$padj),
    if (!is.na(cols$count)) sprintf('d$count <- as.numeric(raw[["%s"]])', cols$count) else "d$count <- NA",
    xline,
    sprintf('d <- d[order(d$padj), ][seq_len(min(%d, nrow(d))), ]', topn),
    "d$term <- factor(d$term, levels = d$term[order(d$x)])",
    'p <- ggplot(d, aes(x, term)) + geom_point(aes(size = count, colour = padj)) +',
    '  scale_colour_gradient(low = "#C1432B", high = "#3B6DB3") + labs(x = xlab, y = NULL) + theme_classic()',
    paste0("p <- p", vline),
    sprintf('ggsave("%s.pdf", p, width = 7.2, height = 6)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
