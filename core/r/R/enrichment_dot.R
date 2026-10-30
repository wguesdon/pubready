# Recipe: dot plot of GSEA / ORA enrichment results, clusterProfiler style.
# Input is an enrichment table; common column names are auto-detected. When an
# NES column is present the plot is a GSEA dot plot (x = NES); otherwise an ORA
# dot plot (x = gene ratio). Dot color encodes -log10 of the adjusted p, size the gene
# count. Descriptive: no test is run here, the enrichment is computed upstream.

.enrich_columns <- function(df, spec) {
  nm <- names(df)
  pick <- function(cands, override = NULL) {
    if (!is.null(override) && nzchar(override)) return(override)
    hit <- cands[cands %in% nm]; if (length(hit)) hit[1] else NA_character_
  }
  list(
    term  = pick(c("Description", "ID", "pathway", "term", "Term", "NAME"), spec$data$label),
    padj  = pick(c("p.adjust", "padj", "qvalue", "q_value", "FDR")),
    praw  = pick(c("pvalue", "p_value", "P.Value", "pval")),
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

# Fill words greedily into lines of at most `width` characters. A word longer
# than the width keeps its own line.
.enrich_greedy <- function(words, width) {
  lines <- words[1]
  for (wd in words[-1]) {
    last <- lines[length(lines)]
    if (nchar(last) + 1 + nchar(wd) <= width) {
      lines[length(lines)] <- paste(last, wd)
    } else {
      lines <- c(lines, wd)
    }
  }
  lines
}

# Wrap one term label at spaces into balanced lines of at most `width`
# characters: the fewest lines that the width allows, each about the same length,
# so no short word is left alone on a line. The Python engine wraps the same way.
.enrich_wrap <- function(s, width = 45) {
  words <- strsplit(trimws(s), "\\s+")[[1]]
  joined <- paste(words, collapse = " ")
  if (length(words) < 2 || nchar(joined) <= width) return(joined)
  n_lines <- length(.enrich_greedy(words, width))
  for (w2 in seq(min(ceiling(nchar(joined) / n_lines), width), width)) {
    lines <- .enrich_greedy(words, w2)
    if (length(lines) <= n_lines) return(paste(lines, collapse = "\n"))
  }
  paste(.enrich_greedy(words, width), collapse = "\n")
}

# The figure size from the wrapped labels: each term row grows with its lines,
# and the label column grows with its longest line. The Python engine uses the
# same sizes.
.enrich_size <- function(labels) {
  parts <- strsplit(labels, "\n", fixed = TRUE)
  n_lines <- lengths(parts)
  longest <- max(nchar(unlist(parts)), 4)
  list(width = max(7.2, 0.3 + 0.075 * longest + 3.2 + 1.6),
       height = max(4.0, 1.6 + sum(0.2 + 0.16 * n_lines)))
}

recipe_enrichment_dot <- function(df, spec) {
  cols <- .enrich_columns(df, spec)
  if (is.na(cols$term)) stop("enrichment_dot: no term/Description column found; set it with --label")
  # A table with only a raw p is drawn and labelled as not adjusted, so a raw p
  # is never shown as an adjusted one (F52).
  adjusted <- !is.na(cols$padj)
  if (!adjusted) cols$padj <- cols$praw
  if (is.na(cols$padj)) {
    stop("enrichment_dot: no p-value column found; give p.adjust, padj, qvalue, FDR or pvalue")
  }
  pw <- .enrich_p_words(adjusted)

  term  <- as.character(df[[cols$term]])
  padj  <- suppressWarnings(as.numeric(df[[cols$padj]]))
  count <- if (!is.na(cols$count)) suppressWarnings(as.numeric(df[[cols$count]])) else rep(NA_real_, length(term))
  # An adjusted p of exactly 0 is the strongest term, not a missing value. It is
  # set to the smallest positive double, as manhattan does, and drawn as a
  # triangle at the top of the colour scale (F51).
  is_zero <- is.finite(padj) & padj == 0
  padj[is_zero] <- .Machine$double.xmin
  gsea_mode <- !is.na(cols$nes)
  neglog_x <- !gsea_mode && is.na(cols$ratio)
  if (gsea_mode) {
    xval <- suppressWarnings(as.numeric(df[[cols$nes]])); xlab_default <- "NES"
  } else if (!is.na(cols$ratio)) {
    xval <- .parse_ratio(df[[cols$ratio]]); xlab_default <- "gene ratio"
  } else {
    xval <- -log10(padj); xlab_default <- pw$axis
  }

  n0 <- length(term)
  keep <- is.finite(padj) & is.finite(xval) & !is.na(term)
  term <- term[keep]; padj <- padj[keep]; count <- count[keep]; xval <- xval[keep]
  is_zero <- is_zero[keep]
  if (length(term) < 1) stop("enrichment_dot: no usable rows")
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing term, %s or x value (%d -> %d).",
                                sum(!keep), pw$noun, n0, length(term)),
    if (any(is_zero)) sprintf("Replaced %d %s value(s) of exactly 0 with the smallest positive double (%s), so that -log10 p is finite.",
                              sum(is_zero), pw$noun, formatC(.Machine$double.xmin, format = "g", digits = 3)))
  if (!adjusted) {
    clean_steps <- c(sprintf("The table holds no adjusted p column, so the raw p column '%s' is shown and labelled as not adjusted.",
                             cols$padj), clean_steps)
  }
  if (length(clean_steps) == 0) clean_steps <- "No cleaning applied; enrichment table used as-is."
  df_used <- data.frame(term = term, x = xval, p_adjust = padj, count = count, stringsAsFactors = FALSE)
  # The colour, and the x of the -log10 p axis, of a p of 0 sits at a ceiling 12%
  # above the largest other value, as in manhattan and volcano.
  shade <- .volcano_display_y(-log10(padj), is_zero, 0.05)
  if (neglog_x) xval <- shade

  topn <- spec$appearance$top_n %||% 15
  idx  <- utils::head(order(padj), topn)
  pd <- data.frame(term = term[idx], x = xval[idx], padj = padj[idx],
                   count = count[idx], shade = shade[idx],
                   shape = ifelse(is_zero[idx], 17, 16), stringsAsFactors = FALSE)
  n_zero_shown <- sum(is_zero[idx])
  pd$term <- factor(pd$term, levels = pd$term[order(pd$x)])
  # Long GO term names would take the whole width and leave no room for the dots,
  # so each label is wrapped at about 45 characters.
  pd$label <- vapply(as.character(pd$term), .enrich_wrap, character(1), USE.NAMES = FALSE)
  pd$label <- factor(pd$label, levels = pd$label[order(pd$x)])
  size <- .enrich_size(as.character(pd$label))
  has_size <- !all(is.na(pd$count))

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = label))
  p <- if (has_size) {
    p + ggplot2::geom_point(ggplot2::aes(size = count, colour = shade, shape = shape))
  } else {
    p + ggplot2::geom_point(ggplot2::aes(colour = shade, shape = shape), size = 3)
  }
  subtitle <- if (n_zero_shown > 0) {
    sprintf("%d term(s) with %s = 0 drawn as triangles", n_zero_shown, pw$noun)
  }
  p <- p +
    # Colour by -log10 of the adjusted p: the values span orders of magnitude, so
    # a linear scale would paint almost every term the same colour.
    ggplot2::scale_colour_gradient(low = "#3B6DB3", high = "#C1432B", name = pw$legend) +
    ggplot2::scale_size_continuous(name = "count") +
    ggplot2::scale_shape_identity() +
    ggplot2::labs(x = spec$appearance$x_label %||% xlab_default, y = NULL,
                  title = spec$appearance$title, subtitle = subtitle) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right", legend.justification = "top",
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(lineheight = 0.9),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))
  if (gsea_mode) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 linewidth = 0.3, colour = "grey50")
  }

  stats_df <- data.frame(term = as.character(pd$term), x = pd$x, p_adjust = pd$padj,
                         count = pd$count, stringsAsFactors = FALSE)
  names(stats_df)[2] <- if (gsea_mode) "NES" else if (!is.na(cols$ratio)) "gene_ratio" else "neg_log10_padj"
  names(stats_df)[3] <- pw$col

  test_meta <- list(name = if (gsea_mode) "GSEA dot plot" else "ORA dot plot",
                    n_terms_total = length(term), n_terms_shown = nrow(pd), x_axis = xlab_default,
                    n_zero_p = sum(is_zero), p_column = cols$padj, p_adjusted = adjusted)
  resolved <- list(family = "enrichment_dot",
                   label = if (gsea_mode) "GSEA dot plot" else "ORA dot plot")
  methods <- .enrich_methods(gsea_mode, length(term), nrow(pd), xlab_default, sum(is_zero),
                             adjusted)
  build_script <- list(fn = .enrich_emit_script,
                       args = list(cols = cols, gsea_mode = gsea_mode, topn = topn,
                                   width = size$width, height = size$height,
                                   adjusted = adjusted))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       width = size$width, height = size$height)
}

# The words for the p column: an adjusted p, or a raw p that the figure, the
# stats table and the methods text call not adjusted (F52). The Python engine
# uses the same words.
.enrich_p_words <- function(adjusted) {
  if (adjusted) {
    list(axis = "-log10 adjusted p", legend = "-log10 p.adjust", noun = "adjusted p", one = "an adjusted p",
         col = "p_adjust", methods = "the adjusted p-value")
  } else {
    list(axis = "-log10 p (not adjusted)", legend = "-log10 p (not adjusted)", noun = "p", one = "a p",
         col = "p_value",
         methods = "the raw p-value, which is not adjusted for multiple testing, because the table holds no adjusted p-value")
  }
}

.enrich_methods <- function(gsea_mode, n_total, n_shown, xlab, n_zero = 0, adjusted = TRUE) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  kind <- if (gsea_mode) "a GSEA" else "an over-representation (ORA)"
  pw <- .enrich_p_words(adjusted)
  zero_txt <- if (n_zero > 0) {
    sprintf("%d term(s) with %s of exactly 0 were set to the smallest positive double (%s) and are drawn as triangles at the top of the colour scale. ",
            n_zero, pw$one, formatC(.Machine$double.xmin, format = "g", digits = 3))
  } else ""
  paste0(
    sprintf("The top %d of %d enriched terms from %s analysis were shown as a dot plot. ",
            n_shown, n_total, kind),
    sprintf(paste("The x-axis is %s, dot color the -log10 of %s, and dot size",
                  "the gene count. "), xlab, pw$methods),
    zero_txt,
    "The enrichment analysis was performed upstream; this figure summarises its output. ",
    sprintf("Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.enrich_emit_script <- function(cols, gsea_mode, topn, width, height, in_name, fig_stub,
                                adjusted = TRUE) {
  pw <- .enrich_p_words(adjusted)
  xline <- if (gsea_mode) {
    sprintf('d$x <- as.numeric(raw[["%s"]]); xlab <- "NES"', cols$nes)
  } else if (!is.na(cols$ratio)) {
    sprintf('d$x <- sapply(strsplit(as.character(raw[["%s"]]), "/"), function(a) as.numeric(a[1])/as.numeric(a[2])); xlab <- "gene ratio"', cols$ratio)
  } else {
    sprintf('d$x <- -log10(d$padj); xlab <- "%s"', pw$axis)
  }
  vline <- if (gsea_mode) '  + geom_vline(xintercept = 0, linetype = "dashed")' else ""
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(term = as.character(raw[["%s"]]), padj = as.numeric(raw[["%s"]]))',
            cols$term, cols$padj),
    "# An adjusted p of exactly 0 is set to the smallest positive double, as in the figure.",
    "d$padj[!is.na(d$padj) & d$padj == 0] <- .Machine$double.xmin",
    if (!is.na(cols$count)) sprintf('d$count <- as.numeric(raw[["%s"]])', cols$count) else "d$count <- NA",
    xline,
    sprintf('d <- d[order(d$padj), ][seq_len(min(%d, nrow(d))), ]', topn),
    "# Wrap each term at about 45 characters, as the figure does.",
    "d$label <- d$term",
    "for (i in seq_len(nrow(d))) d$label[i] <- paste(strwrap(d$term[i], 46), collapse = \"\\n\")",
    "d$label <- factor(d$label, levels = d$label[order(d$x)])",
    'p <- ggplot(d, aes(x, label)) + geom_point(aes(size = count, colour = -log10(padj))) +',
    sprintf('  scale_colour_gradient(low = "#3B6DB3", high = "#C1432B", name = "%s") + labs(x = xlab, y = NULL) + theme_classic()', pw$legend),
    paste0("p <- p", vline),
    sprintf('ggsave("%s.pdf", p, width = %.2f, height = %.2f)', fig_stub, width, height),
    'cat("Reproduced figure.\\n")'
  )
}
