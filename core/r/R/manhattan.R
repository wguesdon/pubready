# Recipe: Manhattan plot and QQ plot of GWAS summary statistics. Input is one
# row per SNP: the chromosome (--group), the base-pair position (--x), the p
# value (--y) and, optionally, the SNP name (--label). The chromosomes are put in
# natural order (a "chr" prefix is dropped; 1..22, then X, Y, XY, MT), and each
# SNP is drawn at its cumulative genomic position with a small gap between the
# chromosomes. The genome-wide line is at 5e-8 and the suggestive line at 1e-5;
# the lead SNP of each chromosome that passes 5e-8 is labelled. The QQ panel
# carries a 95% band from the beta distribution of the uniform order statistics
# and the genomic inflation factor lambda_GC. Every statistic uses all SNPs; for
# the size of the vector outputs, the SNPs with p >= 0.01 are thinned
# deterministically when there are more than 50,000.

.mh_gw         <- 5e-8
.mh_sugg       <- 1e-5
.mh_max_points <- 50000
.mh_keep_below <- 0.01
.mh_min_rest   <- 10000
.mh_gap_frac   <- 0.004
.mh_label_sep  <- 0.03
.mh_special    <- c("X", "Y", "XY", "MT", "M")
.mh_blues      <- c("#3B6DB3", "#1F3A5F")

.mh_chr_name <- function(v) sub("^chr", "", trimws(as.character(v)), ignore.case = TRUE)

# The chromosomes in natural order: the numbers, then X, Y, XY, MT, then any
# other name in C-locale order.
.mh_chr_order <- function(chrs) {
  u <- unique(chrs)
  is_num <- grepl("^[0-9]+$", u)
  sp <- match(toupper(u), .mh_special)
  tier <- ifelse(is_num, 1L, ifelse(is.na(sp), 3L, 2L))
  key <- ifelse(is_num, suppressWarnings(as.numeric(u)), ifelse(is.na(sp), 0, sp))
  u[order(tier, key, ifelse(tier == 3L, u, ""), method = "radix")]
}

# Keep every point with p below keep_below, and every k-th of the rest in the
# given order, so that about max_points are drawn. k = 1 keeps all.
.mh_thin <- function(p, max_points = .mh_max_points, keep_below = .mh_keep_below,
                     min_rest = .mh_min_rest) {
  n <- length(p)
  if (n <= max_points) return(list(keep = rep(TRUE, n), k = 1L))
  low <- p < keep_below
  rest <- which(!low)
  k <- as.integer(ceiling(length(rest) / max(max_points - sum(low), min_rest)))
  keep <- low
  keep[rest[(seq_along(rest) - 1L) %% k == 0L]] <- TRUE
  list(keep = keep, k = k)
}

# Keep a chromosome label only when its centre is far enough from the last kept
# label, so that the labels of the small chromosomes do not collide.
.mh_label_keep <- function(centres, min_sep) {
  keep <- logical(length(centres))
  last <- -Inf
  for (i in seq_along(centres)) {
    if (centres[i] - last >= min_sep) {
      keep[i] <- TRUE
      last <- centres[i]
    }
  }
  keep
}

.mh_fmt_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

recipe_manhattan <- function(df, spec) {
  chr_col <- spec$data$group; pos_col <- spec$data$x; p_col <- spec$data$y
  lab_col <- spec$data$label
  if (is.null(chr_col) || !nzchar(chr_col)) stop("manhattan needs --group (the chromosome column)")
  if (is.null(pos_col) || !nzchar(pos_col)) stop("manhattan needs --x (the base-pair position column)")
  if (is.null(p_col) || !nzchar(p_col)) stop("manhattan needs --y (the p-value column)")
  has_lab <- !is.null(lab_col) && nzchar(lab_col)
  for (cn in c(chr_col, pos_col, p_col, if (has_lab) lab_col)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  raw_chr <- trimws(as.character(df[[chr_col]]))
  chr <- .mh_chr_name(raw_chr)
  pos <- suppressWarnings(as.numeric(df[[pos_col]]))
  p <- suppressWarnings(as.numeric(df[[p_col]]))
  bad <- sum(p < 0 | p > 1, na.rm = TRUE)
  if (bad > 0) {
    stop(sprintf("manhattan: p must lie in (0, 1]; column '%s' holds %d value(s) outside 0 to 1",
                 p_col, bad))
  }
  keep <- !is.na(chr) & nzchar(chr) & is.finite(pos) & is.finite(p)
  had_prefix <- any(grepl("^chr", raw_chr[keep], ignore.case = TRUE))
  chr <- chr[keep]; pos <- pos[keep]; p <- p[keep]
  snp <- if (has_lab) as.character(df[[lab_col]])[keep] else paste0(chr, ":", sprintf("%.0f", pos))
  if (length(p) < 2) stop("manhattan needs at least 2 SNPs with a chromosome, a position and a p value")
  is_zero <- p == 0
  n_zero <- sum(is_zero)
  p[is_zero] <- .Machine$double.xmin
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing chromosome, position or p (%d -> %d).",
                                sum(!keep), n0, length(p)),
    if (had_prefix) "Read the chromosome names without their 'chr' prefix.",
    if (n_zero > 0) sprintf("Replaced %d p value(s) of exactly 0 with the smallest positive double (%s), so that -log10 p is finite.",
                            n_zero, formatC(.Machine$double.xmin, format = "g", digits = 3)))
  if (length(clean_steps) == 0) clean_steps <- "No cleaning applied; input used as-is."

  # Statistics on every SNP.
  n <- length(p)
  lambda <- stats::median(stats::qchisq(p, 1, lower.tail = FALSE)) / stats::qchisq(0.5, 1)
  n_gw <- sum(p < .mh_gw); n_sugg <- sum(p < .mh_sugg)
  topn <- as.integer(spec$appearance$top_n %||% 15)

  lv <- .mh_chr_order(chr)
  ci <- match(chr, lv)
  o <- order(ci, p, seq_len(n))
  first <- o[!duplicated(ci[o])]
  leads <- first[p[first] < .mh_gw]
  leads <- utils::head(leads[order(p[leads], ci[leads], leads)], topn)

  # Cumulative position: each chromosome starts after the previous one plus a gap.
  lo <- as.numeric(tapply(pos, ci, min)); hi <- as.numeric(tapply(pos, ci, max))
  span <- hi - lo
  gap <- .mh_gap_frac * max(sum(span), 1)
  offset <- c(0, cumsum(span + gap))[seq_along(lv)]
  cum <- pos - lo[ci] + offset[ci]
  centre <- offset + span / 2
  xr <- c(min(cum), max(cum))
  lab_keep <- .mh_label_keep(centre, .mh_label_sep * max(diff(xr), 1))

  # A p of exactly 0 has no finite -log10 p of its own; it is drawn as a
  # triangle at a ceiling just above the largest other value.
  thin <- .mh_thin(p)
  draw <- thin$keep
  ydisp <- .mh_display_y(p, is_zero)
  ymax <- max(ydisp, -log10(.mh_gw)) * (if (length(leads)) 1.14 else 1.06)
  pt_size <- if (sum(draw) > 20000) 0.6 else 0.9
  sel <- draw & !is_zero
  mdat <- data.frame(x = cum[sel], y = ydisp[sel],
                     shade = factor((ci[sel] - 1L) %% 2L, levels = c(0L, 1L)))
  zdat <- data.frame(x = cum[is_zero], y = ydisp[is_zero])
  ldat <- data.frame(x = cum[leads], y = ydisp[leads], lab = .mh_short_label(snp[leads]),
                     shape = ifelse(is_zero[leads], 17, 16), stringsAsFactors = FALSE)
  sub <- sprintf("%s SNPs; %s with p < 5e-8; %s with p < 1e-5%s", .mh_fmt_count(n),
                 .mh_fmt_count(n_gw), .mh_fmt_count(n_sugg),
                 if (n_zero > 0) sprintf("; %d with p = 0 (triangles)", n_zero) else "")

  man <- ggplot2::ggplot(mdat, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(colour = shade), size = pt_size, stroke = 0, shape = 16) +
    ggplot2::geom_hline(yintercept = -log10(.mh_sugg), linetype = "dashed",
                        linewidth = 0.4, colour = "grey55") +
    ggplot2::geom_hline(yintercept = -log10(.mh_gw), linetype = "dashed",
                        linewidth = 0.45, colour = "#C1432B") +
    ggplot2::geom_point(data = zdat, colour = "#C1432B", size = 1.8, shape = 17) +
    ggplot2::geom_point(data = ldat, ggplot2::aes(shape = shape), colour = "#C1432B", size = 1.9) +
    ggplot2::scale_shape_identity() +
    ggrepel::geom_text_repel(data = ldat, ggplot2::aes(label = lab), size = 2.8, colour = "black",
                             max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                             segment.colour = "grey40", box.padding = 0.4, nudge_y = 0.4,
                             direction = "both", seed = 1) +
    ggplot2::scale_colour_manual(values = .mh_blues, guide = "none") +
    ggplot2::scale_x_continuous(breaks = centre[lab_keep], labels = lv[lab_keep],
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_extended(6),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::coord_cartesian(ylim = c(0, ymax)) +
    ggplot2::labs(x = spec$appearance$x_label %||% "Chromosome",
                  y = spec$appearance$y_label %||% "-log10 p",
                  title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  qq <- .mh_qq_data(p, is_zero)
  qq_sub <- bquote(lambda[GC] == .(formatC(lambda, format = "f", digits = 3)))
  qlim <- max(qq$exp_max, max(qq$pts$obs)) * 1.04
  qqp <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = qq$band, ggplot2::aes(x = exp, ymin = lo, ymax = hi),
                         fill = "grey85", colour = NA) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#C1432B", linewidth = 0.45) +
    ggplot2::geom_point(data = qq$pts[!qq$pts$zero, , drop = FALSE], ggplot2::aes(x = exp, y = obs),
                        colour = .mh_blues[1], size = pt_size, stroke = 0, shape = 16) +
    ggplot2::geom_point(data = qq$pts[qq$pts$zero, , drop = FALSE], ggplot2::aes(x = exp, y = obs),
                        colour = "#C1432B", size = 1.8, shape = 17) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_extended(5),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_extended(6),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::coord_cartesian(xlim = c(0, qq$exp_max * 1.04), ylim = c(0, qlim)) +
    ggplot2::labs(x = "Expected -log10 p", y = "Observed -log10 p", subtitle = qq_sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  plot <- ggpubr::ggarrange(man, qqp, ncol = 2, widths = c(2.7, 1), align = "h")

  lead_tab <- data.frame(row_type = "lead", chr = lv[ci[leads]], snp = snp[leads], pos = pos[leads],
                         p_value = p[leads], n_snps = NA_integer_, lambda_gc = NA_real_,
                         n_genome_wide = NA_integer_, n_suggestive = NA_integer_,
                         stringsAsFactors = FALSE)
  genome <- data.frame(row_type = "genome", chr = NA_character_, snp = NA_character_, pos = NA_real_,
                       p_value = min(p), n_snps = n, lambda_gc = lambda, n_genome_wide = n_gw,
                       n_suggestive = n_sugg, stringsAsFactors = FALSE)
  stats_df <- rbind(genome, lead_tab)

  test_meta <- list(
    name = "Manhattan and QQ plot", n_snps = n, n_chromosomes = length(lv),
    lambda_gc = lambda, n_genome_wide = n_gw, n_suggestive = n_sugg, min_p = min(p),
    genome_wide_threshold = .mh_gw, suggestive_threshold = .mh_sugg, n_zero_p = n_zero,
    n_leads = length(leads), thinning_step = thin$k, n_drawn = sum(draw))
  resolved <- list(family = "manhattan", label = test_meta$name)
  df_used <- data.frame(chr = chr, pos = pos, p = p, stringsAsFactors = FALSE)
  if (has_lab) df_used$snp <- snp
  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .mh_methods(test_meta),
       build_script = list(fn = .mh_emit_script,
                           args = list(chr_col = chr_col, pos_col = pos_col, p_col = p_col)),
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps, width = 10.5, height = 4.3)
}

# The QQ points and the 95% band of the uniform order statistics,
# Beta(i, n - i + 1), at the drawn points. The sorted points lie on one curve, so
# every point with p < 0.01 is drawn and, for the others, the first point in each
# cell of a grid of .mh_qq_cell in -log10 p on both axes.
.mh_qq_cell <- 0.005

.mh_qq_data <- function(p, is_zero) {
  n <- length(p)
  o <- order(p, method = "radix")
  ps <- p[o]
  zs <- is_zero[o]
  obs <- .mh_display_y(ps, zs)
  i <- seq_len(n)
  expd <- -log10((i - 0.5) / n)
  key <- floor(expd / .mh_qq_cell) * 1e6 + floor(obs / .mh_qq_cell)
  j <- i[ps < .mh_keep_below | !duplicated(key)]
  band <- data.frame(exp = expd[j],
                     lo = -log10(stats::qbeta(0.975, j, n - j + 1)),
                     hi = -log10(stats::qbeta(0.025, j, n - j + 1)))
  list(pts = data.frame(exp = expd[j], obs = obs[j], zero = zs[j]), band = band,
       exp_max = max(expd))
}

# -log10 p for drawing: a p of exactly 0 sits at a ceiling 12% above the
# largest other value (or above the genome-wide line).
.mh_display_y <- function(p, is_zero) {
  y <- -log10(p)
  top <- if (any(!is_zero)) max(y[!is_zero]) else 0
  y[is_zero] <- max(top, -log10(.mh_gw)) * 1.12
  y
}

# A SNP label longer than 20 characters is cut and ends in an ellipsis on the
# figure; the stats table keeps the full name.
.mh_short_label <- function(s) {
  ifelse(nchar(s) > 20, paste0(substr(s, 1, 19), "\u2026"), s)
}

.mh_methods <- function(tm) {
  v <- pkg_versions(c("ggplot2", "ggrepel", "ggpubr"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  thin_txt <- if (tm$thinning_step > 1) {
    sprintf("For the size of the vector outputs, every SNP with p < 0.01 is drawn and one in every %d of the other SNPs in input order (%s of %s points drawn); every statistic uses all SNPs. ",
            tm$thinning_step, .mh_fmt_count(tm$n_drawn), .mh_fmt_count(tm$n_snps))
  } else sprintf("All %s SNPs are drawn in the Manhattan panel. ", .mh_fmt_count(tm$n_snps))
  thin_txt <- paste0(thin_txt, "The QQ panel draws every SNP with p < 0.01 and, for the others, one point per 0.005 by 0.005 cell of -log10 p, because the sorted points lie on one curve. ")
  paste0(
    sprintf("Association results for %s SNPs on %d chromosomes were shown as a Manhattan plot of -log10 p at the cumulative genomic position, with the genome-wide significance line at p = 5e-8 and the suggestive line at p = 1e-5. ",
            .mh_fmt_count(tm$n_snps), tm$n_chromosomes),
    sprintf("%s SNPs reached genome-wide significance and %s passed the suggestive line (smallest p = %s); the lead SNP of each chromosome with p < 5e-8 is labelled (%d shown). ",
            .mh_fmt_count(tm$n_genome_wide), .mh_fmt_count(tm$n_suggestive),
            formatC(tm$min_p, format = "g", digits = 3), tm$n_leads),
    if (tm$n_zero_p > 0) sprintf("%d p value(s) of exactly 0 were set to the smallest positive double for the statistics and are drawn as triangles at the top of both panels. ", tm$n_zero_p) else "",
    sprintf("The QQ plot compares the observed -log10 p with the expected uniform quantiles, (i - 0.5) / n, with a 95%% band from the Beta(i, n - i + 1) distribution of the order statistics. The genomic inflation factor was lambda_GC = %s, the median of the 1-df chi-square statistics over its expected value of 0.455. ",
            formatC(tm$lambda_gc, format = "f", digits = 3)),
    thin_txt,
    "The association test itself was performed upstream; this figure summarises its output. ",
    sprintf("Drawn in R %s with ggplot2 %s, ggrepel %s and ggpubr %s.", rver, v$ggplot2, v$ggrepel, v$ggpubr))
}

.mh_emit_script <- function(chr_col, pos_col, p_col, in_name, fig_stub) {
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('chr <- sub("^chr", "", as.character(raw[["%s"]]), ignore.case = TRUE)', chr_col),
    sprintf('pos <- as.numeric(raw[["%s"]]); p <- as.numeric(raw[["%s"]])', pos_col, p_col),
    "keep <- !is.na(chr) & is.finite(pos) & is.finite(p); chr <- chr[keep]; pos <- pos[keep]; p <- p[keep]",
    "p[p == 0] <- .Machine$double.xmin",
    'num <- grepl("^[0-9]+$", unique(chr)); sp <- match(toupper(unique(chr)), c("X", "Y", "XY", "MT", "M"))',
    "lv <- unique(chr)[order(ifelse(num, 1, ifelse(is.na(sp), 3, 2)), ifelse(num, suppressWarnings(as.numeric(unique(chr))), ifelse(is.na(sp), 0, sp)))]",
    "ci <- match(chr, lv); lo <- tapply(pos, ci, min); span <- tapply(pos, ci, max) - lo",
    "offset <- c(0, cumsum(span + 0.004 * sum(span)))[seq_along(lv)]",
    "d <- data.frame(x = pos - lo[ci] + offset[ci], y = -log10(p), shade = factor(ci %% 2))",
    'lambda <- median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)',
    "man <- ggplot(d, aes(x, y, colour = shade)) + geom_point(size = 0.7, shape = 16) +",
    '  scale_colour_manual(values = c("#3B6DB3", "#1F3A5F"), guide = "none") +',
    '  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", colour = "grey55") +',
    '  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#C1432B") +',
    '  scale_x_continuous(breaks = offset + span / 2, labels = lv) +',
    '  labs(x = "Chromosome", y = "-log10 p") + theme_classic(base_size = 13)',
    "n <- length(p); q <- data.frame(exp = -log10((seq_len(n) - 0.5) / n), obs = -log10(sort(p)))",
    "q$lo <- -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1)); q$hi <- -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1))",
    "qqp <- ggplot(q, aes(exp, obs)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = \"grey85\") +",
    '  geom_abline(colour = "#C1432B") + geom_point(size = 0.7, colour = "#3B6DB3") +',
    '  labs(x = "Expected -log10 p", y = "Observed -log10 p", subtitle = sprintf("lambda_GC = %.3f", lambda)) +',
    "  theme_classic(base_size = 13)",
    "g <- ggarrange(man, qqp, ncol = 2, widths = c(2.7, 1), align = \"h\")",
    sprintf('ggsave("%s.pdf", g, width = 10.5, height = 4.3)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
