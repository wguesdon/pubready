# Recipe: the shape of one numeric variable (--y), overall or per group (--x).
# The figure is a histogram with a Gaussian kernel density (solid) and a normal
# curve with the sample mean and SD (dashed), both scaled to counts so that they
# overlay the bars, and a rug when a panel holds 500 values or fewer. Groups are
# drawn as panels stacked on a shared x axis. The stats table gives the moments,
# the quantiles and the Shapiro-Wilk test of each group. No comparison test is
# run: two_group_compare and multi_group_compare do that.
#
# The bins are the same in both engines: the Freedman-Diaconis width on the
# pooled data, 2 * IQR * n^(-1/3) with type-7 quantiles (Sturges when the IQR is
# zero), anchored at floor(min / width) * width and closed on the left. The
# density uses the bandwidth of stats::bw.nrd0 (Silverman's rule of thumb) and is
# evaluated directly on 512 points from min - 3 bw to max + 3 bw.

.dist_max_bins <- 200
.dist_grid_n <- 512
.dist_rug_max <- 500

# Type-7 quantile, the default of stats::quantile. The Python engine restates
# the same arithmetic, so the bins match to the last bit.
.dist_q7 <- function(v, p) stats::quantile(v, p, type = 7, names = FALSE)

.dist_null_on_error <- function(e) NULL

# Freedman-Diaconis bins on the pooled values, with the Sturges rule when the
# IQR is zero. A long-tailed variable can ask for thousands of bins; the width
# then grows so that 200 bins cover the range.
.dist_bins <- function(v) {
  n <- length(v)
  lo <- min(v)
  hi <- max(v)
  iqr <- .dist_q7(v, 0.75) - .dist_q7(v, 0.25)
  if (iqr > 0) {
    width <- 2 * iqr * n^(-1 / 3)
    rule <- "Freedman-Diaconis"
  } else if (hi > lo) {
    width <- (hi - lo) / ceiling(log2(n) + 1)
    rule <- "Sturges"
  } else {
    width <- 1
    rule <- "unit"
  }
  capped <- hi > lo && (hi - lo) / width > .dist_max_bins
  if (capped) width <- (hi - lo) / .dist_max_bins
  start <- floor(lo / width) * width
  list(width = width, start = start, n_bins = floor((hi - start) / width) + 1,
       rule = rule, capped = capped)
}

# Left-closed bin counts: bin k holds start + (k - 1) * width <= v < start + k * width.
.dist_counts <- function(v, bins) {
  idx <- floor((v - bins$start) / bins$width) + 1
  idx <- pmin(pmax(idx, 1), bins$n_bins)
  tabulate(idx, nbins = bins$n_bins)
}

# Shapiro-Wilk W and P, defined for 3 to 5000 values that are not all equal.
.dist_shapiro <- function(v) {
  n <- length(v)
  if (n < 3 || n > 5000 || stats::sd(v) == 0) return(c(NA_real_, NA_real_))
  st <- tryCatch(stats::shapiro.test(v), error = .dist_null_on_error)
  if (is.null(st)) return(c(NA_real_, NA_real_))
  c(unname(st$statistic), st$p.value)
}

# Count-axis breaks: the extended breaks, whole numbers only, as the Python
# engine's integer tick locator gives them.
.dist_int_breaks <- function(limits) {
  b <- scales::extended_breaks(n = 5)(limits)
  b[abs(b - round(b)) < 1e-9]
}

.dist_bandwidth <- function(v) if (length(v) >= 2) stats::bw.nrd0(v) else NA_real_

# The Gaussian kernel density at one grid point.
.dist_kde_point <- function(g, v, bw) sum(stats::dnorm((g - v) / bw)) / (length(v) * bw)

.dist_p_label <- function(p) {
  if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))
}

# One row of the stats table: the summary of one group.
.dist_group_row <- function(gl, d, bins) {
  v <- d$value[d$group == gl]
  n <- length(v)
  mom <- qc_moment_stats(v)
  sw <- .dist_shapiro(v)
  data.frame(
    group = gl, n = n, mean = mean(v),
    sd = if (n >= 2) stats::sd(v) else NA_real_,
    median = .dist_q7(v, 0.5), q1 = .dist_q7(v, 0.25), q3 = .dist_q7(v, 0.75),
    min = min(v), max = max(v),
    skewness = mom$skewness, excess_kurtosis = mom$kurtosis,
    shapiro_w = sw[1], shapiro_p = sw[2],
    bin_width = bins$width, bandwidth = .dist_bandwidth(v),
    stringsAsFactors = FALSE)
}

# The histogram bars of one group.
.dist_group_hist <- function(gl, d, bins) {
  v <- d$value[d$group == gl]
  k <- seq_len(bins$n_bins)
  data.frame(group = gl, xmin = bins$start + (k - 1) * bins$width,
             xmax = bins$start + k * bins$width, count = .dist_counts(v, bins),
             stringsAsFactors = FALSE)
}

# The density and the normal curve of one group, scaled to counts.
.dist_group_curves <- function(gl, d, bins) {
  v <- d$value[d$group == gl]
  n <- length(v)
  if (n < 2) return(NULL)
  bw <- stats::bw.nrd0(v)
  grid <- seq(min(v) - 3 * bw, max(v) + 3 * bw, length.out = .dist_grid_n)
  to_count <- n * bins$width
  dens <- vapply(grid, .dist_kde_point, numeric(1), v = v, bw = bw)
  out <- data.frame(group = gl, x = grid, y = dens * to_count, curve = "Kernel density",
                    stringsAsFactors = FALSE)
  s <- stats::sd(v)
  if (is.finite(s) && s > 0) {
    out <- rbind(out, data.frame(group = gl, x = grid, y = stats::dnorm(grid, mean(v), s) * to_count,
                                 curve = "Normal fit", stringsAsFactors = FALSE))
  }
  out
}

# The panel title of one group: its name, n and the Shapiro-Wilk P. A long
# title goes on two lines, with the name wrapped at 40 characters.
.dist_strip <- function(i, stats_df, width) {
  r <- stats_df[i, ]
  name <- paste(strwrap(r$group, width = 41), collapse = "\n")
  info <- if (is.na(r$shapiro_p)) sprintf("n = %d", r$n)
          else sprintf("n = %d, Shapiro-Wilk %s", r$n, .dist_p_label(r$shapiro_p))
  if (nchar(name) + nchar(info) + 3 > width) paste0(name, "\n", info) else paste0(name, ": ", info)
}

# The panel grid: one column up to six groups, then two, then three.
.dist_layout <- function(k) {
  ncol <- if (k <= 6) 1 else if (k <= 12) 2 else 3
  nrow <- ceiling(k / ncol)
  panel_h <- if (nrow <= 3) 1.6 else 1.25
  list(ncol = ncol, nrow = nrow,
       width = if (ncol == 1) 5.6 else 1.2 + 3.6 * ncol,
       height = 1.5 + nrow * panel_h,
       strip_chars = if (ncol == 1) 52 else 34)
}

recipe_distribution <- function(df, spec) {
  y <- spec$data$y
  x <- spec$data$x
  has_group <- !is.null(x) && nzchar(x)
  if (is.null(y) || !nzchar(y)) stop("distribution needs --y (the numeric variable)")
  for (cn in c(y, if (has_group) x)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  steps <- character(0)
  n0 <- nrow(df)
  yv <- suppressWarnings(as.numeric(df[[y]]))
  coerced <- sum(is.na(yv) & !is.na(df[[y]]))
  if (coerced > 0) steps <- c(steps, sprintf("Coerced column '%s' to numeric; %d value(s) became NA.", y, coerced))
  grp <- if (has_group) as.character(df[[x]]) else rep("all", n0)
  keep <- is.finite(yv) & !is.na(grp)
  if (sum(!keep) > 0) {
    steps <- c(steps, sprintf("Dropped %d row(s) with a missing value%s (%d -> %d rows).",
                              sum(!keep), if (has_group) " or group" else "", n0, sum(keep)))
  }
  d <- data.frame(group = grp[keep], value = yv[keep], stringsAsFactors = FALSE)
  if (nrow(d) < 2) stop("distribution needs at least 2 finite values")
  if (length(steps) == 0) steps <- "No cleaning applied; input used as-is."
  lv <- if (has_group) sort(unique(d$group), method = "radix") else "all"

  bins <- .dist_bins(d$value)
  stats_df <- do.call(rbind, lapply(lv, .dist_group_row, d = d, bins = bins))
  rownames(stats_df) <- NULL
  hist_df <- do.call(rbind, lapply(lv, .dist_group_hist, d = d, bins = bins))
  curve_df <- do.call(rbind, lapply(lv, .dist_group_curves, d = d, bins = bins))
  lay <- .dist_layout(length(lv))
  strips <- vapply(seq_len(nrow(stats_df)), .dist_strip, character(1), stats_df = stats_df,
                   width = lay$strip_chars)

  p <- .dist_plot(d, hist_df, curve_df, stats_df, strips, lv, lay, bins, has_group, spec, x, y)

  test_meta <- list(
    name = "Shapiro-Wilk normality test per group (descriptive, no comparison)",
    n = nrow(d), n_groups = length(lv),
    bins = list(rule = bins$rule, capped_at_200 = bins$capped, width = bins$width,
                start = bins$start, n_bins = bins$n_bins),
    groups = stats::setNames(lapply(seq_len(nrow(stats_df)), .dist_meta_row, stats_df = stats_df),
                             stats_df$group))
  resolved <- list(family = "distribution", label = test_meta$name)
  qc_panels <- stats::setNames(lapply(lv, .dist_values, d = d), lv)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .dist_methods(stats_df, bins, has_group, spec, x, y),
       build_script = list(fn = .dist_emit_script, args = list(spec = spec, has_group = has_group)),
       label = slugify(y), df_used = d, clean_steps = steps,
       qc = list(quantity = if (has_group) "values of each group" else "values", panels = qc_panels),
       width = if (has_group) lay$width else 5.4, height = if (has_group) lay$height else 4.2)
}

.dist_values <- function(gl, d) d$value[d$group == gl]

.dist_meta_row <- function(i, stats_df) as.list(stats_df[i, setdiff(names(stats_df), "group")])

.dist_plot <- function(d, hist_df, curve_df, stats_df, strips, lv, lay, bins, has_group, spec, x, y) {
  pal <- pub_palette(length(lv), spec$appearance$palette)
  names(pal) <- lv
  hist_df$panel <- factor(strips[match(hist_df$group, lv)], levels = strips)
  if (!is.null(curve_df)) curve_df$panel <- factor(strips[match(curve_df$group, lv)], levels = strips)
  rug_df <- d[d$group %in% stats_df$group[stats_df$n <= .dist_rug_max], , drop = FALSE]
  rug_df$panel <- factor(strips[match(rug_df$group, lv)], levels = strips)

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = hist_df,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = group),
                       colour = "white", linewidth = 0.25, alpha = 0.65)
  if (!is.null(curve_df)) {
    p <- p + ggplot2::geom_line(data = curve_df, ggplot2::aes(x = x, y = y, linetype = curve),
                                colour = "grey15", linewidth = 0.7)
  }
  if (nrow(rug_df) > 0) {
    p <- p + ggplot2::geom_rug(data = rug_df, ggplot2::aes(x = value, colour = group), sides = "b",
                               length = grid::unit(0.04, "npc"), alpha = 0.6, linewidth = 0.3)
  }
  sub <- if (has_group) {
    sprintf("Shapiro-Wilk normality test per group; bin width %s", formatC(bins$width, format = "g", digits = 3))
  } else {
    r <- stats_df[1, ]
    paste0(sprintf("n = %d", r$n),
           if (!is.na(r$shapiro_w)) sprintf(", Shapiro-Wilk W = %.3f, %s", r$shapiro_w, .dist_p_label(r$shapiro_p)) else "",
           if (!is.na(r$skewness)) sprintf("\nskewness = %.2f, excess kurtosis = %.2f", r$skewness, r$excess_kurtosis) else "")
  }
  p <- p +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_linetype_manual(values = c("Kernel density" = "solid", "Normal fit" = "22"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = scales::breaks_extended(n = 6)) +
    ggplot2::scale_y_continuous(breaks = .dist_int_breaks,
                                expand = ggplot2::expansion(mult = c(0.07, 0.06))) +
    ggplot2::labs(x = spec$appearance$x_label %||% spec$appearance$y_label %||% y,
                  y = "Count", title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = grid::unit(2, "lines"),
                   legend.text = ggplot2::element_text(size = 10),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10, hjust = 0, lineheight = 1.05,
                                                      margin = ggplot2::margin(4, 0, 3, 0)),
                   panel.spacing.y = grid::unit(0.6, "lines"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
  if (has_group) {
    # Every panel keeps its own x axis line and ticks, as in the Python engine;
    # only the bottom panels carry the tick labels.
    p <- p + ggplot2::facet_wrap(ggplot2::vars(panel), ncol = lay$ncol, scales = "free_y",
                                 axes = "all_x", axis.labels = "margins")
  }
  p
}

.dist_methods <- function(stats_df, bins, has_group, spec, x, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  scope <- if (has_group) {
    sprintf("in each of %d groups of %s (n = %s)", nrow(stats_df), x,
            paste(sprintf("%d", stats_df$n), collapse = ", "))
  } else sprintf("(n = %d)", stats_df$n[1])
  bin_txt <- switch(bins$rule,
    "Freedman-Diaconis" = "the Freedman-Diaconis rule (2 x IQR x n^(-1/3)) on the pooled values",
    "Sturges" = "the Sturges rule, because the interquartile range of the pooled values is zero",
    "a unit width, because every value is equal")
  cap_txt <- if (isTRUE(bins$capped)) " The rule asked for more than 200 bins, so the width was widened to give 200 bins." else ""
  sw <- stats_df[!is.na(stats_df$shapiro_p), , drop = FALSE]
  sw_txt <- if (nrow(sw) > 0) {
    sprintf(" Normality was assessed with the Shapiro-Wilk test (%s).",
            paste(sprintf("%sW = %.3f, %s", if (has_group) paste0(sw$group, ": ") else "",
                          sw$shapiro_w, vapply(sw$shapiro_p, .dist_p_label, character(1))),
                  collapse = "; "))
  } else " The Shapiro-Wilk test was not run, because it needs 3 to 5000 values that are not all equal."
  cmp_txt <- if (has_group) {
    " No test compared the groups; use two_group_compare or multi_group_compare to test a difference in location."
  } else ""
  paste0(
    sprintf("The distribution of %s was described %s.", y_lab, scope),
    sprintf(" The histogram uses a bin width of %s from %s, with bins anchored at a multiple of the width.",
            formatC(bins$width, format = "g", digits = 4), bin_txt),
    cap_txt,
    " The solid curve is a Gaussian kernel density with the bandwidth of Silverman's rule of thumb (stats::bw.nrd0),",
    " and the dashed curve is a normal distribution with the sample mean and SD; both are scaled to counts.",
    sw_txt,
    " Skewness and excess kurtosis are the moment estimates on standardised values (sum(z^3) / n and sum(z^4) / n - 3, with z scaled by the SD).",
    cmp_txt,
    sprintf(" Computed in R %s and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.dist_emit_script <- function(spec, has_group, in_name, fig_stub) {
  y <- spec$data$y
  x <- spec$data$x
  grp_line <- if (has_group) sprintf('d <- data.frame(group = as.character(raw[["%s"]]), value = as.numeric(raw[["%s"]]))', x, y)
              else sprintf('d <- data.frame(group = "all", value = as.numeric(raw[["%s"]]))', y)
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    grp_line,
    "d <- d[is.finite(d$value) & !is.na(d$group), ]",
    "",
    "# Freedman-Diaconis bins on the pooled values (Sturges when the IQR is zero).",
    "v <- d$value; n <- length(v)",
    "iqr <- diff(quantile(v, c(0.25, 0.75), type = 7, names = FALSE))",
    "w <- if (iqr > 0) 2 * iqr * n^(-1/3) else if (max(v) > min(v)) (max(v) - min(v)) / ceiling(log2(n) + 1) else 1",
    "if (max(v) > min(v) && (max(v) - min(v)) / w > 200) w <- (max(v) - min(v)) / 200",
    "start <- floor(min(v) / w) * w; nb <- floor((max(v) - start) / w) + 1",
    "",
    "bars <- NULL; curves <- NULL",
    "for (g in sort(unique(d$group), method = \"radix\")) {",
    "  vg <- d$value[d$group == g]",
    "  idx <- pmin(pmax(floor((vg - start) / w) + 1, 1), nb)",
    "  bars <- rbind(bars, data.frame(group = g, xmin = start + (seq_len(nb) - 1) * w,",
    "                                 xmax = start + seq_len(nb) * w, count = tabulate(idx, nb)))",
    "  if (length(vg) < 2) next",
    "  bw <- bw.nrd0(vg)",
    "  grid <- seq(min(vg) - 3 * bw, max(vg) + 3 * bw, length.out = 512)",
    "  dens <- rowSums(dnorm(outer(grid, vg, \"-\") / bw)) / (length(vg) * bw)",
    "  curves <- rbind(curves,",
    "    data.frame(group = g, x = grid, y = dens * length(vg) * w, curve = \"Kernel density\"),",
    "    data.frame(group = g, x = grid, y = dnorm(grid, mean(vg), sd(vg)) * length(vg) * w, curve = \"Normal fit\"))",
    "  if (length(vg) >= 3 && length(vg) <= 5000) print(shapiro.test(vg))",
    "}",
    "",
    "p <- ggplot() +",
    "  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = group),",
    "            colour = \"white\", linewidth = 0.25, alpha = 0.65) +",
    "  geom_line(data = curves, aes(x = x, y = y, linetype = curve), colour = \"grey15\", linewidth = 0.7) +",
    "  scale_linetype_manual(values = c(\"Kernel density\" = \"solid\", \"Normal fit\" = \"22\"), name = NULL) +",
    sprintf('  labs(x = "%s", y = "Count") +', spec$appearance$x_label %||% spec$appearance$y_label %||% y),
    "  theme_classic(base_size = 13) +",
    "  theme(legend.position = \"bottom\", strip.background = element_blank()) +",
    "  guides(fill = \"none\")",
    if (has_group) "p <- p + facet_wrap(~ group, ncol = 1, scales = \"free_y\")" else character(0),
    "",
    sprintf('ggsave("%s.pdf", p, width = 5.6, height = %s)', fig_stub,
            if (has_group) "1.5 + 1.6 * length(unique(d$group))" else "4.2"),
    'cat("Reproduced figure.\\n")'
  )
}
