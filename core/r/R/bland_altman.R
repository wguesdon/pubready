# Recipe: Bland-Altman agreement between two methods that measure the same
# quantity on the same subjects. Input is one row per subject with the reading of
# method A in --x and of method B in --y. Each difference is A - B and each mean
# is (A + B) / 2. The figure plots the difference against the mean with the bias
# (the mean difference) and the 95% limits of agreement (bias +/- 1.96 SD), each
# with a shaded 95% confidence interval.
#
# The bias CI is t-based. The limits-of-agreement CI uses the standard error of
# Bland and Altman (1999), SD * sqrt(1/n + 1.96^2 / (2 (n - 1))). A regression
# of the difference on the mean tests for a proportional bias, and the normality
# of the differences, which the limits assume, goes to the QC panel.

.ba_stats <- function(a, b) {
  d <- a - b; m <- (a + b) / 2; n <- length(d)
  bias <- mean(d); s <- stats::sd(d)
  z <- stats::qnorm(0.975); tq <- stats::qt(0.975, n - 1)
  se_bias <- s / sqrt(n)
  se_loa <- s * sqrt(1 / n + z^2 / (2 * (n - 1)))
  fit <- summary(stats::lm(d ~ m))$coefficients
  list(n = n, bias = bias, bias_ci_low = bias - tq * se_bias, bias_ci_high = bias + tq * se_bias,
       sd_diff = s, loa_low = bias - z * s, loa_high = bias + z * s,
       loa_low_ci_low = bias - z * s - tq * se_loa, loa_low_ci_high = bias - z * s + tq * se_loa,
       loa_high_ci_low = bias + z * s - tq * se_loa, loa_high_ci_high = bias + z * s + tq * se_loa,
       prop_bias_slope = fit["m", "Estimate"], prop_bias_p = fit["m", "Pr(>|t|)"],
       shapiro_p_diff = shapiro_safe(d)$p)
}

# One label in the right margin, just outside the panel, level with its line.
# The plot turns clipping off so that the label can sit outside the panel.
.ba_margin_label <- function(i, lines) {
  grob <- grid::textGrob(lines$lab[i], x = grid::unit(1, "npc") + grid::unit(5, "pt"),
                         just = "left", gp = grid::gpar(fontsize = 9, col = "grey20"))
  ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf,
                             ymin = lines$y[i], ymax = lines$y[i])
}

recipe_bland_altman <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y
  if (is.null(x) || is.null(y)) {
    stop("bland_altman needs --x (method A) and --y (method B), one row per subject")
  }
  for (cn in c(x, y)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  n0 <- nrow(df)
  a <- suppressWarnings(as.numeric(df[[x]]))
  b <- suppressWarnings(as.numeric(df[[y]]))
  keep <- is.finite(a) & is.finite(b)
  a <- a[keep]; b <- b[keep]
  if (length(a) < 3) stop("bland_altman needs at least 3 complete pairs")
  clean_steps <- if (sum(!keep) > 0) {
    sprintf("Dropped %d row(s) with a missing or non-numeric reading in '%s' or '%s' (%d -> %d).",
            sum(!keep), x, y, n0, length(a))
  } else "No cleaning applied; input used as-is."

  s <- .ba_stats(a, b)
  pd <- data.frame(mean = (a + b) / 2, diff = a - b)
  pal <- spec$appearance$palette %||% pubready_palette
  bands <- data.frame(
    ymin = c(s$loa_low_ci_low, s$bias_ci_low, s$loa_high_ci_low),
    ymax = c(s$loa_low_ci_high, s$bias_ci_high, s$loa_high_ci_high))
  lines <- data.frame(y = c(s$loa_low, s$bias, s$loa_high),
                      lty = c("dashed", "solid", "dashed"),
                      lab = c(paste("-1.96 SD:", sig_text(s$loa_low)), paste("Bias:", sig_text(s$bias)),
                              paste("+1.96 SD:", sig_text(s$loa_high))))

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = mean, y = diff)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = bands$ymin, ymax = bands$ymax,
                      fill = "grey80", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = lines$y, linetype = lines$lty, colour = "grey25",
                        linewidth = 0.5) +
    ggplot2::geom_point(colour = pal[1], size = 2.2, alpha = 0.85) +
    # The bias and limit labels go in the right margin, next to their lines, so
    # that they never cover a point.
    lapply(seq_len(nrow(lines)), .ba_margin_label, lines = lines) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 6)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = spec$appearance$x_label %||% sprintf("Mean of %s and %s", x, y),
                  y = spec$appearance$y_label %||% sprintf("Difference (%s - %s)", x, y),
                  title = spec$appearance$title,
                  subtitle = sprintf("n = %d; proportional bias slope %.3f, P = %s",
                                     s$n, s$prop_bias_slope, sig_text(s$prop_bias_p, 2))) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   plot.margin = ggplot2::margin(10, 12 + 5.2 * max(nchar(lines$lab)), 10, 10))

  stats_df <- data.frame(method_a = x, method_b = y, s, stringsAsFactors = FALSE)
  test_meta <- c(list(name = "Bland-Altman agreement analysis"), s)
  resolved <- list(family = "bland_altman", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .ba_methods(s, x, y, spec),
       build_script = list(fn = .ba_emit_script, args = list(spec = spec)),
       label = slugify(sprintf("%s_vs_%s", x, y)),
       df_used = data.frame(a = a, b = b, mean = pd$mean, difference = pd$diff),
       clean_steps = clean_steps,
       qc = list(quantity = "differences between the methods",
                 panels = stats::setNames(list(pd$diff), sprintf("%s - %s", x, y))),
       width = 6.4, height = 4.4)
}

.ba_methods <- function(s, x, y, spec) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  paste0(
    sprintf("Agreement between %s and %s was assessed by the method of Bland and Altman (n = %d pairs). ", x, y, s$n),
    sprintf("The mean difference (%s - %s) was %s (95%% CI %s to %s) with a standard deviation of %s. ",
            x, y, sig_text(s$bias), sig_text(s$bias_ci_low), sig_text(s$bias_ci_high), sig_text(s$sd_diff)),
    sprintf("The 95%% limits of agreement were %s (95%% CI %s to %s) and %s (95%% CI %s to %s), ",
            sig_text(s$loa_low), sig_text(s$loa_low_ci_low), sig_text(s$loa_low_ci_high),
            sig_text(s$loa_high), sig_text(s$loa_high_ci_low), sig_text(s$loa_high_ci_high)),
    "computed as the mean difference plus or minus 1.96 standard deviations, with confidence intervals from the standard error of Bland and Altman (1999). ",
    sprintf("A linear regression of the difference on the mean gave a slope of %.3f (P = %s) as a test of proportional bias. ",
            s$prop_bias_slope, sig_text(s$prop_bias_p, 2)),
    sprintf("The normality of the differences, which the limits assume, was checked with the Shapiro-Wilk test (P = %s). ",
            sig_text(s$shapiro_p_diff, 2)),
    sprintf("Computed in R %s and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.ba_emit_script <- function(spec, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('a <- as.numeric(raw[["%s"]]); b <- as.numeric(raw[["%s"]])', x, y),
    "keep <- is.finite(a) & is.finite(b); a <- a[keep]; b <- b[keep]",
    "d <- data.frame(mean = (a + b) / 2, diff = a - b)",
    "bias <- mean(d$diff); s <- sd(d$diff)",
    "loa <- bias + c(-1, 1) * qnorm(0.975) * s",
    'cat(sprintf("bias %.4f, SD %.4f, limits %.4f to %.4f\\n", bias, s, loa[1], loa[2]))',
    "print(summary(lm(diff ~ mean, data = d)))",
    "p <- ggplot(d, aes(mean, diff)) + geom_point() +",
    '  geom_hline(yintercept = bias) + geom_hline(yintercept = loa, linetype = "dashed") +',
    sprintf('  labs(x = "Mean of %s and %s", y = "Difference (%s - %s)") + theme_classic()', x, y, x, y),
    sprintf('ggsave("%s.pdf", p, width = 6.4, height = 4.4)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
