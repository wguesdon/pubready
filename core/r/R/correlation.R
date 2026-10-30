# Recipe: correlation between two continuous variables.
# Chooses Pearson or Spearman from the normality of each variable (unless the
# caller forces the method), draws a scatter (with a linear fit and CI band for
# Pearson only, because a rank correlation describes no line), and annotates the coefficient and p on the plot. Returns the plot, a tidy stats
# table, and test metadata.

# Resolve which correlation method to run from the request and the normality of
# each variable. "auto" uses Pearson only when both variables look normal.
resolve_correlation_method <- function(method, x_normal, y_normal) {
  method <- tolower(method %||% "auto")
  if (method == "auto") {
    m <- if (isTRUE(x_normal) && isTRUE(y_normal)) "pearson" else "spearman"
  } else if (method %in% c("pearson", "spearman", "kendall")) {
    m <- method
  } else {
    stop(sprintf("unknown correlation method '%s' (use pearson, spearman, or kendall)", method))
  }
  label <- switch(m,
    pearson  = "Pearson correlation",
    spearman = "Spearman rank correlation",
    kendall  = "Kendall rank correlation")
  estimate_name <- switch(m, pearson = "r", spearman = "rho", kendall = "tau")
  list(method = m, family = "correlation", label = label, estimate_name = estimate_name)
}

# The subtitle: the method, the coefficient symbol and the P value, for example
# "Pearson correlation, r = 0.71, P = 2.1e-10". Python writes the same text.
cor_label <- function(method_label, estimate_name, estimate, p) {
  sprintf("%s, %s = %.2f, P = %s", method_label, estimate_name, estimate, sprintf("%.2g", p))
}

# The --test values that correlation knows (F63).
.correlation_tests <- c("auto", "pearson", "spearman", "kendall")

recipe_correlation <- function(df, spec) {
  check_test_method("correlation", spec$test$method, .correlation_tests)
  x <- spec$data$x
  y <- spec$data$y
  if (is.null(x) || is.null(y)) {
    stop("correlation needs --x and --y (two numeric columns)")
  }
  if (!x %in% names(df)) stop(sprintf("column '%s' not found in data", x))
  if (!y %in% names(df)) stop(sprintf("column '%s' not found in data", y))

  # Clean inline: both columns numeric, drop non-finite pairs.
  n0 <- nrow(df)
  xv <- suppressWarnings(as.numeric(df[[x]]))
  yv <- suppressWarnings(as.numeric(df[[y]]))
  keep <- is.finite(xv) & is.finite(yv)
  xv <- xv[keep]; yv <- yv[keep]
  df_used <- data.frame(a = xv, b = yv)
  names(df_used) <- c(x, y)
  clean_steps <- if (sum(!keep) > 0) {
    sprintf("Dropped %d row(s) with a non-numeric or missing value in '%s' or '%s' (%d -> %d).",
            sum(!keep), x, y, n0, length(xv))
  } else "No cleaning applied; input used as-is."
  if (length(xv) < 3) stop("correlation needs at least 3 complete pairs")
  .cor_check_constant(xv, yv, x, y)

  x_norm <- shapiro_safe(xv)
  y_norm <- shapiro_safe(yv)
  resolved <- resolve_correlation_method(spec$test$method, x_norm$normal, y_norm$normal)

  ct <- stats::cor.test(xv, yv, method = resolved$method)
  estimate  <- unname(ct$estimate)
  ci        <- if (!is.null(ct$conf.int)) as.numeric(ct$conf.int) else c(NA_real_, NA_real_)
  statistic <- if (!is.null(ct$statistic)) unname(ct$statistic) else NA_real_
  dfree     <- if (!is.null(ct$parameter)) unname(ct$parameter) else NA_real_

  pal <- spec$appearance$palette %||% pubready_palette
  point_col <- pal[1]

  # The method, the coefficient and the P value go in a subtitle, in the style of
  # the other recipes and in the same words as the Python engine.
  # A least-squares line describes a Pearson correlation only. A rank correlation
  # measures a monotonic association, so its scatter gets no fitted line.
  p <- ggplot2::ggplot(df_used, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
    ggplot2::geom_point(colour = point_col, size = 2, alpha = 0.8)
  if (resolved$method == "pearson") {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                                  colour = "#333333", fill = "grey75", linewidth = 0.6)
  }
  p <- p +
    ggplot2::labs(
      x        = spec$appearance$x_label %||% x,
      y        = spec$appearance$y_label %||% y,
      title    = spec$appearance$title,
      subtitle = cor_label(resolved$label, resolved$estimate_name, estimate, ct$p.value)) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  ylim <- spec$appearance$y_limits
  ymin <- if (length(ylim) >= 1) ylim[[1]] else NULL
  ymax <- if (length(ylim) >= 2) ylim[[2]] else NULL
  if (!is.null(ymin) || !is.null(ymax)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax))
  }

  stats_df <- data.frame(
    x             = x,
    y             = y,
    n             = length(xv),
    method        = resolved$label,
    estimate_name = resolved$estimate_name,
    estimate      = estimate,
    statistic     = statistic,
    df            = dfree,
    p_value       = unname(ct$p.value),
    conf_low      = ci[1],
    conf_high     = ci[2],
    shapiro_p_x   = x_norm$p,
    shapiro_p_y   = y_norm$p,
    stringsAsFactors = FALSE)

  test_meta <- list(
    name      = resolved$label,
    estimate  = list(name = resolved$estimate_name, value = estimate),
    statistic = statistic,
    df        = dfree,
    p_value   = unname(ct$p.value),
    conf_int  = list(low = ci[1], high = ci[2]),
    assumptions = list(
      normality = list(test = "shapiro-wilk", p_x = x_norm$p, p_y = y_norm$p)))

  methods <- .cor_methods(resolved, stats_df, spec)
  build_script <- list(fn = .cor_emit_script,
                       args = list(spec = spec, resolved = resolved, point_col = point_col))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(sprintf("%s_vs_%s", y, x)),
       df_used = df_used, clean_steps = clean_steps,
       qc = list(quantity = "variable values",
                 panels = stats::setNames(list(xv, yv), c(x, y))),
       width = 4.6, height = 4.2)
}

# A correlation divides by the standard deviation of each variable, so a
# constant column stops the run and is named. The Python engine has the same rule.
.cor_check_constant <- function(xv, yv, x, y) {
  const <- c(x, y)[c(is_zero_variance(xv), is_zero_variance(yv))]
  if (length(const) == 1) {
    stop(sprintf("correlation cannot run because the column %s is constant (all values are equal)",
                 quote_join(const)))
  }
  if (length(const) == 2) {
    stop(sprintf("correlation cannot run because the columns %s are constant (all values are equal)",
                 quote_join(const)))
  }
  invisible(NULL)
}

.cor_methods <- function(resolved, stats_df, spec) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  x_lab <- spec$appearance$x_label %||% stats_df$x[1]
  y_lab <- spec$appearance$y_label %||% stats_df$y[1]
  auto_txt <- if (tolower(spec$test$method %||% "auto") == "auto") {
    " The method was chosen from the Shapiro-Wilk normality of each variable."
  } else ""
  est <- formatC(stats_df$estimate[1], format = "f", digits = 2)
  pv  <- sprintf("%.2g", stats_df$p_value[1])
  paste0(
    sprintf("The association between %s and %s was assessed with a %s (n = %d).",
            x_lab, y_lab, resolved$label, stats_df$n[1]),
    auto_txt,
    sprintf(" The correlation was %s = %s (P = %s)%s",
            resolved$estimate_name, est, pv, .cor_line_text(resolved$method)),
    " Significance was set at P < 0.05.",
    sprintf(" Rendered in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.cor_emit_script <- function(spec, resolved, point_col, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  xlab <- spec$appearance$x_label %||% x
  ylab <- spec$appearance$y_label %||% y
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('x <- suppressWarnings(as.numeric(raw[["%s"]]))', x),
    sprintf('y <- suppressWarnings(as.numeric(raw[["%s"]]))', y),
    "keep <- is.finite(x) & is.finite(y)",
    "df <- data.frame(x = x[keep], y = y[keep])",
    sprintf('ct <- cor.test(df$x, df$y, method = "%s")', resolved$method),
    "print(ct)",
    sprintf('sub <- sprintf("%s, %s = %%.2f, P = %%.2g", ct$estimate, ct$p.value)',
            resolved$label, resolved$estimate_name),
    "",
    "p <- ggplot(df, aes(x = x, y = y)) +",
    sprintf('  geom_point(colour = "%s", size = 2, alpha = 0.8) +', point_col),
    if (resolved$method == "pearson") {
      '  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#333333", fill = "grey75") +'
    },
    sprintf('  labs(x = "%s", y = "%s", subtitle = sub) + theme_classic() +', xlab, ylab),
    '  theme(plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"))',
    sprintf('ggsave("%s.pdf", p, width = 4.6, height = 4.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}

# The end of the result sentence: a Pearson scatter shows the least-squares line
# and its band, and a rank correlation shows the points only.
.cor_line_text <- function(method) {
  if (method == "pearson") {
    ", with a linear fit and 95% confidence band shown."
  } else ". No line is drawn, because a rank correlation measures a monotonic association and not a straight line."
}
