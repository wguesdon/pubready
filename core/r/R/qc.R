# Normality QC. A QQ panel for each quantity a recipe's normality check runs on,
# with the Shapiro-Wilk p, skewness, and excess kurtosis printed on the panel.
# Two consumers use it: the bundle writer (qc_normality_<ts>.png) and
# `figkit diagnose`, so the test choice can be shown to the scientist, not just
# asserted. The QQ plot is the tiebreaker the Shapiro p-value cannot be on its
# own: Shapiro is underpowered at small n and over-rejects at large n.

# Sample skewness and excess kurtosis on the finite values of a vector, plus the
# Shapiro-Wilk p from the shared guarded helper.
qc_moment_stats <- function(v) {
  v <- v[is.finite(v)]
  n <- length(v)
  out <- list(n = n, shapiro_p = NA_real_, skewness = NA_real_, kurtosis = NA_real_)
  if (n < 3) return(out)
  out$shapiro_p <- shapiro_safe(v)$p
  s <- stats::sd(v)
  if (is.finite(s) && s > 0) {
    z <- (v - mean(v)) / s
    out$skewness <- sum(z^3) / n
    out$kurtosis <- sum(z^4) / n - 3
  }
  out
}

# qc = list(quantity = <label>, panels = named list of numeric vectors).
# Returns list(plot = ggplot, panels = data.frame of per-panel stats, quantity).
qc_normality <- function(qc) {
  nm <- names(qc$panels)
  rows <- lapply(nm, function(k) {
    s <- qc_moment_stats(qc$panels[[k]])
    data.frame(panel = k, n = s$n, shapiro_p = s$shapiro_p,
               skewness = s$skewness, kurtosis = s$kurtosis,
               stringsAsFactors = FALSE)
  })
  panels_df <- do.call(rbind, rows)

  fmt <- function(x, g = FALSE) {
    if (is.na(x)) "NA" else formatC(x, format = if (g) "g" else "f", digits = 2)
  }
  strip <- function(k) {
    r <- panels_df[panels_df$panel == k, ]
    sprintf("%s  (n = %d)\nShapiro p = %s | skew = %s | kurtosis = %s",
            k, r$n, fmt(r$shapiro_p, TRUE), fmt(r$skewness), fmt(r$kurtosis))
  }
  long <- do.call(rbind, lapply(nm, function(k) {
    v <- qc$panels[[k]][is.finite(qc$panels[[k]])]
    if (length(v) < 1) return(NULL)
    data.frame(panel = strip(k), value = v, stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot(long, ggplot2::aes(sample = value)) +
    ggplot2::stat_qq(colour = pubplot_palette[1], size = 1.2, alpha = 0.75) +
    ggplot2::stat_qq_line(colour = "#333333", linewidth = 0.5) +
    ggplot2::facet_wrap(~ panel, scales = "free") +
    ggplot2::labs(
      x = "Theoretical quantiles", y = "Sample quantiles",
      title = sprintf("Normality QC: %s", qc$quantity),
      subtitle = "Points near the line support normality; mild, judged deviation is usually acceptable."
    ) +
    theme_pubplot(12) +
    ggplot2::theme(
      strip.text    = ggplot2::element_text(size = 8.5, lineheight = 1.05),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30")
    )
  list(plot = p, panels = panels_df, quantity = qc$quantity)
}

# Size a QC figure to the number of panels (up to three across).
qc_dims <- function(n_panels) {
  ncol <- min(max(n_panels, 1), 3)
  nrow <- ceiling(n_panels / ncol)
  list(width = 2.9 * ncol + 0.4, height = 2.9 * nrow + 1.0)
}
