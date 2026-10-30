# Recipe: meta-analysis forest plot. Input is one row per study with an effect
# estimate (--y), its 95% CI (--ci_low, --ci_high) or its standard error on the
# analysis scale (--se), and a study label (--label). A ratio measure (odds,
# risk or hazard ratio) is pooled on the log scale and a difference measure on
# its own scale. --measure ratio|difference sets the scale; auto picks ratio when
# every value is positive and the CIs are more symmetric on the log scale.
#
# The studies are pooled by inverse variance with metafor: a fixed-effect model
# and a DerSimonian-Laird random-effects model, with Cochran's Q, I^2 and tau^2
# for the heterogeneity. The Python engine uses statsmodels combine_effects,
# which gives the same numbers.

.fm_z <- function() stats::qnorm(0.975)

# Median relative asymmetry of the CIs around the estimates.
.fm_asymmetry <- function(est, lo, hi) {
  stats::median(abs((hi - est) - (est - lo)) / (hi - lo))
}

.fm_resolve_measure <- function(measure, est, lo, hi) {
  measure <- tolower(measure %||% "auto")
  if (measure %in% c("ratio", "difference")) return(measure)
  if (measure != "auto") stop(sprintf("unknown --measure '%s' (use ratio or difference)", measure))
  if (is.null(lo) || any(c(est, lo, hi) <= 0)) return("difference")
  if (.fm_asymmetry(log(est), log(lo), log(hi)) < .fm_asymmetry(est, lo, hi)) "ratio" else "difference"
}

.fm_back <- function(v, measure) if (measure == "ratio") exp(v) else v

# The study labels: the --label column, with "Study <row>" for an empty cell, so
# that no row reads "NA".
.fm_labels <- function(v, n) {
  out <- if (is.null(v)) rep(NA_character_, n) else trimws(as.character(v))
  empty <- is.na(out) | !nzchar(out)
  out[empty] <- sprintf("Study %d", which(empty))
  out
}

# Why one study row cannot be pooled, or NA when it can (F77): its CI bounds are
# reversed, its CI has zero width, its CI does not hold its estimate, or its
# standard error is 0 or below. The Python engine has the same rule.
.fm_study_problem <- function(i, est, lo, hi, se) {
  if (!is.null(se)) {
    if (se[i] <= 0) return(sprintf("its standard error is %s", format(se[i])))
    return(NA_character_)
  }
  if (lo[i] > hi[i]) return(sprintf("its CI bounds are reversed (%s to %s)", format(lo[i]), format(hi[i])))
  if (lo[i] == hi[i]) return(sprintf("its CI has zero width (%s to %s)", format(lo[i]), format(hi[i])))
  if (est[i] < lo[i] || est[i] > hi[i]) {
    return(sprintf("its CI (%s to %s) does not hold its estimate %s", format(lo[i]), format(hi[i]),
                   format(est[i])))
  }
  NA_character_
}

# Stop and name each study that cannot be pooled.
.fm_check_studies <- function(labels, est, lo, hi, se) {
  why <- vapply(seq_along(est), .fm_study_problem, character(1), est = est, lo = lo, hi = hi,
                se = se)
  bad <- which(!is.na(why))
  if (length(bad) == 0) return(invisible(TRUE))
  msg <- paste(sprintf("study '%s': %s", labels[bad], why[bad])[seq_len(min(5, length(bad)))],
               collapse = "; ")
  if (length(bad) > 5) msg <- sprintf("%s; and %d more", msg, length(bad) - 5)
  stop(sprintf("forest_meta cannot pool %d study row(s): %s", length(bad), msg))
}

recipe_forest_meta <- function(df, spec) {
  y <- spec$data$y; lab <- spec$data$label
  lo_col <- spec$data$ci_low; hi_col <- spec$data$ci_high; se_col <- spec$data$se
  if (is.null(y)) stop("forest_meta needs --y (the effect estimate of each study)")
  has_ci <- !is.null(lo_col) && !is.null(hi_col)
  if (!has_ci && is.null(se_col)) {
    stop("forest_meta needs --ci_low and --ci_high (the 95% CI), or --se (the standard error)")
  }
  cols <- c(y, if (has_ci) c(lo_col, hi_col) else se_col, if (!is.null(lab)) lab)
  for (cn in cols) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  n0 <- nrow(df)
  est <- suppressWarnings(as.numeric(df[[y]]))
  lo <- if (has_ci) suppressWarnings(as.numeric(df[[lo_col]])) else NULL
  hi <- if (has_ci) suppressWarnings(as.numeric(df[[hi_col]])) else NULL
  se_in <- if (!has_ci) suppressWarnings(as.numeric(df[[se_col]])) else NULL
  labels <- .fm_labels(if (!is.null(lab)) df[[lab]] else NULL, n0)
  keep <- is.finite(est) & (if (has_ci) is.finite(lo) & is.finite(hi) else is.finite(se_in))
  est <- est[keep]; labels <- labels[keep]
  if (has_ci) { lo <- lo[keep]; hi <- hi[keep] } else se_in <- se_in[keep]
  .fm_check_studies(labels, est, lo, hi, se_in)
  if (length(est) < 2) stop("forest_meta needs at least 2 studies")

  measure <- .fm_resolve_measure(spec$data$measure, est, lo, hi)
  if (measure == "ratio" && any(c(est, lo, hi) <= 0)) {
    stop("a ratio measure needs every estimate and CI bound above zero")
  }
  tr <- if (measure == "ratio") log else identity
  yi <- tr(est)
  sei <- if (has_ci) (tr(hi) - tr(lo)) / (2 * .fm_z()) else se_in
  if (!has_ci) { lo <- .fm_back(yi - .fm_z() * sei, measure); hi <- .fm_back(yi + .fm_z() * sei, measure) }
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d study row(s) with a missing estimate or interval.", sum(!keep)),
    sprintf("Pooled %d of %d study row(s) as a %s measure%s.", length(est), n0, measure,
            if (measure == "ratio") " on the log scale" else ""),
    if (has_ci) "Derived each standard error from the 95% CI as (upper - lower) / (2 x 1.96) on the analysis scale." else NULL)

  fe <- metafor::rma(yi = yi, sei = sei, method = "FE")
  re <- metafor::rma(yi = yi, sei = sei, method = "DL")
  w_fe <- 1 / sei^2
  w_re <- 1 / (sei^2 + re$tau2)
  pooled <- data.frame(
    row_type = "pooled", label = c("Fixed effect", "Random effects (DerSimonian-Laird)"),
    estimate = .fm_back(c(fe$b[1], re$b[1]), measure),
    ci_low = .fm_back(c(fe$ci.lb, re$ci.lb), measure),
    ci_high = .fm_back(c(fe$ci.ub, re$ci.ub), measure),
    yi = c(fe$b[1], re$b[1]), se = c(fe$se, re$se),
    weight_fixed_pct = c(100, NA), weight_random_pct = c(NA, 100),
    z = c(fe$zval, re$zval), p_value = c(fe$pval, re$pval),
    q = re$QE, q_df = re$k - 1, q_p = re$QEp, i2 = re$I2, tau2 = re$tau2,
    stringsAsFactors = FALSE)
  studies <- data.frame(
    row_type = "study", label = labels, estimate = est, ci_low = lo, ci_high = hi,
    yi = yi, se = sei, weight_fixed_pct = 100 * w_fe / sum(w_fe),
    weight_random_pct = 100 * w_re / sum(w_re), z = NA_real_, p_value = NA_real_,
    q = NA_real_, q_df = NA_real_, q_p = NA_real_, i2 = NA_real_, tau2 = NA_real_,
    stringsAsFactors = FALSE)
  stats_df <- rbind(studies, pooled)

  p <- .fm_plot(studies, pooled, measure, spec)
  test_meta <- list(
    name = "Inverse-variance meta-analysis (fixed effect and DerSimonian-Laird random effects)",
    measure = measure, k = length(est),
    fixed = list(estimate = pooled$estimate[1], ci_low = pooled$ci_low[1], ci_high = pooled$ci_high[1],
                 p_value = pooled$p_value[1]),
    random = list(estimate = pooled$estimate[2], ci_low = pooled$ci_low[2], ci_high = pooled$ci_high[2],
                  p_value = pooled$p_value[2]),
    heterogeneity = list(q = re$QE, df = re$k - 1, p_value = re$QEp, i2 = re$I2, tau2 = re$tau2))
  resolved <- list(family = "meta", label = test_meta$name)
  list(plot = p$plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .fm_methods(pooled, measure, length(est), y, spec),
       build_script = list(fn = .fm_emit_script,
                           args = list(spec = spec, measure = measure, has_ci = has_ci)),
       label = slugify(y),
       df_used = data.frame(label = labels, estimate = est, ci_low = lo, ci_high = hi,
                            yi = yi, se = sei, stringsAsFactors = FALSE),
       clean_steps = clean_steps, width = p$width, height = p$height)
}

.fm_diamond <- function(i, pooled, ypos) {
  data.frame(id = i, x = c(pooled$ci_low[i], pooled$estimate[i], pooled$ci_high[i], pooled$estimate[i]),
             y = ypos[i] + c(0, 0.32, 0, -0.32))
}

.fm_plot <- function(studies, pooled, measure, spec) {
  k <- nrow(studies)
  sy <- rev(seq_len(k)) + 2
  py <- c(1.5, 0.5)
  pal <- spec$appearance$palette %||% pubready_palette
  sd <- data.frame(y = sy, est = studies$estimate, lo = studies$ci_low, hi = studies$ci_high,
                   w = studies$weight_random_pct)
  dia <- do.call(rbind, lapply(1:2, .fm_diamond, pooled = pooled, ypos = py))
  null_x <- if (measure == "ratio") 1 else 0
  txt <- data.frame(
    y = c(sy, py),
    est = c(sprintf("%.2f [%.2f, %.2f]", studies$estimate, studies$ci_low, studies$ci_high),
            sprintf("%.2f [%.2f, %.2f]", pooled$estimate, pooled$ci_low, pooled$ci_high)),
    w = c(sprintf("%.1f%%", studies$weight_random_pct), "", ""))
  het <- sprintf("Heterogeneity: Q(%d) = %.1f, P = %s; I2 = %.1f%%; tau2 = %.3g",
                 as.integer(pooled$q_df[2]), pooled$q[2], sprintf("%.2g", pooled$q_p[2]),
                 pooled$i2[2], pooled$tau2[2])
  ylim <- c(0, max(sy) + 0.8)
  brk <- c(sy, py)
  # A long study label wraps at 30 characters, and the figure grows with the
  # longest line, so the panel keeps its width (F81). Python has the same rule.
  blab <- c(.est_wrap(studies$label, 30), "Fixed effect", "Random effects (DL)")
  size <- .fm_size(blab, k)

  forest <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = null_x, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = pooled$estimate[2], linetype = "dotted", colour = pal[2],
                        linewidth = 0.4) +
    ggplot2::geom_errorbarh(data = sd, ggplot2::aes(y = y, xmin = lo, xmax = hi), height = 0,
                            linewidth = 0.5, colour = "grey20") +
    ggplot2::geom_point(data = sd, ggplot2::aes(x = est, y = y, size = w), shape = 15,
                        colour = pal[1]) +
    ggplot2::geom_polygon(data = dia, ggplot2::aes(x = x, y = y, group = id),
                          fill = pal[2], colour = "grey20", linewidth = 0.3) +
    ggplot2::scale_size_area(max_size = 5, guide = "none") +
    ggplot2::scale_y_continuous(breaks = brk, labels = blab, limits = ylim,
                                expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::labs(x = spec$appearance$x_label %||%
                    sprintf("%s (95%% CI%s)", spec$data$y, if (measure == "ratio") ", log scale" else ""),
                  y = NULL, title = spec$appearance$title, subtitle = het) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(hjust = 0, size = 10),
                   plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey30"))
  if (measure == "ratio") {
    # The same 1-2-5 ticks as the Python engine, with at least one tick on each
    # side of the null line, and the axis widened to show every tick.
    lo <- min(c(studies$ci_low, pooled$ci_low)); hi <- max(c(studies$ci_high, pooled$ci_high))
    ticks <- log_ticks_ref(lo, hi, 1)
    forest <- forest +
      ggplot2::scale_x_log10(breaks = ticks, limits = c(min(lo, ticks), max(hi, ticks)),
                             labels = scales::label_number(drop0trailing = TRUE, big.mark = ""))
  }

  table <- ggplot2::ggplot(txt, ggplot2::aes(y = y)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = est), hjust = 0, size = 3.3) +
    ggplot2::geom_text(ggplot2::aes(x = 1, label = w), hjust = 1, size = 3.3) +
    ggplot2::annotate("text", x = c(0, 1), y = max(sy) + 0.8,
                      label = c("Estimate [95% CI]", "Weight"), hjust = c(0, 1),
                      fontface = "bold", size = 3.3) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(add = 0.2)) +
    ggplot2::theme_void() +
    ggplot2::labs(subtitle = " ", x = " ") +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 9.5),
                   axis.title.x = ggplot2::element_text(size = 13, colour = "white"),
                   axis.text.x = ggplot2::element_text(size = 11, colour = "white"),
                   plot.margin = ggplot2::margin(10, 10, 10, 0))
  list(plot = ggpubr::ggarrange(forest, table, ncol = 2,
                                widths = c(size$width - 8.2 / 3.6, 8.2 / 3.6), align = "h"),
       height = size$height, width = size$width)
}

.fm_methods <- function(pooled, measure, k, y, spec) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions(c("metafor", "ggplot2"))
  scale_txt <- if (measure == "ratio") "on the log scale and back-transformed" else "on the original scale"
  paste0(
    sprintf("%d studies were pooled by inverse-variance weighting %s, with the standard error of each study derived from its 95%% confidence interval or given directly. ",
            k, scale_txt),
    sprintf("The fixed-effect estimate of %s was %.3f (95%% CI %.3f to %.3f, P = %s), ",
            y, pooled$estimate[1], pooled$ci_low[1], pooled$ci_high[1], sprintf("%.2g", pooled$p_value[1])),
    sprintf("and the DerSimonian-Laird random-effects estimate was %.3f (95%% CI %.3f to %.3f, P = %s). ",
            pooled$estimate[2], pooled$ci_low[2], pooled$ci_high[2], sprintf("%.2g", pooled$p_value[2])),
    sprintf("Heterogeneity was assessed with Cochran's Q (Q = %.2f, df = %d, P = %s), I2 = %.1f%% and tau2 = %.4f%s. ",
            pooled$q[2], as.integer(pooled$q_df[2]), sprintf("%.2g", pooled$q_p[2]), pooled$i2[2], pooled$tau2[2],
            .fm_truncation_text(pooled$q[2], pooled$q_df[2])),
    "Squares are sized by the random-effects weight and the diamonds show the pooled estimates. ",
    sprintf("Computed in R %s with metafor %s and drawn with ggplot2 %s.", rver, v$metafor, v$ggplot2))
}

.fm_emit_script <- function(spec, measure, has_ci, in_name, fig_stub) {
  y <- spec$data$y; lab <- spec$data$label %||% ""
  tr <- if (measure == "ratio") "log" else "identity"
  se_line <- if (has_ci) {
    sprintf('sei <- (%s(raw[["%s"]]) - %s(raw[["%s"]])) / (2 * qnorm(0.975))',
            tr, spec$data$ci_high, tr, spec$data$ci_low)
  } else sprintf('sei <- raw[["%s"]]', spec$data$se)
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(metafor))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('yi <- %s(raw[["%s"]])', tr, y),
    se_line,
    "fe <- rma(yi = yi, sei = sei, method = \"FE\"); print(fe)",
    "re <- rma(yi = yi, sei = sei, method = \"DL\"); print(re)",
    sprintf('pdf("%s.pdf", width = 8.2, height = 1.6 + 0.3 * (length(yi) + 3))', fig_stub),
    sprintf('forest(re, slab = %s, atransf = %s)',
            if (nzchar(lab)) sprintf('raw[["%s"]]', lab) else "NULL",
            if (measure == "ratio") "exp" else "NULL"),
    "invisible(dev.off())",
    'cat("Reproduced figure.\\n")'
  )
}

# The DerSimonian-Laird tau2 is truncated at zero only when Q is below its
# degrees of freedom, so the methods paragraph names the truncation only then.
.fm_truncation_text <- function(q, q_df) {
  if (is.finite(q) && is.finite(q_df) && q < q_df) {
    "; Q was below its degrees of freedom, so the DerSimonian-Laird tau2 was truncated at zero and I2 is 0%"
  } else ""
}

# The figure size: the label column grows with the longest line of the wrapped
# labels, at least the 2.21 inches of the default layout, and each extra line of
# a wrapped label adds 0.2 inch of height (F81). The Python engine has the same
# rule.
.fm_size <- function(labels, k) {
  parts <- strsplit(labels, "\n", fixed = TRUE)
  label_in <- max(0.27 * 8.2, 0.25 + 0.075 * max(nchar(unlist(parts))))
  list(width = 8.2 - 0.27 * 8.2 + label_in,
       height = 1.6 + 0.3 * (k + 3) + 0.2 * sum(lengths(parts) - 1))
}
