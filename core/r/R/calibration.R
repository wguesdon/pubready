# Recipe: calibration of a risk prediction model. Do the predicted probabilities
# match the observed event rates?
# The outcome is --y and the predicted probabilities are --x, one column per
# model, separated by commas. The event level is --positive, or the second level
# in sorted order, as in roc_curve. For each model the recipe reports the Brier
# score, the calibration intercept (calibration-in-the-large: the intercept of a
# logistic model of the outcome with logit(p) as an offset), the calibration
# slope (the coefficient of logit(p) in a logistic model), each with a 95% Wald
# CI, the Hosmer-Lemeshow test and the ratio of observed to expected events.
#
# The Hosmer-Lemeshow groups are those of ResourceSelection::hoslem.test: the
# breaks are the unique type-7 quantiles of p at 0, 1/g, ..., 1 (g = 10), and
# cut(include.lowest = TRUE) assigns the subjects. Tied quantiles merge, so a
# model with few distinct values has fewer groups, and df = groups - 2. The
# Python engine restates the same rules, quantile for quantile.

.cal_g <- 10L
.cal_clip <- 1e-8

.cal_levels <- function(v) sort(unique(v), method = "radix")

# The Hosmer-Lemeshow group of each subject, 1 to the count of groups.
.cal_groups <- function(p, g) {
  qq <- unique(stats::quantile(p, probs = seq(0, 1, 1 / g), names = FALSE, type = 7))
  if (length(qq) < 2) return(rep(1L, length(p)))
  as.integer(cut(p, breaks = qq, include.lowest = TRUE))
}

# Per-group counts for the test and the figure: size, events, expected events
# and expected non-events, mean prediction, observed proportion and its 95%
# Wilson interval. Empty groups are left out.
.cal_group_table <- function(p, y, g) {
  grp <- .cal_groups(p, g)
  ids <- sort(unique(grp))
  n <- as.vector(tabulate(grp, nbins = max(ids))[ids])
  o1 <- as.vector(tapply(y, grp, sum)[as.character(ids)])
  e1 <- as.vector(tapply(p, grp, sum)[as.character(ids)])
  e0 <- as.vector(tapply(1 - p, grp, sum)[as.character(ids)])
  obs <- o1 / n
  z <- stats::qnorm(0.975)
  centre <- (obs + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(obs * (1 - obs) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  data.frame(group = ids, n = n, events = o1, expected = e1, expected_non = e0,
             mean_pred = e1 / n, observed = obs, ci_low = centre - half, ci_high = centre + half)
}

.cal_hosmer_lemeshow <- function(tab) {
  stat <- sum((tab$events - tab$expected)^2 / tab$expected) +
    sum(((tab$n - tab$events) - tab$expected_non)^2 / tab$expected_non)
  df <- nrow(tab) - 2
  p <- if (df > 0 && is.finite(stat)) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_
  list(statistic = stat, df = df, p_value = p, groups = nrow(tab))
}

# Coefficient, standard error and 95% Wald CI of one term of a logistic model.
.cal_wald <- function(fit, term) {
  co <- summary(fit)$coefficients
  est <- unname(co[term, "Estimate"]); se <- unname(co[term, "Std. Error"])
  z <- stats::qnorm(0.975)
  list(estimate = est, se = se, ci_low = est - z * se, ci_high = est + z * se)
}

.cal_model <- function(m, d, case) {
  p <- d[[m]]
  y <- as.integer(d$outcome == case)
  lp <- stats::qlogis(pmin(pmax(p, .cal_clip), 1 - .cal_clip))
  dd <- data.frame(y = y, lp = lp)
  ctl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
  fit_int <- stats::glm(y ~ 1 + offset(lp), family = stats::binomial(), data = dd, control = ctl)
  fit_slope <- stats::glm(y ~ lp, family = stats::binomial(), data = dd, control = ctl)
  # A second fit that starts at the estimate reads the standard error at the
  # estimate itself, not at the last but one iteration, as the Python engine does.
  fit_int <- stats::glm(y ~ 1 + offset(lp), family = stats::binomial(), data = dd, control = ctl,
                        start = stats::coef(fit_int))
  fit_slope <- stats::glm(y ~ lp, family = stats::binomial(), data = dd, control = ctl,
                          start = stats::coef(fit_slope))
  int <- .cal_wald(fit_int, "(Intercept)")
  slope <- .cal_wald(fit_slope, "lp")
  tab <- .cal_group_table(p, y, .cal_g)
  hl <- .cal_hosmer_lemeshow(tab)
  sm <- stats::lowess(p, y, f = 2 / 3, iter = 0, delta = 0.01 * diff(range(p)))
  list(model = m, n = length(y), n_events = sum(y), brier = mean((p - y)^2),
       intercept = int, slope = slope, hl = hl, oe_ratio = sum(y) / sum(p),
       groups = tab, smooth = data.frame(x = sm$x, y = pmin(pmax(sm$y, 0), 1)),
       p = p, y = y)
}

.cal_row <- function(r, case) {
  data.frame(model = r$model, case_level = case, n = r$n, n_events = r$n_events,
             brier = r$brier,
             intercept = r$intercept$estimate, intercept_ci_low = r$intercept$ci_low,
             intercept_ci_high = r$intercept$ci_high,
             slope = r$slope$estimate, slope_ci_low = r$slope$ci_low,
             slope_ci_high = r$slope$ci_high,
             hl_statistic = r$hl$statistic, hl_df = r$hl$df, hl_p = r$hl$p_value,
             hl_groups = r$hl$groups, oe_ratio = r$oe_ratio, stringsAsFactors = FALSE)
}

.cal_p_text <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) "< 0.001" else sprintf("= %.2g", p)
}

# A legend entry longer than .cal_wrap characters puts the model name on its own
# line, so the legend stays within the figure width.
.cal_wrap <- 62

.cal_legend_label <- function(r) {
  st <- sprintf("slope %.2f, intercept %.2f, Brier %.3f, HL P %s", r$slope$estimate,
                r$intercept$estimate, r$brier, .cal_p_text(r$hl$p_value))
  one <- sprintf("%s: %s", r$model, st)
  if (nchar(one) > .cal_wrap) sprintf("%s\n%s", r$model, st) else one
}

# The layout below the calibration panel: one band per model for the spike
# histogram of the predictions, events up and non-events down.
.cal_band_base <- function(k) -0.17 - 0.24 * (k - 1)
.cal_band_bottom <- function(n_models) .cal_band_base(n_models) - 0.11

.cal_spikes <- function(k, fits, top) {
  r <- fits[[k]]
  bin <- pmin(floor(r$p * 100), 99)
  up <- tabulate(bin[r$y == 1] + 1, nbins = 100)
  down <- tabulate(bin[r$y == 0] + 1, nbins = 100)
  x <- (seq_len(100) - 0.5) / 100
  base <- .cal_band_base(k)
  rows <- rbind(data.frame(x = x, y0 = base, y1 = base + 0.09 * up / top, n = up),
                data.frame(x = x, y0 = base, y1 = base - 0.09 * down / top, n = down))
  rows <- rows[rows$n > 0, , drop = FALSE]
  rows$model <- r$model
  rows
}

.cal_spike_max <- function(r) {
  bin <- pmin(floor(r$p * 100), 99)
  max(tabulate(bin[r$y == 1] + 1, nbins = 100), tabulate(bin[r$y == 0] + 1, nbins = 100))
}

.cal_points <- function(k, fits, n_models) {
  r <- fits[[k]]
  off <- (k - (n_models + 1) / 2) * 0.008
  data.frame(model = r$model, x = r$groups$mean_pred + off, y = r$groups$observed,
             lo = r$groups$ci_low, hi = r$groups$ci_high, stringsAsFactors = FALSE)
}

.cal_smooth <- function(r) data.frame(model = r$model, x = r$smooth$x, y = r$smooth$y,
                                      stringsAsFactors = FALSE)

recipe_calibration <- function(df, spec) {
  y <- spec$data$y
  x_raw <- spec$data$x
  if (is.null(y) || is.null(x_raw)) {
    stop("calibration needs --y (the binary outcome) and --x (one column of predicted probabilities, or several separated by commas)")
  }
  models <- trimws(strsplit(x_raw, ",")[[1]])
  models <- models[nzchar(models)]
  for (cn in c(y, models)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(outcome = as.character(df[[y]]), stringsAsFactors = FALSE)
  for (m in models) d[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  keep <- !is.na(d$outcome) & stats::complete.cases(d[, models, drop = FALSE])
  d <- d[keep, , drop = FALSE]
  for (m in models) {
    if (any(d[[m]] < 0 | d[[m]] > 1)) {
      stop(sprintf("column '%s' holds values outside [0, 1]; calibration needs predicted probabilities", m))
    }
  }
  lv <- .cal_levels(d$outcome)
  if (length(lv) != 2) {
    stop(sprintf("calibration needs exactly 2 outcome levels in '%s'; found %d (%s)",
                 y, length(lv), paste(lv, collapse = ", ")))
  }
  case <- as.character(spec$data$positive %||% lv[2])
  if (!case %in% lv) {
    stop(sprintf("--positive '%s' is not a level of '%s' (%s)", case, y, paste(lv, collapse = ", ")))
  }
  control <- setdiff(lv, case)
  n_case <- sum(d$outcome == case)
  if (n_case < 2 || nrow(d) - n_case < 2) stop("calibration needs at least 2 events and 2 non-events")
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing outcome or prediction, so every model is read on the same subjects.", sum(!keep)),
    sprintf("Kept %d of %d row(s): %d event(s) (%s) and %d non-event(s) (%s).",
            nrow(d), n0, n_case, case, nrow(d) - n_case, control),
    sprintf("Predictions were clipped to [%g, 1 - %g] before the logit transform.", .cal_clip, .cal_clip))

  fits <- lapply(models, .cal_model, d = d, case = case)
  stats_df <- do.call(rbind, lapply(fits, .cal_row, case = case))
  rownames(stats_df) <- NULL

  pal <- pub_palette(length(models), spec$appearance$palette)
  n_models <- length(models)
  top <- max(vapply(fits, .cal_spike_max, numeric(1)))
  spikes <- do.call(rbind, lapply(seq_len(n_models), .cal_spikes, fits = fits, top = top))
  pts <- do.call(rbind, lapply(seq_len(n_models), .cal_points, fits = fits, n_models = n_models))
  # A Wilson bound can pass 1 or 0 by a rounding error; the y scale would then
  # drop the bar, so clamp the bounds to the probability range.
  pts$hi <- pmin(pts$hi, 1)
  pts$lo <- pmax(pts$lo, 0)
  sm <- do.call(rbind, lapply(fits, .cal_smooth))
  spikes$model <- factor(spikes$model, levels = models)
  pts$model <- factor(pts$model, levels = models)
  sm$model <- factor(sm$model, levels = models)
  bases <- data.frame(y = .cal_band_base(seq_len(n_models)))
  labels <- data.frame(y = c(bases$y + 0.045, bases$y - 0.045),
                       label = rep(c("events", "non-events"), each = n_models))
  ylow <- .cal_band_bottom(n_models)
  dims <- .cal_dims(n_models, vapply(fits, .cal_legend_label, character(1)))
  sub <- sprintf("%d subjects, %d events (%s = %s)", nrow(d), n_case, y, case)

  p <- ggplot2::ggplot() +
    ggplot2::annotate("segment", x = 0, xend = 1, y = 0, yend = 1, linetype = "dashed",
                      colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_segment(data = bases, ggplot2::aes(x = 0, xend = 1, y = y, yend = y),
                          colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_segment(data = spikes, ggplot2::aes(x = x, xend = x, y = y0, yend = y1,
                                                      colour = model),
                          linewidth = 0.45, show.legend = FALSE) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = -0.025, y = y, label = label),
                       hjust = 1, size = 2.6, colour = "grey35") +
    ggplot2::geom_path(data = sm, ggplot2::aes(x = x, y = y, colour = model), linewidth = 0.8) +
    ggplot2::geom_errorbar(data = pts, ggplot2::aes(x = x, ymin = lo, ymax = hi, colour = model),
                           width = 0.015, linewidth = 0.45, show.legend = FALSE) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = model), size = 2,
                        show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, models), breaks = models,
                                 labels = vapply(fits, .cal_legend_label, character(1))) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                                labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                                expand = ggplot2::expansion(add = 0.02), oob = scales::oob_keep) +
    ggplot2::scale_y_continuous(limits = c(ylow, 1), breaks = seq(0, 1, 0.25),
                                labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
                                expand = ggplot2::expansion(add = c(0, 0.03)), oob = scales::oob_keep) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::guides(x = ggplot2::guide_axis(cap = "both"), y = ggplot2::guide_axis(cap = "both"),
                    colour = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::labs(x = spec$appearance$x_label %||% "Predicted probability",
                  y = spec$appearance$y_label %||% "Observed proportion",
                  colour = NULL, title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 9),
                   legend.key.height = ggplot2::unit(0.18, "in"),
                   legend.key.spacing.y = ggplot2::unit(0.02, "in"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   axis.title.y = ggplot2::element_text(hjust = .cal_ylab_pos(ylow)),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  test_meta <- list(
    name = "Calibration (Hosmer-Lemeshow test, calibration intercept and slope)",
    case_level = case, control_level = control, n = nrow(d), n_events = n_case,
    groups_requested = .cal_g,
    models = stats::setNames(lapply(fits, .cal_meta), models))
  resolved <- list(family = "calibration", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .cal_methods(fits, case, control, n_case, nrow(d), y),
       build_script = list(fn = .cal_emit_script,
                           args = list(spec = spec, models = models, case = case, dims = dims)),
       label = slugify(paste(c(y, models), collapse = "_")),
       df_used = d, clean_steps = clean_steps,
       groups = do.call(rbind, lapply(fits, .cal_groups_frame)),
       width = dims$width, height = dims$height)
}

.cal_groups_frame <- function(r) cbind(model = r$model, r$groups, stringsAsFactors = FALSE)

.cal_meta <- function(r) {
  list(n = r$n, n_events = r$n_events, brier = r$brier,
       intercept = r$intercept$estimate, slope = r$slope$estimate,
       hosmer_lemeshow = list(statistic = r$hl$statistic, df = r$hl$df, p_value = r$hl$p_value,
                              groups = r$hl$groups),
       oe_ratio = r$oe_ratio)
}

# The y title sits at the middle of the 0 to 1 part of the axis, not of the
# whole panel, which also holds the histogram bands.
.cal_ylab_pos <- function(ylow) (0.5 - ylow) / (1.03 - ylow)

# Figure size: the panel is square in data units, so the height grows with the
# histogram bands and the legend lines.
.cal_dims <- function(n_models, labels) {
  span <- 1.03 - .cal_band_bottom(n_models)
  lines <- sum(1 + grepl("\n", labels, fixed = TRUE))
  list(width = 5.2, height = 4.3 * span / 1.04 + 1.35 + 0.2 * lines)
}

.cal_methods <- function(fits, case, control, n_case, n, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions(c("ggplot2"))
  per <- paste(vapply(fits, .cal_model_text, character(1)), collapse = "; ")
  paste0(
    sprintf("The calibration of the predicted probabilities was assessed against the observed outcome %s = %s (events, n = %d) versus %s (non-events, n = %d). ",
            y, case, n_case, control, n - n_case),
    "For each model, the Brier score is the mean squared difference between the prediction and the outcome. ",
    "The calibration intercept (calibration-in-the-large) is the intercept of a logistic regression of the outcome with logit(p) as an offset, and the calibration slope is the coefficient of logit(p) in a logistic regression, each with a 95% Wald CI; ideal values are 0 and 1. ",
    sprintf("The Hosmer-Lemeshow test grouped the subjects at the deciles of the prediction (g = %d; tied deciles merge, df = groups - 2). ", .cal_g),
    "The figure shows the observed proportion against the mean prediction in each group with a 95% Wilson CI, a lowess smooth of the outcome on the prediction (f = 2/3, without robustness iterations, which would treat the events of a binary outcome as outliers), and the spike histogram of the predictions. ",
    per, ".",
    sprintf(" Computed in R %s with stats::glm and stats::lowess and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.cal_model_text <- function(r) {
  sprintf("%s: Brier score %.3f, calibration intercept %.2f (95%% CI %.2f to %.2f), slope %.2f (95%% CI %.2f to %.2f), Hosmer-Lemeshow chi-square(%d) = %.2f, P %s, observed/expected %.2f",
          r$model, r$brier, r$intercept$estimate, r$intercept$ci_low, r$intercept$ci_high,
          r$slope$estimate, r$slope$ci_low, r$slope$ci_high, as.integer(r$hl$df),
          r$hl$statistic, .cal_p_text(r$hl$p_value), r$oe_ratio)
}

.cal_emit_script <- function(spec, models, case, dims, in_name, fig_stub) {
  y <- spec$data$y
  mk <- paste0("c(", paste(sprintf('"%s"', models), collapse = ", "), ")")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ResourceSelection); library(ggplot2) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE, stringsAsFactors = FALSE)', in_name),
    sprintf("models <- %s", mk),
    sprintf('raw <- raw[!is.na(raw[["%s"]]) & complete.cases(raw[, models, drop = FALSE]), ]', y),
    sprintf('yb <- as.integer(as.character(raw[["%s"]]) == "%s")', y, case),
    "curves <- NULL",
    "for (m in models) {",
    "  p <- raw[[m]]",
    "  lp <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))",
    "  cat(m, \"Brier\", mean((p - yb)^2), \"\\n\")",
    "  print(coef(summary(glm(yb ~ 1 + offset(lp), family = binomial))))",
    "  print(coef(summary(glm(yb ~ lp, family = binomial))))",
    "  print(hoslem.test(yb, p, g = 10))",
    "  sm <- lowess(p, yb, f = 2 / 3, iter = 0, delta = 0.01 * diff(range(p)))",
    "  curves <- rbind(curves, data.frame(model = m, x = sm$x, y = sm$y))",
    "}",
    "g <- ggplot(curves, aes(x, y, colour = model)) + geom_abline(linetype = \"dashed\") + geom_path() +",
    "  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +",
    '  labs(x = "Predicted probability", y = "Observed proportion") + theme_classic(base_size = 13)',
    sprintf('ggsave("%s.pdf", g, width = %.2f, height = %.2f)', fig_stub, dims$width, dims$height),
    'cat("Reproduced figure.\\n")'
  )
}
