# Recipe: standard curve interpolation, for an ELISA, a protein assay or any
# assay read against standards. Rows with a concentration in --x are standards.
# Rows with an empty concentration are unknowns, named by --id (or by row
# number), whose concentration is read off the fitted curve.
#
# --test auto|linear|4pl picks the curve. auto fits the four-parameter logistic
# of dose_response (same parameterisation, same starting values and the same
# fitting routine) when the standards have at least 5 distinct positive
# concentrations and the fit converges, and a straight line otherwise. The 4PL is
# fitted on log10 concentration, so a zero standard (a blank) is left out of it.
#
# An unknown is interpolated by inverting the fitted curve. An unknown whose
# signal lies outside the signal range of the curve between the lowest and the
# highest standard is flagged "below range" or "above range" (its concentration
# would fall below the lowest or above the highest standard) and gets no
# concentration: nothing is extrapolated. Each standard level is also read back
# through the curve, and its recovery (back-calculated over nominal, in %) is
# the QC of the fit.

# Mean signal of one standard level. Module level so vapply can call it by name.
.stc_level_mean <- function(u, conc, sig) mean(sig[conc == u])
.stc_level_sd <- function(u, conc, sig) {
  v <- sig[conc == u]
  if (length(v) > 1) stats::sd(v) else NA_real_
}
.stc_level_n <- function(u, conc) sum(conc == u)

# The fitted signal at a concentration, for either model.
.stc_predict <- function(model, conc) {
  p <- model$parameters
  if (model$type == "4pl") {
    p$bottom + (p$top - p$bottom) / (1 + 10^((p$log_ec50 - log10(conc)) * p$hill))
  } else {
    p$intercept + p$slope * conc
  }
}

# The concentration of each signal, read off the curve, with its flag. The
# signal range is that of the curve between the lowest and the highest standard;
# outside it the flag says on which side of the standards the concentration lies.
.stc_invert <- function(model, sig) {
  ends <- .stc_predict(model, model$conc_range)
  lo <- min(ends); hi <- max(ends)
  inside <- sig >= lo & sig <= hi
  # A signal past the end of the lowest standard means a concentration below it.
  low_side <- (sig - ends[1]) * (ends[2] - ends[1]) < 0
  flag <- ifelse(inside, "in range", ifelse(low_side, "below range", "above range"))
  p <- model$parameters
  conc <- rep(NA_real_, length(sig))
  if (model$type == "4pl") {
    s <- sig[inside]
    conc[inside] <- 10^(p$log_ec50 - log10((p$top - p$bottom) / (s - p$bottom) - 1) / p$hill)
  } else {
    conc[inside] <- (sig[inside] - p$intercept) / p$slope
  }
  list(concentration = conc, flag = flag)
}

# Fit the 4PL to the positive standards with the dose_response routine.
.stc_fit_4pl <- function(conc, sig) {
  keep <- conc > 0
  lx <- log10(conc[keep]); y <- sig[keep]
  fr <- .drc_fit_curve(lx, y)
  if (is.null(fr$fit)) return(NULL)
  cf <- .drc_canonical(stats::coef(fr$fit))
  fitted <- .drc_predict(cf, lx)
  list(type = "4pl", label = "four-parameter logistic",
       parameters = list(bottom = cf[["bottom"]], top = cf[["top"]],
                         ec50 = 10^cf[["log_ec50"]], log_ec50 = cf[["log_ec50"]],
                         hill = cf[["hill"]]),
       r_squared = 1 - sum((y - fitted)^2) / sum((y - mean(y))^2),
       n_fit = length(y), conc_range = range(conc[keep]),
       start = fr$start, algorithm = fr$algorithm)
}

.stc_fit_linear <- function(conc, sig) {
  fit <- stats::lm(sig ~ conc)
  cf <- unname(stats::coef(fit))
  list(type = "linear", label = "linear",
       parameters = list(intercept = cf[1], slope = cf[2]),
       r_squared = summary(fit)$r.squared, n_fit = length(sig), conc_range = range(conc))
}

# The recovery of each positive standard level, read back through the curve.
.stc_recovery <- function(model, conc, sig) {
  lv <- sort(unique(conc[conc > 0]))
  m <- vapply(lv, .stc_level_mean, numeric(1), conc = conc, sig = sig)
  back <- .stc_invert(model, m)$concentration
  rec <- 100 * back / lv
  ok <- is.finite(rec)
  list(levels = lapply(seq_along(lv), .stc_recovery_row, lv = lv, m = m, back = back, rec = rec),
       n_levels = length(lv), n_read_back = sum(ok),
       min_pct = if (any(ok)) min(rec[ok]) else NA_real_,
       max_pct = if (any(ok)) max(rec[ok]) else NA_real_,
       n_within_80_120 = sum(ok & rec >= 80 & rec <= 120))
}

.stc_recovery_row <- function(i, lv, m, back, rec) {
  list(concentration = lv[i], mean_signal = m[i],
       back_calculated = if (is.finite(back[i])) back[i] else NA_real_,
       recovery_pct = if (is.finite(rec[i])) rec[i] else NA_real_)
}

# Ticks of the 1-2-5 sequence on a log axis, or of the 1-10 sequence when the
# range spans more than 2.5 decades, so the labels never crowd.
.stc_log_breaks <- function(lo, hi) {
  mult <- if (log10(hi / lo) > 2.5) 1 else c(1, 2, 5)
  seq125 <- as.vector(outer(mult, 10^(-6:9)))
  seq125 <- sort(seq125)
  inside <- seq125[seq125 >= lo * (1 - 1e-9) & seq125 <= hi * (1 + 1e-9)]
  if (length(inside) < 2) {
    below <- seq125[seq125 < lo]; above <- seq125[seq125 > hi]
    inside <- c(utils::tail(below, 1), inside, utils::head(above, 1))
  }
  inside
}

.stc_fmt <- function(v) sprintf("%.3g", v)

recipe_standard_curve <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; id <- spec$data$id
  if (is.null(x) || is.null(y)) {
    stop("standard_curve needs --x (the concentration, empty for an unknown) and --y (the signal)")
  }
  has_id <- !is.null(id) && nzchar(id)
  for (cn in c(x, y, if (has_id) id)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }
  method <- tolower(spec$test$method %||% "auto")
  if (!method %in% c("auto", "linear", "4pl")) {
    stop(sprintf("standard_curve --test must be auto, linear or 4pl; got '%s'", method))
  }

  n0 <- nrow(df)
  raw_conc <- df[[x]]
  blank <- is.na(raw_conc) | trimws(as.character(raw_conc)) == ""
  conc <- suppressWarnings(as.numeric(as.character(raw_conc)))
  sig <- suppressWarnings(as.numeric(df[[y]]))
  ids <- if (has_id) as.character(df[[id]]) else sprintf("row %d", seq_len(n0))
  bad_sig <- !is.finite(sig)
  bad_conc <- !blank & (!is.finite(conc) | conc < 0)
  keep <- !bad_sig & !bad_conc
  clean_steps <- c(
    if (sum(bad_sig) > 0) sprintf("Dropped %d row(s) with a missing or non-numeric signal in '%s'.", sum(bad_sig), y),
    if (sum(bad_conc & !bad_sig) > 0) sprintf("Dropped %d row(s) with a concentration that is not a number of zero or more.", sum(bad_conc & !bad_sig)),
    sprintf("Kept %d of %d row(s): %d standard(s) and %d unknown(s).", sum(keep),
            n0, sum(keep & !blank), sum(keep & blank)))
  is_std <- keep & !blank
  is_unk <- keep & blank
  s_conc <- conc[is_std]; s_sig <- sig[is_std]
  n_pos_levels <- length(unique(s_conc[s_conc > 0]))

  model <- NULL
  if (method %in% c("auto", "4pl")) {
    enough <- if (method == "auto") n_pos_levels >= 5 else n_pos_levels >= 4 && sum(s_conc > 0) >= 5
    if (method == "4pl" && !enough) {
      stop(sprintf("standard_curve --test 4pl needs at least 4 distinct positive standard concentrations and 5 positive standards; found %d and %d",
                   n_pos_levels, sum(s_conc > 0)))
    }
    if (enough) model <- .stc_fit_4pl(s_conc, s_sig)
    if (method == "4pl" && is.null(model)) stop("the four-parameter logistic did not converge on the standards")
  }
  if (is.null(model)) {
    if (length(unique(s_conc)) < 2 || length(s_conc) < 3) {
      stop(sprintf("standard_curve needs at least 3 standards at 2 or more concentrations; found %d at %d",
                   length(s_conc), length(unique(s_conc))))
    }
    model <- .stc_fit_linear(s_conc, s_sig)
    if (model$parameters$slope == 0) stop("the linear standard curve is flat, so no unknown can be read off it")
  }
  if (model$type == "4pl" && any(s_conc == 0)) {
    clean_steps <- c(clean_steps, sprintf("Left %d standard(s) at a concentration of zero out of the four-parameter logistic fit, because it is fitted on log10 concentration.",
                                          sum(s_conc == 0)))
  }

  inv <- .stc_invert(model, sig[is_unk])
  stats_df <- data.frame(id = ids[is_unk], signal = sig[is_unk],
                         concentration = inv$concentration, flag = inv$flag,
                         stringsAsFactors = FALSE)
  recovery <- .stc_recovery(model, s_conc, s_sig)
  counts <- list(n_unknowns = nrow(stats_df), n_in_range = sum(stats_df$flag == "in range"),
                 n_below = sum(stats_df$flag == "below range"),
                 n_above = sum(stats_df$flag == "above range"))
  ends <- .stc_predict(model, model$conc_range)
  test_meta <- c(
    list(name = "Standard curve interpolation", model = model$label,
         model_choice = if (method == "auto") "auto" else "requested",
         parameters = model$parameters, r_squared = model$r_squared,
         n_standards = length(s_conc), n_fitted = model$n_fit,
         n_levels = length(unique(s_conc)),
         concentration_range = model$conc_range, signal_range = sort(ends)),
    counts, list(recovery = recovery))

  fig <- .stc_plot(model, s_conc, s_sig, stats_df, counts, spec, x, y)
  df_used <- data.frame(id = ids[keep], role = ifelse(blank[keep], "unknown", "standard"),
                        concentration = conc[keep], signal = sig[keep], stringsAsFactors = FALSE)
  list(plot = fig, stats = stats_df, test_meta = test_meta,
       resolved = list(family = "standard_curve", label = test_meta$name),
       methods = .stc_methods(model, test_meta, method, spec, x, y),
       build_script = list(fn = .stc_emit_script, args = list(spec = spec, model = model)),
       label = slugify(y), df_used = df_used, clean_steps = clean_steps,
       width = 5.2, height = 4.8)
}

.stc_plot <- function(model, s_conc, s_sig, unk, counts, spec, x, y) {
  ap <- spec$appearance
  pal <- pub_palette(2, ap$palette)
  is_4pl <- model$type == "4pl"
  use <- if (is_4pl) s_conc > 0 else rep(TRUE, length(s_conc))
  lv <- sort(unique(s_conc[use]))
  lev <- data.frame(conc = lv,
                    mean = vapply(lv, .stc_level_mean, numeric(1), conc = s_conc, sig = s_sig),
                    sd = vapply(lv, .stc_level_sd, numeric(1), conc = s_conc, sig = s_sig),
                    n = vapply(lv, .stc_level_n, integer(1), conc = s_conc))
  rng <- model$conc_range
  grid <- if (is_4pl) 10^seq(log10(rng[1]), log10(rng[2]), length.out = 200) else
    seq(rng[1], rng[2], length.out = 200)
  curve <- data.frame(conc = grid, sig = .stc_predict(model, grid))
  reps <- data.frame(conc = s_conc[use], sig = s_sig[use])
  inr <- unk[unk$flag == "in range", , drop = FALSE]
  std_lab <- "Standard (mean \u00b1 SD)"
  unk_lab <- "Unknown (interpolated)"

  # Many unknowns would hide the standards and turn the drop lines into a band,
  # so a large set gets smaller diamonds and fainter lines.
  many <- nrow(inr) > 50
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = curve, ggplot2::aes(x = conc, y = sig), colour = pal[1], linewidth = 0.8)
  if (nrow(inr)) {
    p <- p +
      ggplot2::geom_segment(data = inr, ggplot2::aes(x = concentration, xend = concentration,
                                                     y = signal, yend = -Inf),
                            colour = "grey60", linewidth = if (many) 0.2 else 0.3,
                            linetype = "dashed", alpha = if (many) 0.3 else 0.8) +
      ggplot2::geom_point(data = inr, ggplot2::aes(x = concentration, y = signal,
                                                   shape = unk_lab, colour = unk_lab),
                          fill = "white", size = if (many) 1.4 else 2.2,
                          stroke = if (many) 0.45 else 0.7)
  }
  p <- p +
    ggplot2::geom_point(data = reps, ggplot2::aes(x = conc, y = sig), colour = pal[1], size = 1.2,
                        alpha = 0.4)
  if (any(lev$n > 1)) {
    p <- p + ggplot2::geom_errorbar(data = lev[lev$n > 1, , drop = FALSE],
                                    ggplot2::aes(x = conc, ymin = mean - sd, ymax = mean + sd),
                                    colour = pal[1], width = 0, linewidth = 0.5)
  }
  p <- p +
    ggplot2::geom_point(data = lev, ggplot2::aes(x = conc, y = mean, shape = std_lab, colour = std_lab),
                        size = 2.4)
  sub <- sprintf("%s fit, R\u00b2 = %.4f\n%d of %d unknowns interpolated; %d below range, %d above range",
                 cap_first(model$label), model$r_squared, counts$n_in_range, counts$n_unknowns,
                 counts$n_below, counts$n_above)
  p <- p +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 23), c(std_lab, unk_lab)),
                                breaks = c(std_lab, unk_lab), name = NULL) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, c(std_lab, unk_lab)),
                                 breaks = c(std_lab, unk_lab), name = NULL) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_extended(n = 6)) +
    ggplot2::labs(x = ap$x_label %||% x, y = ap$y_label %||% y, title = ap$title, subtitle = sub) +
    pub_theme(ap$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   legend.position = "bottom", legend.text = ggplot2::element_text(size = 10))
  if (is_4pl) {
    p <- p + ggplot2::scale_x_log10(breaks = .stc_log_breaks(rng[1], rng[2]),
                                    labels = scales::label_number(drop0trailing = TRUE, big.mark = ""))
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.04)))
  }
  p
}

.stc_methods <- function(model, meta, method, spec, x, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  x_lab <- spec$appearance$x_label %||% x
  p <- model$parameters
  rec <- meta$recovery
  choice <- if (method == "auto") {
    "The model was chosen automatically: the four-parameter logistic when the standards have at least five distinct positive concentrations and the fit converges, and a straight line otherwise. "
  } else "The model was set by the analyst. "
  fit_txt <- if (model$type == "4pl") {
    sprintf("The four-parameter logistic, signal = bottom + (top - bottom) / (1 + 10^((log10 EC50 - log10 concentration) x Hill slope)), was fitted by nonlinear least squares on log10 concentration to %d standards: bottom = %s, top = %s, EC50 = %s, Hill slope = %.3f, R2 = %.4f. ",
            model$n_fit, .stc_fmt(p$bottom), .stc_fmt(p$top), .stc_fmt(p$ec50), p$hill, model$r_squared)
  } else {
    sprintf("A straight line, signal = intercept + slope x concentration, was fitted by least squares to %d standards: intercept = %s, slope = %s, R2 = %.4f. ",
            model$n_fit, .stc_fmt(p$intercept), .stc_fmt(p$slope), model$r_squared)
  }
  rec_txt <- if (rec$n_read_back > 0) {
    sprintf("Read back through the curve, the standard levels had a recovery of %.1f%% to %.1f%% of the nominal concentration, and %d of %d levels fell within 80%% to 120%%. ",
            rec$min_pct, rec$max_pct, rec$n_within_80_120, rec$n_levels)
  } else "No standard level could be read back through the curve. "
  paste0(
    sprintf("Concentrations (%s) were interpolated from the signal (%s) on a standard curve of %d standards at %d concentrations. ",
            x_lab, y_lab, meta$n_standards, meta$n_levels),
    choice, fit_txt,
    sprintf("Each unknown was interpolated by inverting the fitted curve. An unknown whose signal lies outside the signal range of the curve between the lowest and the highest standard (%s to %s) was flagged below range or above range and given no concentration, so no value was extrapolated. ",
            .stc_fmt(meta$signal_range[1]), .stc_fmt(meta$signal_range[2])),
    sprintf("Of %d unknowns, %d were in range, %d below range and %d above range. ",
            meta$n_unknowns, meta$n_in_range, meta$n_below, meta$n_above),
    rec_txt,
    sprintf("Fitted in R %s with %s and drawn with ggplot2 %s.", rver,
            if (model$type == "4pl") "stats::nls" else "stats::lm", v$ggplot2))
}

.stc_emit_script <- function(spec, model, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  fit_lines <- if (model$type == "4pl") {
    c(sprintf("st <- list(bottom = %.17g, top = %.17g, log_ec50 = %.17g, hill = %g)",
              model$start$bottom, model$start$top, model$start$log_ec50, model$start$hill),
      "pos <- std[std$conc > 0, ]; pos$lx <- log10(pos$conc)",
      "fit <- nls(sig ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill)), data = pos, start = st,",
      sprintf('           control = nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12), algorithm = "%s")',
              model$algorithm),
      "cf <- coef(fit); if (cf[['bottom']] > cf[['top']]) cf <- c(bottom = cf[['top']], top = cf[['bottom']], log_ec50 = cf[['log_ec50']], hill = -cf[['hill']])",
      "f <- function(c) cf[['bottom']] + (cf[['top']] - cf[['bottom']]) / (1 + 10^((cf[['log_ec50']] - log10(c)) * cf[['hill']]))",
      "finv <- function(s) 10^(cf[['log_ec50']] - log10((cf[['top']] - cf[['bottom']]) / (s - cf[['bottom']]) - 1) / cf[['hill']])",
      "rng <- range(pos$conc)")
  } else {
    c("fit <- lm(sig ~ conc, data = std); cf <- coef(fit)",
      "f <- function(c) cf[[1]] + cf[[2]] * c",
      "finv <- function(s) (s - cf[[1]]) / cf[[2]]",
      "rng <- range(std$conc)")
  }
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(conc = suppressWarnings(as.numeric(as.character(raw[["%s"]]))), sig = as.numeric(raw[["%s"]]))', x, y),
    "d$row <- seq_len(nrow(d))",
    "std <- d[is.finite(d$conc) & is.finite(d$sig), ]",
    "unk <- d[is.na(d$conc) & is.finite(d$sig), ]",
    fit_lines,
    "print(summary(fit))",
    "ends <- f(rng)",
    "inside <- unk$sig >= min(ends) & unk$sig <= max(ends)",
    "unk$concentration <- ifelse(inside, finv(unk$sig), NA)",
    "print(unk)",
    "grid <- data.frame(conc = seq(rng[1], rng[2], length.out = 200)); grid$sig <- f(grid$conc)",
    "p <- ggplot() + geom_point(data = std, aes(conc, sig)) + geom_line(data = grid, aes(conc, sig)) +",
    "  geom_point(data = unk[inside, ], aes(concentration, sig), shape = 23) +",
    sprintf('  %slabs(x = "%s", y = "%s") + theme_classic()', if (model$type == "4pl") "scale_x_log10() + " else "",
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    sprintf('ggsave("%s.pdf", p, width = 5.2, height = 4.8)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
