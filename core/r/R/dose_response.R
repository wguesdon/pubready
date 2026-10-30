# Recipe: dose-response curves. A four-parameter logistic (4PL) is fitted to each
# curve by least squares,
#   y = bottom + (top - bottom) / (1 + 10^((log10(EC50) - log10(dose)) * hill)),
# which is the variable-slope model of GraphPad Prism. A positive Hill slope
# rises with the dose and a negative one falls, so a falling curve reports its
# midpoint as an IC50. The midpoint carries a 95% CI from the standard error of
# log10(EC50), which makes the interval asymmetric on the dose scale. With two or
# more curves (--group), an extra-sum-of-squares F test asks whether one shared
# EC50 fits the data as well as one EC50 per curve.
#
# Both engines compute the same starting values and minimise the same residual
# sum of squares, so they converge to the same estimates.
#
# A curve whose mean responses vary by less than 10% of the response range, or
# whose responses do not differ between the doses in a one-way ANOVA, is flat,
# and it gets no fit (.drc_flat_check). A curve whose fit fails, whose fitted
# span is below 10% of the response range, whose EC50 lies outside the tested
# doses or whose CI is not finite is reported as "not determined"
# (.drc_status_reason). Both are drawn as their mean +/- SEM without a curve,
# report EC50 NA and stay out of the shared-midpoint F test.

# Mean response at one dose. Module level so vapply can call it by name.
.drc_mean_at <- function(u, lx, y) mean(y[lx == u])

.drc_null <- function(e) NULL

# Starting values: the plateaus from the lowest and highest mean response, the
# midpoint at the dose whose mean is nearest half way, and a Hill slope of +1 or
# -1 from the sign of the correlation between log dose and response.
.drc_start <- function(lx, y) {
  ux <- sort(unique(lx))
  means <- vapply(ux, .drc_mean_at, numeric(1), lx = lx, y = y)
  bottom <- min(means); top <- max(means)
  hill <- if (stats::cor(lx, y) < 0) -1 else 1
  log_ec50 <- ux[which.min(abs(means - (bottom + top) / 2))]
  list(bottom = bottom, top = top, log_ec50 = log_ec50, hill = hill)
}

.drc_formula <- function() {
  y ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill))
}

.drc_control <- function() {
  stats::nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)
}

# Fit one curve. Gauss-Newton first, then the PORT routine if it fails.
.drc_fit_curve <- function(lx, y) {
  d <- data.frame(lx = lx, y = y)
  st <- .drc_start(lx, y)
  algorithm <- "default"
  fit <- tryCatch(stats::nls(.drc_formula(), data = d, start = st, control = .drc_control()),
                  error = .drc_null)
  if (is.null(fit)) {
    algorithm <- "port"
    fit <- tryCatch(stats::nls(.drc_formula(), data = d, start = st,
                               control = .drc_control(), algorithm = "port"),
                    error = .drc_null)
  }
  list(fit = fit, start = st, algorithm = algorithm)
}

# The 4PL is symmetric under swapping the plateaus and negating the slope. Report
# the form with bottom below top, so the sign of hill gives the direction.
.drc_canonical <- function(cf) {
  if (cf[["bottom"]] > cf[["top"]]) {
    cf <- c(bottom = cf[["top"]], top = cf[["bottom"]],
            log_ec50 = cf[["log_ec50"]], hill = -cf[["hill"]])
  }
  cf
}

.drc_predict <- function(cf, lx) {
  cf[["bottom"]] + (cf[["top"]] - cf[["bottom"]]) /
    (1 + 10^((cf[["log_ec50"]] - lx) * cf[["hill"]]))
}

# The share of the response range below which a fitted curve counts as flat.
.drc_min_span <- function() 0.10

.drc_flat_reason <- function() "the mean responses vary by less than 10% of the data range"

# The check before the fit. A curve whose mean responses vary by less than 10%
# of the range of the mean responses of all curves has no midpoint to estimate.
# The check uses no fitted value, so both engines decide the same way whatever
# their optimisers do on flat data. A single curve is never flat by this rule,
# because its own range is the reference.
.drc_is_flat <- function(lx, y, span_ref) {
  means <- vapply(unique(lx), .drc_mean_at, numeric(1), lx = lx, y = y)
  is.finite(span_ref) && span_ref > 0 && diff(range(means)) < .drc_min_span() * span_ref
}

.drc_anova_reason <- function() {
  "the responses do not differ between the doses (one-way ANOVA, P >= 0.05)"
}

# The P value of a one-way ANOVA of the responses of one curve across its doses,
# or NA when no dose has a replicate, as the ANOVA then has no residual degree
# of freedom. The sums of squares are written out, so the Python engine computes
# the same value.
.drc_dose_anova_p <- function(lx, y) {
  u <- unique(lx)
  k <- length(u); n <- length(y)
  if (k < 2 || n - k < 1) return(NA_real_)
  fitted <- vapply(u, .drc_mean_at, numeric(1), lx = lx, y = y)[match(lx, u)]
  ssw <- sum((y - fitted)^2)
  ssb <- sum((fitted - mean(y))^2)
  if (ssw == 0) return(if (ssb > 0) 0 else NA_real_)
  stats::pf((ssb / (k - 1)) / (ssw / (n - k)), k - 1, n - k, lower.tail = FALSE)
}

# The checks before the fit, in order: the 10% range rule of .drc_is_flat, then
# the one-way ANOVA across the doses, which also covers a single curve. Returns
# the reason of a flat curve, or NA when the curve goes to the fit.
.drc_flat_check <- function(lx, y, span_ref) {
  if (.drc_is_flat(lx, y, span_ref)) return(.drc_flat_reason())
  p <- .drc_dose_anova_p(lx, y)
  if (is.finite(p) && p >= 0.05) return(.drc_anova_reason())
  NA_character_
}

# The fit result of a flat curve: no fit, the starting values for its direction
# and the reason of the check before the fit.
.drc_flat_result <- function(lx, y, reason) {
  list(fit = NULL, start = .drc_start(lx, y), algorithm = NA_character_, flat = TRUE,
       reason = reason)
}

.drc_na_se <- function(e) c(bottom = NA_real_, top = NA_real_, log_ec50 = NA_real_, hill = NA_real_)

# Why a curve is "not determined", or NA when its fit stands. The rule is the
# same in both engines: the fit did not converge, the fitted span (top minus
# bottom) is below 10% of the range of the mean responses of all curves, the
# EC50 lies outside the tested doses, or the 95% CI of the EC50 is not finite.
.drc_status_reason <- function(converged, bottom, top, log_ec50, lx, ci_low, ci_high,
                               span_ref) {
  if (!converged) return("the fit did not converge")
  if (is.finite(span_ref) && span_ref > 0 && abs(top - bottom) < .drc_min_span() * span_ref) {
    return("the response span is below 10% of the data range")
  }
  if (!is.finite(log_ec50) || log_ec50 < min(lx) || log_ec50 > max(lx)) {
    return("the EC50 lies outside the tested doses")
  }
  if (!is.finite(ci_low) || !is.finite(ci_high) || ci_low <= 0) {
    return("the 95% CI of the EC50 is not finite")
  }
  NA_character_
}

# Summarise one curve as one row of the stats table. A curve whose fit fails the
# rule of .drc_status_reason keeps its sizes and direction but no fitted value.
.drc_curve_row <- function(g, lx, y, fitres, span_ref) {
  fit <- fitres$fit
  base <- data.frame(group = g, n = length(y), doses = length(unique(lx)),
                     midpoint_type = if (fitres$start$hill < 0) "IC50" else "EC50",
                     bottom = NA_real_, top = NA_real_, hill = NA_real_,
                     log_midpoint = NA_real_, se_log_midpoint = NA_real_, midpoint = NA_real_,
                     midpoint_ci_low = NA_real_, midpoint_ci_high = NA_real_,
                     r_squared = NA_real_, rss = NA_real_, df_resid = NA_integer_,
                     fit_status = "not determined",
                     status_reason = if (isTRUE(fitres$flat)) fitres$reason else "the fit did not converge",
                     stringsAsFactors = FALSE)
  if (is.null(fit)) return(base)
  cf <- .drc_canonical(stats::coef(fit))
  se <- tryCatch(suppressWarnings(sqrt(diag(stats::vcov(fit)))), error = .drc_na_se)
  dfr <- stats::df.residual(fit)
  rss <- sum(stats::residuals(fit)^2)
  tq <- stats::qt(0.975, dfr)
  ci_low <- 10^(cf[["log_ec50"]] - tq * se[["log_ec50"]])
  ci_high <- 10^(cf[["log_ec50"]] + tq * se[["log_ec50"]])
  reason <- .drc_status_reason(TRUE, cf[["bottom"]], cf[["top"]], cf[["log_ec50"]], lx,
                               ci_low, ci_high, span_ref)
  if (!is.na(reason)) {
    base$status_reason <- reason
    return(base)
  }
  data.frame(
    group = g, n = length(y), doses = length(unique(lx)),
    midpoint_type = if (cf[["hill"]] < 0) "IC50" else "EC50",
    bottom = cf[["bottom"]], top = cf[["top"]], hill = cf[["hill"]],
    log_midpoint = cf[["log_ec50"]], se_log_midpoint = unname(se[["log_ec50"]]),
    midpoint = 10^cf[["log_ec50"]], midpoint_ci_low = ci_low, midpoint_ci_high = ci_high,
    r_squared = 1 - rss / sum((y - mean(y))^2),
    rss = rss, df_resid = as.integer(dfr), fit_status = "fitted",
    status_reason = NA_character_, stringsAsFactors = FALSE)
}

# Residual sum of squares of one curve at a fixed log EC50 and Hill slope. With
# both fixed, the 4PL is linear in the two plateaus, so they come from an exact
# least-squares solve (variable projection). u is log(|hill|) and sgn its sign.
.drc_rss_fixed <- function(u, sgn, lx, y, log_ec50) {
  s <- 1 / (1 + 10^((log_ec50 - lx) * sgn * exp(u)))
  if (max(s) - min(s) < 1e-10) return(sum((y - mean(y))^2))
  sum(qr.resid(qr(cbind(1 - s, s)), y)^2)
}

# The best residual sum of squares of curve i at a fixed log EC50: Brent's method
# over log(|hill|) in [log 0.05, log 20], with the sign of the separate fit.
.drc_curve_rss_at <- function(i, log_ec50, d, rows) {
  sel <- d$group == rows$group[i]
  sgn <- if (rows$hill[i] < 0) -1 else 1
  stats::optimize(.drc_rss_fixed, c(log(0.05), log(20)), sgn = sgn, lx = d$lx[sel],
                  y = d$y[sel], log_ec50 = log_ec50, tol = 1e-10)$objective
}

# Residual sum of squares of the shared model at one fixed log EC50.
.drc_rss_at <- function(log_ec50, d, rows) {
  sum(vapply(seq_len(nrow(rows)), .drc_curve_rss_at, numeric(1), log_ec50 = log_ec50,
             d = d, rows = rows))
}

# Extra-sum-of-squares F test of one shared log EC50 against one per curve. The
# plateaus and the slopes stay free per curve in both models. The shared fit is a
# profile: the residual sum of squares over a grid of 31 log EC50 values across
# the doses, then Brent's method in the interval around the best grid point. A
# joint fit from one start can stop in a local minimum when the curves sit far
# apart; the profile finds the global one, and both engines run the same steps.
.drc_shared_test <- function(d, rows) {
  grid <- seq(min(d$lx), max(d$lx), length.out = 31)
  prof <- vapply(grid, .drc_rss_at, numeric(1), d = d, rows = rows)
  if (!any(is.finite(prof))) return(NULL)
  i <- which.min(prof)
  opt <- stats::optimize(.drc_rss_at, c(grid[max(1, i - 1)], grid[min(length(grid), i + 1)]),
                         d = d, rows = rows, tol = 1e-10)
  rss_sh <- min(opt$objective, prof[i])
  log_sh <- if (opt$objective <= prof[i]) opt$minimum else grid[i]
  rss_sep <- sum(rows$rss); df_sep <- sum(rows$df_resid)
  df1 <- nrow(rows) - 1
  fval <- ((rss_sh - rss_sep) / df1) / (rss_sep / df_sep)
  list(statistic = fval, df1 = df1, df2 = df_sep,
       p_value = stats::pf(fval, df1, df_sep, lower.tail = FALSE),
       shared_ec50 = 10^log_sh, rss_shared = rss_sh, rss_separate = rss_sep)
}

# Mean and standard error of the response at each dose of each curve.
.drc_dose_summary <- function(d) {
  key <- paste(d$group, d$lx, sep = "\r")
  first <- !duplicated(key)
  out <- d[first, c("group", "lx"), drop = FALSE]
  out$mean <- as.numeric(tapply(d$y, key, mean)[key[first]])
  out$sem <- as.numeric(tapply(d$y, key, .drc_sem)[key[first]])
  out$n <- as.integer(tapply(d$y, key, length)[key[first]])
  out
}

.drc_sem <- function(v) if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_

.drc_curve_grid <- function(g, rows, d) {
  r <- rows[rows$group == g, ]
  lxg <- d$lx[d$group == g]
  grid <- seq(min(lxg), max(lxg), length.out = 200)
  cf <- c(bottom = r$bottom, top = r$top, log_ec50 = r$log_midpoint, hill = r$hill)
  data.frame(group = g, dose = 10^grid, y = .drc_predict(cf, grid), stringsAsFactors = FALSE)
}

recipe_dose_response <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; grp <- spec$data$group
  if (is.null(x) || is.null(y)) {
    stop("dose_response needs --x (the dose or concentration) and --y (the response)")
  }
  has_group <- !is.null(grp) && nzchar(grp)
  for (cn in c(x, y, if (has_group) grp)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(dose = suppressWarnings(as.numeric(df[[x]])),
                  y = suppressWarnings(as.numeric(df[[y]])),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  ok <- is.finite(d$dose) & is.finite(d$y) & !is.na(d$group)
  n_missing <- sum(!ok)
  d <- d[ok, , drop = FALSE]
  nonpos <- d$dose <= 0
  d <- d[!nonpos, , drop = FALSE]
  clean_steps <- c(
    if (n_missing > 0) sprintf("Dropped %d row(s) with a missing or non-numeric dose, response or group.", n_missing),
    if (sum(nonpos) > 0) sprintf("Dropped %d row(s) with a dose of zero or below, because the dose axis is logarithmic.", sum(nonpos)),
    sprintf("Kept %d of %d row(s).", nrow(d), n0))
  d$lx <- log10(d$dose)
  # A response on a scale far from 1, such as 1e12 RLU, keeps nls from
  # converging. The fit runs on the response divided by a power of ten, and the
  # plateaus and the sums of squares go back to the original scale (F80).
  rs <- .drc_response_scale(d$y)
  y_orig <- d$y
  d$y <- d$y / rs

  g_lv <- unique(d$group)
  # The span rule compares each fitted span with the range of the mean
  # responses of all curves.
  pts <- .drc_dose_summary(d)
  span_ref <- diff(range(pts$mean))
  rows <- list(); starts <- list(); algos <- list()
  for (g in g_lv) {
    sel <- d$group == g
    if (length(unique(d$lx[sel])) < 4 || sum(sel) < 5) {
      stop(sprintf("dose_response needs at least 4 distinct doses and 5 points per curve; curve '%s' has %d dose(s) and %d point(s)",
                   g, length(unique(d$lx[sel])), sum(sel)))
    }
    flat_reason <- .drc_flat_check(d$lx[sel], d$y[sel], span_ref)
    fr <- if (!is.na(flat_reason)) {
      .drc_flat_result(d$lx[sel], d$y[sel], flat_reason)
    } else {
      .drc_fit_curve(d$lx[sel], d$y[sel])
    }
    rows[[g]] <- .drc_curve_row(g, d$lx[sel], d$y[sel], fr, span_ref)
    if (rows[[g]]$fit_status == "fitted") {
      starts[[g]] <- fr$start
      algos[[g]] <- fr$algorithm
    }
  }
  stats_df <- do.call(rbind, rows)
  rownames(stats_df) <- NULL
  det <- stats_df$fit_status == "fitted"
  det_lv <- stats_df$group[det]
  nd_lv <- stats_df$group[!det]

  # A curve that is not determined stays out of the shared-midpoint test.
  shared <- if (length(det_lv) >= 2) {
    .drc_shared_test(d[d$group %in% det_lv, , drop = FALSE], stats_df[det, , drop = FALSE])
  } else NULL
  stats_df <- .drc_unscale_rows(stats_df, rs)
  shared <- .drc_unscale_shared(shared, rs)
  d$y <- y_orig
  pts$mean <- pts$mean * rs
  pts$sem <- pts$sem * rs

  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = length(g_lv))
  mid_name <- .drc_mid_name(stats_df)
  leg <- ifelse(det, sprintf("%s (%s = %s)", g_lv, stats_df$midpoint_type, sig_text(stats_df$midpoint)),
                sprintf("%s (not determined)", g_lv))
  # A long curve name wraps, and the figure grows with the longest line, so the
  # legend does not squeeze the panel to nothing (F81). Python has the same rule.
  leg <- .est_wrap(leg, 28)
  names(pal) <- g_lv; names(leg) <- g_lv

  curves <- do.call(rbind, lapply(det_lv, .drc_curve_grid, rows = stats_df, d = d))
  pts$dose <- 10^pts$lx
  replicated <- any(pts$n > 1)
  pts$group <- factor(pts$group, levels = g_lv)

  p <- ggplot2::ggplot()
  if (!is.null(curves) && nrow(curves) > 0) {
    curves$group <- factor(curves$group, levels = g_lv)
    p <- p + ggplot2::geom_line(data = curves, ggplot2::aes(x = dose, y = y, colour = group),
                                linewidth = 0.8)
  }
  if (replicated) {
    p <- p +
      ggplot2::geom_errorbar(data = pts[pts$n > 1, , drop = FALSE],
                             ggplot2::aes(x = dose, ymin = mean - sem, ymax = mean + sem,
                                          colour = group), width = 0.05, linewidth = 0.4)
  }
  p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = dose, y = mean, colour = group), size = 2)
  subtitle <- .drc_subtitle(stats_df, shared, mid_name, has_group)
  p <- p +
    ggplot2::scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    ggplot2::scale_colour_manual(values = pal, labels = leg, breaks = g_lv,
                                 guide = ggplot2::guide_legend(
                                   override.aes = list(linetype = ifelse(det, 1, 0)))) +
    ggplot2::labs(x = spec$appearance$x_label %||% x, y = spec$appearance$y_label %||% y,
                  colour = if (has_group) grp else NULL,
                  title = spec$appearance$title, subtitle = subtitle) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   legend.position = if (has_group) "right" else "none",
                   legend.key.spacing.y = ggplot2::unit(4, "pt"))

  test_meta <- list(
    name = "Four-parameter logistic fit",
    n_curves = length(g_lv),
    n_determined = length(det_lv),
    not_determined = as.list(nd_lv),
    curves = stats::setNames(lapply(seq_len(nrow(stats_df)), .drc_meta_row, rows = stats_df), g_lv),
    shared_ec50_test = shared)
  resolved <- list(family = "dose_response", label = test_meta$name)
  methods <- .drc_methods(stats_df, shared, d, spec, x, y, replicated, sum(nonpos), rs)
  build_script <- list(fn = .drc_emit_script,
                       args = list(spec = spec, starts = starts, algos = algos,
                                   has_group = has_group, rs = rs))

  df_used <- data.frame(group = d$group, dose = d$dose, response = d$y, stringsAsFactors = FALSE)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = slugify(y),
       df_used = df_used, clean_steps = clean_steps,
       width = if (has_group) .drc_width(c(leg, grp)) else 4.6, height = 4.2)
}

# IC50 when every determined curve falls, EC50 otherwise.
.drc_mid_name <- function(rows) {
  mids <- rows$midpoint_type[rows$fit_status == "fitted"]
  if (length(mids) == 0) mids <- rows$midpoint_type
  if (all(mids == "IC50")) "IC50" else "EC50"
}

# The subtitle: the shared-midpoint test, or the midpoint of the one determined
# curve, and a second line that names the curves left out as not determined.
.drc_subtitle <- function(rows, shared, mid_name, has_group) {
  det <- rows$fit_status == "fitted"
  first <- if (!is.null(shared)) {
    sprintf("Shared %s F test: F(%d, %d) = %.2f, P = %s", mid_name, shared$df1, shared$df2,
            shared$statistic, sig_text(shared$p_value, 2))
  } else if (any(det)) {
    r <- rows[which(det)[1], ]
    sprintf("%s%s = %s (95%% CI %s to %s)", if (has_group) paste0(r$group, " ") else "",
            r$midpoint_type, sig_text(r$midpoint), sig_text(r$midpoint_ci_low),
            sig_text(r$midpoint_ci_high))
  } else if (has_group) {
    sprintf("%s not determined for any curve", mid_name)
  } else {
    sprintf("%s not determined: %s", mid_name, rows$status_reason[1])
  }
  nd <- rows$group[!det]
  if (!has_group || length(nd) == 0 || !any(det)) return(first)
  second <- if (length(nd) <= 2) {
    sprintf("Not determined, left out of the test: %s", paste(nd, collapse = ", "))
  } else {
    sprintf("%d curves not determined, left out of the test", length(nd))
  }
  paste(first, second, sep = "\n")
}

.drc_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(midpoint_type = r$midpoint_type, fit_status = r$fit_status, status_reason = r$status_reason,
       midpoint = r$midpoint, ci_low = r$midpoint_ci_low, ci_high = r$midpoint_ci_high,
       hill = r$hill, bottom = r$bottom, top = r$top, r_squared = r$r_squared)
}

.drc_methods <- function(rows, shared, d, spec, x, y, replicated, n_nonpos = 0, rs = 1) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  x_lab <- spec$appearance$x_label %||% x
  det <- rows$fit_status == "fitted"
  fr <- rows[det, , drop = FALSE]
  fits <- if (nrow(fr) > 0) {
    paste0(paste(sprintf("%s%s = %s (95%% CI %s to %s, Hill slope %.2f, R2 = %.3f)",
                         if (nrow(rows) > 1) paste0(fr$group, ": ") else "",
                         fr$midpoint_type,
                         sig_text(fr$midpoint),
                         sig_text(fr$midpoint_ci_low),
                         sig_text(fr$midpoint_ci_high),
                         fr$hill, fr$r_squared), collapse = "; "), ".")
  } else ""
  nd <- rows[!det, , drop = FALSE]
  nd_txt <- if (nrow(nd) > 0) {
    paste0(" ", paste(sprintf("%s was not determined, because %s.",
                              if (nrow(rows) > 1) paste0("The ", nd$midpoint_type, " of ", nd$group)
                              else paste0("The ", nd$midpoint_type),
                              nd$status_reason), collapse = " "),
           " A curve that is not determined is drawn as its mean +/- SEM without a fit",
           if (nrow(rows) > 1) " and is left out of the shared-midpoint F test." else ".")
  } else ""
  shared_txt <- if (!is.null(shared)) {
    sprintf(" An extra-sum-of-squares F test compared one shared %s with one %s per curve over the %d determined curves, with the shared fit profiled over the midpoint and each Hill slope held between 0.05 and 20 in magnitude: F(%d, %d) = %.2f, P = %s.",
            .drc_mid_name(rows), .drc_mid_name(rows), sum(det), shared$df1, shared$df2,
            shared$statistic, sig_text(shared$p_value, 2))
  } else if (nrow(rows) > 1 && sum(det) < 2) {
    " The shared-midpoint F test was not run, because fewer than two curves were determined."
  } else if (nrow(rows) > 1) {
    " The shared-midpoint F test could not be computed, because the shared fit did not converge."
  } else ""
  paste0(
    sprintf("%s was modelled as a function of %s with a four-parameter logistic curve (variable Hill slope) fitted by nonlinear least squares on log10 dose (n = %d points, %d curve(s)). ",
            cap_first(y_lab), x_lab, nrow(d), nrow(rows)),
    "The 95% confidence interval of the midpoint was computed on the log10 scale from the standard error of the fit and back-transformed. ",
    .drc_points_text(replicated),
    .drc_nonpos_text(n_nonpos),
    .drc_scale_text(rs),
    "A curve whose mean responses varied by less than 10% of the range of the mean responses of all curves, or whose responses did not differ between the doses in a one-way ANOVA (P >= 0.05), was not fitted, and its midpoint was reported as not determined. ",
    "A midpoint was also reported as not determined when the fit did not converge, when the fitted span (top minus bottom) was below 10% of the range of the mean responses, when the midpoint lay outside the tested doses, or when its 95% CI was not finite. ",
    fits, nd_txt, shared_txt,
    sprintf(" Fitted in R %s with stats::nls and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.drc_emit_script <- function(spec, starts, algos, has_group, in_name, fig_stub, rs = 1) {
  x <- spec$data$x; y <- spec$data$y; grp <- spec$data$group
  start_lines <- vapply(names(starts), .drc_start_literal, character(1), starts = starts)
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(dose = as.numeric(raw[["%s"]]), y = as.numeric(raw[["%s"]]),', x, y),
    sprintf('                group = %s)', if (has_group) sprintf('as.character(raw[["%s"]])', grp) else '"all"'),
    "d <- d[is.finite(d$dose) & is.finite(d$y) & d$dose > 0, ]",
    "d$lx <- log10(d$dose)",
    "# Only the curves with a determined midpoint are refitted; the others are",
    "# drawn as points, as in the figure.",
    "starts <- list(",
    paste0("  ", start_lines, c(rep(",", max(length(start_lines) - 1, 0)), "")[seq_along(start_lines)]),
    ")",
    sprintf("algos <- c(%s)", paste(sprintf('"%s" = "%s"', names(algos), unlist(algos)), collapse = ", ")),
    "# The fit runs on the response divided by a power of ten, as in the figure.",
    sprintf("rs <- %s; d$ys <- d$y / rs", format(rs, scientific = TRUE)),
    "fml <- ys ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill))",
    "ctrl <- nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)",
    "curves <- NULL",
    "for (g in names(starts)) {",
    "  fit <- nls(fml, data = d[d$group == g, ], start = starts[[g]], control = ctrl,",
    "             algorithm = algos[[g]])",
    "  print(summary(fit))",
    "  grid <- seq(min(d$lx[d$group == g]), max(d$lx[d$group == g]), length.out = 200)",
    "  curves <- rbind(curves, data.frame(group = g, dose = 10^grid,",
    "                                     y = predict(fit, data.frame(lx = grid)) * rs))",
    "}",
    "p <- ggplot() + geom_point(data = d, aes(dose, y, colour = group), alpha = 0.6) +",
    sprintf('  scale_x_log10() + labs(x = "%s", y = "%s") + theme_classic()',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    "if (!is.null(curves)) p <- p + geom_line(data = curves, aes(dose, y, colour = group))",
    sprintf('ggsave("%s.pdf", p, width = %s, height = 4.2)', fig_stub,
            if (has_group) "6.2" else "4.6"),
    'cat("Reproduced figure.\\n")'
  )
}

.drc_start_literal <- function(g, starts) {
  s <- starts[[g]]
  sprintf('"%s" = list(bottom = %.10g, top = %.10g, log_ec50 = %.10g, hill = %g)',
          g, s$bottom, s$top, s$log_ec50, s$hill)
}

# The figure draws error bars only when a dose has replicates, so the methods
# paragraph says what each point shows.
.drc_points_text <- function(replicated) {
  if (replicated) "Points show the mean +/- SEM of the replicates at each dose. "
  else "Each point is the single response at its dose. "
}

# The rows with a dose of zero or below have no place on the log dose axis. The
# methods text gives their count, as the data log does (F76).
.drc_nonpos_text <- function(n) {
  if (n > 0) {
    sprintf("%d row(s) with a dose of zero or below were left out of the fit and the figure, because the dose axis is logarithmic. ", n)
  } else ""
}

# The power of ten that the response is divided by before the fit: 10^e for a
# largest absolute response of the order 10^e when e is 3 or more or -3 or
# less, else 1 (F80). The Python engine has the same rule.
.drc_response_scale <- function(y) {
  m <- max(abs(y))
  if (!is.finite(m) || m == 0) return(1)
  e <- floor(log10(m))
  if (abs(e) >= 3) 10^e else 1
}

# The plateaus and the residual sum of squares of the stats table back on the
# scale of the response.
.drc_unscale_rows <- function(rows, rs) {
  rows$bottom <- rows$bottom * rs
  rows$top <- rows$top * rs
  rows$rss <- rows$rss * rs^2
  rows
}

.drc_unscale_shared <- function(shared, rs) {
  if (is.null(shared)) return(NULL)
  shared$rss_shared <- shared$rss_shared * rs^2
  shared$rss_separate <- shared$rss_separate * rs^2
  shared
}

.drc_scale_text <- function(rs) {
  if (rs == 1) return("")
  sprintf("The responses were divided by %s for the fit, so that the optimiser works near 1, and the plateaus are reported on the original scale. ",
          sprintf("%g", rs))
}

# The width of a figure with a legend: the 4.6 inch panel and a legend column
# that grows with the longest line of the wrapped labels and the legend title,
# at least 1.6 inches (F81). The Python engine has the same rule.
.drc_width <- function(labels) {
  longest <- max(nchar(unlist(strsplit(labels, "\n", fixed = TRUE))))
  4.6 + max(1.6, 0.2 + 0.065 * longest)
}
