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

# Summarise one fitted curve as one row of the stats table.
.drc_curve_row <- function(g, lx, y, fitres) {
  fit <- fitres$fit
  cf <- .drc_canonical(stats::coef(fit))
  se <- sqrt(diag(stats::vcov(fit)))
  dfr <- stats::df.residual(fit)
  rss <- sum(stats::residuals(fit)^2)
  tq <- stats::qt(0.975, dfr)
  data.frame(
    group = g, n = length(y), doses = length(unique(lx)),
    midpoint = if (cf[["hill"]] < 0) "IC50" else "EC50",
    bottom = cf[["bottom"]], top = cf[["top"]], hill = cf[["hill"]],
    log_ec50 = cf[["log_ec50"]], se_log_ec50 = unname(se[["log_ec50"]]),
    ec50 = 10^cf[["log_ec50"]],
    ec50_ci_low = 10^(cf[["log_ec50"]] - tq * se[["log_ec50"]]),
    ec50_ci_high = 10^(cf[["log_ec50"]] + tq * se[["log_ec50"]]),
    r_squared = 1 - rss / sum((y - mean(y))^2),
    rss = rss, df_resid = dfr, stringsAsFactors = FALSE)
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
  cf <- c(bottom = r$bottom, top = r$top, log_ec50 = r$log_ec50, hill = r$hill)
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

  g_lv <- unique(d$group)
  rows <- list(); starts <- list(); algos <- list()
  for (g in g_lv) {
    sel <- d$group == g
    if (length(unique(d$lx[sel])) < 4 || sum(sel) < 5) {
      stop(sprintf("dose_response needs at least 4 distinct doses and 5 points per curve; curve '%s' has %d dose(s) and %d point(s)",
                   g, length(unique(d$lx[sel])), sum(sel)))
    }
    fr <- .drc_fit_curve(d$lx[sel], d$y[sel])
    if (is.null(fr$fit)) {
      stop(sprintf("the four-parameter logistic did not converge for curve '%s'", g))
    }
    rows[[g]] <- .drc_curve_row(g, d$lx[sel], d$y[sel], fr)
    starts[[g]] <- fr$start
    algos[[g]] <- fr$algorithm
  }
  stats_df <- do.call(rbind, rows)
  rownames(stats_df) <- NULL

  shared <- if (length(g_lv) >= 2) .drc_shared_test(d, stats_df) else NULL

  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = length(g_lv))
  mid_name <- if (all(stats_df$midpoint == "IC50")) "IC50" else "EC50"
  leg <- sprintf("%s (%s = %s)", g_lv, stats_df$midpoint,
                 sprintf("%.3g", stats_df$ec50))
  names(pal) <- g_lv; names(leg) <- g_lv

  curves <- do.call(rbind, lapply(g_lv, .drc_curve_grid, rows = stats_df, d = d))
  pts <- .drc_dose_summary(d)
  pts$dose <- 10^pts$lx
  replicated <- any(pts$n > 1)
  curves$group <- factor(curves$group, levels = g_lv)
  pts$group <- factor(pts$group, levels = g_lv)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = dose, y = y, colour = group),
                       linewidth = 0.8)
  if (replicated) {
    p <- p +
      ggplot2::geom_errorbar(data = pts[pts$n > 1, , drop = FALSE],
                             ggplot2::aes(x = dose, ymin = mean - sem, ymax = mean + sem,
                                          colour = group), width = 0.05, linewidth = 0.4) +
      ggplot2::geom_point(data = pts, ggplot2::aes(x = dose, y = mean, colour = group), size = 2)
  } else {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = dose, y = mean, colour = group), size = 2)
  }
  subtitle <- if (!is.null(shared)) {
    sprintf("Shared %s F test: F(%d, %d) = %.2f, P = %s", mid_name, shared$df1, shared$df2,
            shared$statistic, sprintf("%.2g", shared$p_value))
  } else {
    sprintf("%s = %s (95%% CI %s to %s)", stats_df$midpoint[1],
            sprintf("%.3g", stats_df$ec50[1]),
            sprintf("%.3g", stats_df$ec50_ci_low[1]),
            sprintf("%.3g", stats_df$ec50_ci_high[1]))
  }
  p <- p +
    ggplot2::scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    ggplot2::scale_colour_manual(values = pal, labels = leg, breaks = g_lv) +
    ggplot2::labs(x = spec$appearance$x_label %||% x, y = spec$appearance$y_label %||% y,
                  colour = if (has_group) grp else NULL,
                  title = spec$appearance$title, subtitle = subtitle) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   legend.position = if (has_group) "right" else "none")

  test_meta <- list(
    name = "Four-parameter logistic fit",
    n_curves = length(g_lv),
    curves = stats::setNames(lapply(seq_len(nrow(stats_df)), .drc_meta_row, rows = stats_df), g_lv),
    shared_ec50_test = shared)
  resolved <- list(family = "dose_response", label = test_meta$name)
  methods <- .drc_methods(stats_df, shared, d, spec, x, y)
  build_script <- list(fn = .drc_emit_script,
                       args = list(spec = spec, starts = starts, algos = algos,
                                   has_group = has_group))

  df_used <- data.frame(group = d$group, dose = d$dose, response = d$y, stringsAsFactors = FALSE)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = slugify(y),
       df_used = df_used, clean_steps = clean_steps,
       width = if (has_group) 6.2 else 4.6, height = 4.2)
}

.drc_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(midpoint = r$midpoint, ec50 = r$ec50, ci_low = r$ec50_ci_low, ci_high = r$ec50_ci_high,
       hill = r$hill, bottom = r$bottom, top = r$top, r_squared = r$r_squared)
}

.drc_methods <- function(rows, shared, d, spec, x, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  x_lab <- spec$appearance$x_label %||% x
  fits <- paste(sprintf("%s %s = %s (95%% CI %s to %s, Hill slope %.2f, R2 = %.3f)",
                        if (nrow(rows) > 1) paste0(rows$group, ":") else "",
                        rows$midpoint,
                        sprintf("%.3g", rows$ec50),
                        sprintf("%.3g", rows$ec50_ci_low),
                        sprintf("%.3g", rows$ec50_ci_high),
                        rows$hill, rows$r_squared), collapse = "; ")
  shared_txt <- if (!is.null(shared)) {
    sprintf(" An extra-sum-of-squares F test compared one shared %s with one %s per curve, with the shared fit profiled over the midpoint and each Hill slope held between 0.05 and 20 in magnitude: F(%d, %d) = %.2f, P = %s.",
            rows$midpoint[1], rows$midpoint[1], shared$df1, shared$df2, shared$statistic,
            sprintf("%.2g", shared$p_value))
  } else if (nrow(rows) > 1) {
    " The shared-midpoint F test could not be computed, because the shared fit did not converge."
  } else ""
  paste0(
    sprintf("%s was modelled as a function of %s with a four-parameter logistic curve (variable Hill slope) fitted by nonlinear least squares on log10 dose (n = %d points, %d curve(s)). ",
            cap_first(y_lab), x_lab, nrow(d), nrow(rows)),
    "The 95% confidence interval of the midpoint was computed on the log10 scale from the standard error of the fit and back-transformed. ",
    fits, ".", shared_txt,
    sprintf(" Fitted in R %s with stats::nls and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.drc_emit_script <- function(spec, starts, algos, has_group, in_name, fig_stub) {
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
    "starts <- list(",
    paste0("  ", start_lines, c(rep(",", length(start_lines) - 1), "")),
    ")",
    sprintf("algos <- c(%s)", paste(sprintf('"%s" = "%s"', names(algos), unlist(algos)), collapse = ", ")),
    "fml <- y ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill))",
    "ctrl <- nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)",
    "curves <- NULL",
    "for (g in names(starts)) {",
    "  fit <- nls(fml, data = d[d$group == g, ], start = starts[[g]], control = ctrl,",
    "             algorithm = algos[[g]])",
    "  print(summary(fit))",
    "  grid <- seq(min(d$lx[d$group == g]), max(d$lx[d$group == g]), length.out = 200)",
    "  curves <- rbind(curves, data.frame(group = g, dose = 10^grid,",
    "                                     y = predict(fit, data.frame(lx = grid))))",
    "}",
    "p <- ggplot() + geom_line(data = curves, aes(dose, y, colour = group)) +",
    "  geom_point(data = d, aes(dose, y, colour = group), alpha = 0.6) +",
    sprintf('  scale_x_log10() + labs(x = "%s", y = "%s") + theme_classic()',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
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
