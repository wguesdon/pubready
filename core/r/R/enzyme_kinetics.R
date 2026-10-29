# Recipe: enzyme kinetics. The Michaelis-Menten equation
#   v = Vmax * S / (Km + S)
# is fitted to each curve (--group, for example one curve per inhibitor
# concentration) by nonlinear least squares. Each parameter carries its
# asymptotic standard error and the 95% CI estimate +/- qt(0.975, n - 2) * SE,
# which is the interval GraphPad Prism reports. With two or more curves, two
# extra-sum-of-squares F tests compare one shared Km (a Vmax per curve) and one
# shared Vmax (a Km per curve) with the separate fits.
#
# Both engines start from the same values and minimise the same residual sum of
# squares, so they converge to the same estimates.

.enz_null <- function(e) NULL

.enz_mean_at <- function(u, s, v) mean(v[s == u])

.enz_control <- function() stats::nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)

# Starting values: Vmax0 is the largest rate, and Km0 is the positive substrate
# concentration whose mean rate is nearest Vmax0 / 2.
.enz_start <- function(s, v) {
  us <- sort(unique(s[s > 0]))
  means <- vapply(us, .enz_mean_at, numeric(1), s = s, v = v)
  vmax <- max(v)
  list(vmax = vmax, km = us[which.min(abs(means - vmax / 2))])
}

.enz_predict <- function(s, vmax, km) vmax * s / (km + s)

# Fit one curve: Gauss-Newton first, then the PORT routine if it fails.
.enz_fit_curve <- function(s, v) {
  d <- data.frame(s = s, v = v)
  st <- .enz_start(s, v)
  algorithm <- "default"
  fit <- tryCatch(stats::nls(v ~ vmax * s / (km + s), data = d, start = st,
                             control = .enz_control()), error = .enz_null)
  if (is.null(fit)) {
    algorithm <- "port"
    fit <- tryCatch(stats::nls(v ~ vmax * s / (km + s), data = d, start = st,
                               control = .enz_control(), algorithm = "port"),
                    error = .enz_null)
  }
  list(fit = fit, start = st, algorithm = algorithm)
}

.enz_curve_row <- function(g, s, v, fit) {
  cf <- summary(fit)$coefficients
  rss <- sum(stats::residuals(fit)^2)
  dfr <- length(v) - 2
  tq <- stats::qt(0.975, dfr)
  vm <- cf["vmax", "Estimate"]; vse <- cf["vmax", "Std. Error"]
  km <- cf["km", "Estimate"]; kse <- cf["km", "Std. Error"]
  data.frame(group = g, n = length(v), concentrations = length(unique(s)),
             vmax = vm, vmax_se = vse, vmax_ci_low = vm - tq * vse, vmax_ci_high = vm + tq * vse,
             km = km, km_se = kse, km_ci_low = km - tq * kse, km_ci_high = km + tq * kse,
             residual_sd = sqrt(rss / dfr), r_squared = 1 - rss / sum((v - mean(v))^2),
             rss = rss, df_resid = dfr, stringsAsFactors = FALSE)
}

# One shared-parameter model fitted to all curves, from the separate estimates.
# shared = "km": one Km and a Vmax per curve; shared = "vmax": the reverse.
.enz_shared_fit <- function(d, rows, shared) {
  d$gi <- match(d$group, rows$group)
  if (shared == "km") {
    st <- list(vmax = rows$vmax, km = stats::median(rows$km))
    fml <- v ~ vmax[gi] * s / (km + s)
  } else {
    st <- list(vmax = stats::median(rows$vmax), km = rows$km)
    fml <- v ~ vmax * s / (km[gi] + s)
  }
  fit <- tryCatch(stats::nls(fml, data = d, start = st, control = .enz_control()),
                  error = .enz_null)
  if (is.null(fit)) {
    fit <- tryCatch(stats::nls(fml, data = d, start = st, control = .enz_control(),
                               algorithm = "port"), error = .enz_null)
  }
  fit
}

# Extra-sum-of-squares F test of one shared parameter against separate fits.
.enz_shared_test <- function(d, rows, shared) {
  fit <- .enz_shared_fit(d, rows, shared)
  if (is.null(fit)) return(NULL)
  rss_sh <- sum(stats::residuals(fit)^2)
  rss_sep <- sum(rows$rss); df_sep <- sum(rows$df_resid)
  df1 <- nrow(rows) - 1
  fval <- ((rss_sh - rss_sep) / df1) / (rss_sep / df_sep)
  cf <- stats::coef(fit)
  list(statistic = fval, df1 = df1, df2 = df_sep,
       p_value = stats::pf(fval, df1, df_sep, lower.tail = FALSE),
       shared_value = unname(cf[[if (shared == "km") "km" else "vmax"]]),
       rss_shared = rss_sh, rss_separate = rss_sep)
}

# Mean and standard error of the rate at each concentration of each curve.
.enz_summary <- function(d) {
  key <- paste(d$group, d$s, sep = "\r")
  first <- !duplicated(key)
  out <- d[first, c("group", "s"), drop = FALSE]
  out$mean <- as.numeric(tapply(d$v, key, mean)[key[first]])
  out$sem <- as.numeric(tapply(d$v, key, .enz_sem)[key[first]])
  out$n <- as.integer(tapply(d$v, key, length)[key[first]])
  out
}

.enz_sem <- function(v) if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_

.enz_curve_grid <- function(i, rows, smax) {
  grid <- seq(0, smax, length.out = 200)
  data.frame(group = rows$group[i], s = grid, v = .enz_predict(grid, rows$vmax[i], rows$km[i]),
             stringsAsFactors = FALSE)
}

.enz_fmt_test <- function(name, t) {
  if (is.null(t)) return(sprintf("%s: not computed", name))
  p <- if (t$p_value < 0.001) "P < 0.001" else sprintf("P = %.2g", t$p_value)
  sprintf("%s: F(%d, %d) = %.2f, %s", name, t$df1, t$df2, t$statistic, p)
}

.enz_legend_label <- function(i, rows) {
  sprintf("%s\nVmax %s, Km %s", .est_wrap(rows$group[i], 18), sprintf("%.3g", rows$vmax[i]),
          sprintf("%.3g", rows$km[i]))
}

recipe_enzyme_kinetics <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; grp <- spec$data$group
  if (is.null(x) || is.null(y)) {
    stop("enzyme_kinetics needs --x (the substrate concentration) and --y (the rate)")
  }
  has_group <- !is.null(grp) && nzchar(grp)
  for (cn in c(x, y, if (has_group) grp)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(s = suppressWarnings(as.numeric(df[[x]])),
                  v = suppressWarnings(as.numeric(df[[y]])),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  ok <- is.finite(d$s) & is.finite(d$v) & !is.na(d$group)
  n_missing <- sum(!ok)
  d <- d[ok, , drop = FALSE]
  neg <- d$s < 0
  d <- d[!neg, , drop = FALSE]
  clean_steps <- c(
    if (n_missing > 0) sprintf("Dropped %d row(s) with a missing or non-numeric substrate, rate or group.", n_missing),
    if (sum(neg) > 0) sprintf("Dropped %d row(s) with a negative substrate concentration.", sum(neg)),
    sprintf("Kept %d of %d row(s).", nrow(d), n0))

  g_lv <- unique(d$group)
  rows <- list(); starts <- list(); algos <- list()
  for (g in g_lv) {
    sel <- d$group == g
    n_conc <- length(unique(d$s[sel & d$s > 0]))
    if (n_conc < 3 || sum(sel) < 3) {
      stop(sprintf("enzyme_kinetics needs at least 3 distinct positive substrate concentrations and 3 points per curve; curve '%s' has %d and %d",
                   g, n_conc, sum(sel)))
    }
    fr <- .enz_fit_curve(d$s[sel], d$v[sel])
    if (is.null(fr$fit)) stop(sprintf("the Michaelis-Menten fit did not converge for curve '%s'", g))
    rows[[g]] <- .enz_curve_row(g, d$s[sel], d$v[sel], fr$fit)
    starts[[g]] <- fr$start
    algos[[g]] <- fr$algorithm
  }
  stats_df <- do.call(rbind, rows)
  rownames(stats_df) <- NULL

  shared_km <- shared_vmax <- NULL
  if (length(g_lv) >= 2) {
    shared_km <- .enz_shared_test(d, stats_df, "km")
    shared_vmax <- .enz_shared_test(d, stats_df, "vmax")
  }

  ap <- spec$appearance
  pal <- stats::setNames(pub_palette(length(g_lv), ap$palette), g_lv)
  leg <- stats::setNames(vapply(seq_along(g_lv), .enz_legend_label, character(1), rows = stats_df), g_lv)
  # Each legend key is as tall as the longest label, so multi-line labels never overlap.
  key_lines <- max(lengths(regmatches(leg, gregexpr("\n", leg)))) + 1
  curves <- do.call(rbind, lapply(seq_along(g_lv), .enz_curve_grid, rows = stats_df,
                                  smax = max(d$s)))
  pts <- .enz_summary(d)
  replicated <- any(pts$n > 1)
  curves$group <- factor(curves$group, levels = g_lv)
  pts$group <- factor(pts$group, levels = g_lv)
  raw <- d
  raw$group <- factor(raw$group, levels = g_lv)

  p <- ggplot2::ggplot()
  if (replicated) {
    p <- p +
      ggplot2::geom_point(data = raw, ggplot2::aes(x = s, y = v, colour = group), size = 1.1,
                          alpha = 0.35) +
      ggplot2::geom_errorbar(data = pts[pts$n > 1, , drop = FALSE],
                             ggplot2::aes(x = s, ymin = mean - sem, ymax = mean + sem, colour = group),
                             width = 0, linewidth = 0.5)
  }
  p <- p +
    ggplot2::geom_line(data = curves, ggplot2::aes(x = s, y = v, colour = group), linewidth = 0.8) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = s, y = mean, colour = group), size = 2)
  subtitle <- if (length(g_lv) >= 2) {
    paste(.enz_fmt_test("Shared Km", shared_km), .enz_fmt_test("Shared Vmax", shared_vmax), sep = "\n")
  } else {
    r <- stats_df[1, ]
    sprintf("Vmax = %s (95%% CI %s to %s)\nKm = %s (95%% CI %s to %s)",
            sprintf("%.3g", r$vmax), sprintf("%.3g", r$vmax_ci_low), sprintf("%.3g", r$vmax_ci_high),
            sprintf("%.3g", r$km), sprintf("%.3g", r$km_ci_low), sprintf("%.3g", r$km_ci_high))
  }
  p <- p +
    ggplot2::scale_colour_manual(values = pal, labels = leg, breaks = g_lv) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.04))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(x = ap$x_label %||% x, y = ap$y_label %||% y,
                  colour = if (has_group) grp else NULL, title = ap$title, subtitle = subtitle) +
    pub_theme(ap$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   legend.position = if (has_group) "right" else "none",
                   legend.key.height = ggplot2::unit(0.42 * key_lines, "cm"),
                   legend.key.spacing.y = ggplot2::unit(5, "pt"),
                   legend.title = ggplot2::element_text(size = 11),
                   legend.text = ggplot2::element_text(size = 9.5))

  test_meta <- list(
    name = "Michaelis-Menten fit", n_curves = length(g_lv),
    curves = stats::setNames(lapply(seq_len(nrow(stats_df)), .enz_meta_row, rows = stats_df), g_lv),
    shared_km_test = shared_km, shared_vmax_test = shared_vmax)
  df_used <- data.frame(group = d$group, substrate = d$s, rate = d$v, stringsAsFactors = FALSE)
  list(plot = p, stats = stats_df, test_meta = test_meta,
       resolved = list(family = "enzyme_kinetics", label = test_meta$name),
       methods = .enz_methods(stats_df, shared_km, shared_vmax, nrow(d), spec, x, y, grp, replicated),
       build_script = list(fn = .enz_emit_script,
                           args = list(spec = spec, starts = starts, algos = algos,
                                       has_group = has_group)),
       label = slugify(y), df_used = df_used, clean_steps = clean_steps,
       width = if (has_group) 6.6 else 4.8,
       height = max(4.4, 1.4 + 0.165 * length(g_lv) * key_lines))
}

.enz_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(vmax = r$vmax, vmax_ci_low = r$vmax_ci_low, vmax_ci_high = r$vmax_ci_high,
       km = r$km, km_ci_low = r$km_ci_low, km_ci_high = r$km_ci_high,
       r_squared = r$r_squared)
}

.enz_methods <- function(rows, shared_km, shared_vmax, n_points, spec, x, y, grp, replicated) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  x_lab <- spec$appearance$x_label %||% x
  multi <- nrow(rows) > 1
  fits <- paste(sprintf("%sVmax = %s (95%% CI %s to %s), Km = %s (95%% CI %s to %s), R2 = %.3f",
                        if (multi) paste0(rows$group, ": ") else "",
                        sprintf("%.3g", rows$vmax), sprintf("%.3g", rows$vmax_ci_low),
                        sprintf("%.3g", rows$vmax_ci_high), sprintf("%.3g", rows$km),
                        sprintf("%.3g", rows$km_ci_low), sprintf("%.3g", rows$km_ci_high),
                        rows$r_squared), collapse = "; ")
  shared_txt <- if (multi) {
    paste0(" Extra-sum-of-squares F tests compared the separate fits with one shared Km (a Vmax for each curve) and with one shared Vmax (a Km for each curve): ",
           .enz_test_text("shared Km", shared_km), "; ", .enz_test_text("shared Vmax", shared_vmax), ".")
  } else ""
  paste0(
    sprintf("%s was modelled as a function of %s with the Michaelis-Menten equation v = Vmax S / (Km + S), fitted by nonlinear least squares to %s (n = %d points). ",
            cap_first(y_lab), x_lab,
            if (multi) sprintf("each level of %s (%d curves)", grp, nrow(rows)) else "one curve",
            n_points),
    "The starting values were the largest rate for Vmax and, for Km, the substrate concentration whose mean rate is nearest half of it. ",
    "Each 95% confidence interval is the estimate plus or minus t(0.975, n - 2) times its asymptotic standard error. ",
    .enz_points_text(replicated),
    fits, ".", shared_txt,
    sprintf(" Fitted in R %s with stats::nls and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.enz_test_text <- function(name, t) {
  if (is.null(t)) return(sprintf("%s not computed, because the shared fit did not converge", name))
  sprintf("%s F(%d, %d) = %.2f, P = %s", name, t$df1, t$df2, t$statistic, sprintf("%.2g", t$p_value))
}

.enz_emit_script <- function(spec, starts, algos, has_group, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y; grp <- spec$data$group
  start_lines <- vapply(names(starts), .enz_start_literal, character(1), starts = starts)
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(s = as.numeric(raw[["%s"]]), v = as.numeric(raw[["%s"]]),', x, y),
    sprintf('                group = %s)', if (has_group) sprintf('as.character(raw[["%s"]])', grp) else '"all"'),
    "d <- d[is.finite(d$s) & is.finite(d$v) & !is.na(d$group) & d$s >= 0, ]",
    "starts <- list(",
    paste0("  ", start_lines, c(rep(",", length(start_lines) - 1), "")),
    ")",
    sprintf("algos <- c(%s)", paste(sprintf('"%s" = "%s"', names(algos), unlist(algos)), collapse = ", ")),
    "ctrl <- nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)",
    "curves <- NULL",
    "for (g in names(starts)) {",
    "  fit <- nls(v ~ vmax * s / (km + s), data = d[d$group == g, ], start = starts[[g]],",
    "             control = ctrl, algorithm = algos[[g]])",
    "  print(summary(fit))",
    "  grid <- seq(0, max(d$s), length.out = 200)",
    "  curves <- rbind(curves, data.frame(group = g, s = grid, v = predict(fit, data.frame(s = grid))))",
    "}",
    "p <- ggplot() + geom_line(data = curves, aes(s, v, colour = group)) +",
    "  geom_point(data = d, aes(s, v, colour = group), alpha = 0.6) +",
    sprintf('  expand_limits(x = 0, y = 0) + labs(x = "%s", y = "%s") + theme_classic()',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    sprintf('ggsave("%s.pdf", p, width = %s, height = 4.4)', fig_stub, if (has_group) "6.6" else "4.8"),
    'cat("Reproduced figure.\\n")'
  )
}

.enz_start_literal <- function(g, starts) {
  s <- starts[[g]]
  sprintf('"%s" = list(vmax = %.17g, km = %.17g)', g, s$vmax, s$km)
}

# The figure draws the replicates and the error bars only when a concentration
# has replicates, so the methods paragraph says what each point shows.
.enz_points_text <- function(replicated) {
  if (replicated) {
    "Large points show the mean +/- SEM of the replicates at each concentration, and small points show each replicate. "
  } else "Each point is the single rate at its concentration. "
}
