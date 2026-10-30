# Recipe: cumulative incidence of competing risks.
# Input is one row per subject: the follow-up time (--time), the status
# (--event: 0 = censored, any other value = the event type, or cause, that ended
# the follow-up) and an optional grouping column (--x). The cause of interest is
# --cause (default 1). Kaplan-Meier overestimates the probability of an event
# when another event can prevent it, so each cause gets its Aalen-Johansen
# cumulative incidence function (CIF) in each group, from cmprsk::cuminc, and the
# groups are compared for each cause with Gray's k-sample test (rho = 0).
#
# The figure has one panel per cause, the cause of interest first, with a step
# curve per group and the number at risk below. A subject is at risk at time t
# when its follow-up time is t or more. The x breaks, the three tabulated time
# points and the y range come from the rules below, which the Python engine
# restates, so both engines draw the same axes.

# A step of 1, 2, 2.5 or 5 times a power of ten that cuts `span` into about
# `target` intervals.
.cif_nice_step <- function(span, target) {
  if (!is.finite(span) || span <= 0) return(1)
  raw <- span / target
  mag <- 10^floor(log10(raw))
  for (m in c(1, 2, 2.5, 5, 10)) {
    if (m * mag >= raw * (1 - 1e-9)) return(m * mag)
  }
  10 * mag
}

# The x breaks: multiples of the step from 0 up to the largest follow-up time.
.cif_breaks <- function(tmax) {
  step <- .cif_nice_step(tmax, 6)
  list(step = step, breaks = (0:floor(tmax / step + 1e-9)) * step)
}

# A group label longer than this is shown in the legend only, and the risk table
# marks its row with a short line in the group colour.
.cif_key_chars <- 14

# The three tabulated time points: 25%, 50% and 75% of the largest follow-up
# time, each rounded to a multiple of half the x step.
.cif_time_points <- function(tmax, step) {
  half <- step / 2
  floor(c(0.25, 0.5, 0.75) * tmax / half + 0.5) * half
}

# Codes as strings, as the Python engine reads them ("1", not "1.0").
.cif_codes <- function(v) {
  out <- trimws(as.character(v))
  out[is.na(v)] <- NA_character_
  out
}

# Group labels by the rule that both engines share (level_labels): a whole number
# reads "80", not "80.0".
.cif_group_labels <- function(v) {
  out <- trimws(level_labels(v))
  out[is.na(v)] <- NA_character_
  out
}

# The CIF of one curve at time t, as cmprsk::timepoints reads it: the value of
# the right-continuous step function, or NA past the last follow-up time of the
# group.
.cif_value_at <- function(t, time, est) {
  if (t > max(time)) return(NA_real_)
  est[max(which(time <= t))]
}

# Gray's test of one cause and the curves of each group, from cmprsk. The other
# causes are pooled into code 2, which leaves the CIF of the cause and Gray's
# test unchanged.
.cif_fit_cause <- function(cause, d) {
  st <- ifelse(d$status == "0", 0L, ifelse(d$status == cause, 1L, 2L))
  fit <- cmprsk::cuminc(d$time, st, d$group, cencode = 0)
  lv <- levels(d$group)
  curves <- do.call(rbind, lapply(seq_along(lv), .cif_curve_frame, fit = fit, lv = lv,
                                  cause = cause))
  test <- list(statistic = NA_real_, df = NA_real_, p_value = NA_real_)
  if (!is.null(fit$Tests)) {
    stat <- unname(fit$Tests["1", "stat"])
    df <- unname(fit$Tests["1", "df"])
    # cmprsk reports -1 when the variance matrix is singular.
    if (is.finite(stat) && stat >= 0) {
      test <- list(statistic = stat, df = df,
                   p_value = stats::pchisq(stat, df, lower.tail = FALSE))
    }
  }
  list(cause = cause, curves = curves, test = test)
}

# The step corners of one group, cause code 1, as cmprsk returns them.
.cif_curve_frame <- function(j, fit, lv, cause) {
  z <- fit[[j]]
  data.frame(cause = cause, group = lv[j], time = z$time, est = z$est,
             stringsAsFactors = FALSE)
}

.cif_at_risk <- function(g, d, breaks) {
  tg <- d$time[d$group == g]
  vapply(breaks, .cif_count_at_risk, numeric(1), tg = tg)
}

.cif_count_at_risk <- function(b, tg) sum(tg >= b)

.cif_stats_rows <- function(k, fits, d, tps) {
  f <- fits[[k]]
  lv <- levels(d$group)
  rows <- lapply(lv, .cif_stats_row, f = f, d = d, tps = tps)
  do.call(rbind, rows)
}

.cif_stats_row <- function(g, f, d, tps) {
  cv <- f$curves[f$curves$group == g, , drop = FALSE]
  vals <- vapply(tps, .cif_value_at, numeric(1), time = cv$time, est = cv$est)
  row <- data.frame(cause = f$cause, group = g, n = sum(d$group == g),
                    n_events = sum(d$group == g & d$status == f$cause),
                    stringsAsFactors = FALSE)
  for (i in seq_along(tps)) row[[sprintf("cif_%s", .cif_fmt_time(tps[i]))]] <- vals[i]
  row$gray_statistic <- f$test$statistic
  row$gray_df <- f$test$df
  row$gray_p <- f$test$p_value
  row
}

.cif_fmt_time <- function(t) sprintf("%g", t)

.cif_p_label <- function(test) {
  if (!is.finite(test$p_value)) return("One group; no test")
  if (test$p_value < 0.001) return("Gray's test P < 0.001")
  sprintf("Gray's test P = %.2g", test$p_value)
}

# The legend columns: as many as fit across the figure width, matched to the
# Python engine.
.cif_legend_ncol <- function(lv, width) {
  per <- 0.085 * max(nchar(lv)) + 0.7
  ncol <- min(length(lv), max(1, floor((width - 0.5) / per)))
  # Both engines fill the legend column by column, so the column count is the
  # one that the row count implies.
  as.integer(ceiling(length(lv) / ceiling(length(lv) / ncol)))
}

.cif_panel <- function(k, fits, d, brks, xmax, ytop, ystep, pal, spec, ecol, xlab,
                       legend_title, ncol_leg) {
  f <- fits[[k]]
  lv <- levels(d$group)
  cv <- f$curves
  cv$group <- factor(cv$group, levels = lv)
  xexp <- ggplot2::expansion(mult = c(0.02, 0.04))
  xstep <- if (length(brks) > 1) brks[2] - brks[1] else 1
  ybrks <- (0:floor(ytop / ystep + 1e-9)) * ystep
  p <- ggplot2::ggplot(cv, ggplot2::aes(x = time, y = est, colour = group)) +
    ggplot2::geom_path(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, lv), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = brks, labels = .cif_tick_labels(brks, xstep),
                                limits = c(0, xmax), expand = xexp) +
    ggplot2::scale_y_continuous(breaks = ybrks, labels = .cif_tick_labels(ybrks, ystep),
                                limits = c(0, ytop), expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::labs(x = xlab, y = spec$appearance$y_label %||% "Cumulative incidence",
                  colour = legend_title, title = sprintf("%s = %s", ecol, f$cause),
                  subtitle = .cif_p_label(f$test)) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size = 11),
                   legend.title = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30")) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = ncol_leg, title.position = "top",
                                                   title.hjust = 0.5))

  risk <- data.frame(group = factor(rep(lv, each = length(brks)), levels = rev(lv)),
                     time = rep(brks, length(lv)),
                     n = as.vector(vapply(lv, .cif_at_risk, numeric(length(brks)),
                                          d = d, breaks = brks)))
  keyed <- max(nchar(lv)) > .cif_key_chars
  tb <- ggplot2::ggplot(risk, ggplot2::aes(x = time, y = group, label = n, colour = group)) +
    ggplot2::geom_text(size = 3.2)
  if (keyed) {
    keys <- data.frame(group = factor(lv, levels = rev(lv)), x = -0.11 * xmax, xend = -0.05 * xmax)
    tb <- tb + ggplot2::geom_segment(data = keys, ggplot2::aes(x = x, xend = xend, y = group,
                                                               yend = group, colour = group),
                                     inherit.aes = FALSE, linewidth = 1.2)
  }
  tb <- tb +
    ggplot2::coord_cartesian(xlim = c(0, xmax), clip = "off") +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, lv)) +
    ggplot2::scale_x_continuous(breaks = brks, expand = xexp) +
    ggplot2::labs(x = NULL, y = NULL, title = "Number at risk") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = if (keyed) ggplot2::element_blank() else
                     ggplot2::element_text(size = 9, colour = "black",
                                           margin = ggplot2::margin(r = 12)),
                   plot.title = ggplot2::element_text(size = 10, hjust = 0),
                   plot.margin = ggplot2::margin(0, 14, 4, 10),
                   legend.position = "none")
  list(curve = p, table = tb)
}

.cif_column <- function(panel, ng) {
  ggpubr::ggarrange(panel$curve + ggplot2::theme(legend.position = "none"), panel$table,
                    ncol = 1, nrow = 2, heights = c(3.2, 0.55 + 0.24 * ng), align = "v")
}

recipe_cumulative_incidence <- function(df, spec) {
  tcol <- spec$data$time; ecol <- spec$data$event; gcol <- spec$data$x
  if (is.null(tcol) || is.null(ecol)) {
    stop("cumulative_incidence needs --time and --event (0 = censored, other values = the cause)")
  }
  has_group <- !is.null(gcol) && nzchar(gcol)
  for (cn in c(tcol, ecol, if (has_group) gcol)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(time = suppressWarnings(as.numeric(df[[tcol]])),
                  status = .cif_codes(df[[ecol]]),
                  group = if (has_group) .cif_group_labels(df[[gcol]]) else "all",
                  stringsAsFactors = FALSE)
  keep <- is.finite(d$time) & d$time >= 0 & !is.na(d$status) & !is.na(d$group)
  d <- d[keep, , drop = FALSE]
  if (nrow(d) < 2) stop("cumulative_incidence needs at least 2 subjects with a time and a status")
  causes <- sort(unique(d$status[d$status != "0"]), method = "radix")
  if (length(causes) == 0) stop(sprintf("column '%s' holds no event (every value is 0 = censored)", ecol))

  given <- spec$data$cause
  clean_steps <- character(0)
  if (!is.null(given) && nzchar(given)) {
    cause <- .cif_codes(given)
    if (!cause %in% causes) {
      stop(sprintf("--cause '%s' is not an event code of '%s' (%s)", cause, ecol,
                   paste(causes, collapse = ", ")))
    }
  } else if ("1" %in% causes) {
    cause <- "1"
  } else {
    cause <- causes[1]
    clean_steps <- c(clean_steps, sprintf(
      "No --cause was given and no event has code 1, so the cause of interest is the first code in sorted order, %s.",
      cause))
  }
  causes <- c(cause, setdiff(causes, cause))
  # The groups sort as numbers when every label is a number, else in code-point
  # order, by the level rule that both engines share (level_order).
  lv <- level_order(d$group)
  d$group <- factor(d$group, levels = lv)
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing or negative time, a missing status or a missing group.",
                                sum(!keep)),
    sprintf("Kept %d of %d row(s): %s.", nrow(d), n0,
            paste(sprintf("%d censored", sum(d$status == "0")),
                  paste(sprintf("%d with %s = %s", vapply(causes, .cif_count_code, numeric(1), d = d),
                                ecol, causes), collapse = ", "), sep = ", ")),
    clean_steps)

  fits <- lapply(causes, .cif_fit_cause, d = d)
  tmax <- max(d$time)
  ax <- .cif_breaks(tmax)
  tps <- .cif_time_points(tmax, ax$step)
  stats_df <- do.call(rbind, lapply(seq_along(fits), .cif_stats_rows, fits = fits, d = d, tps = tps))
  rownames(stats_df) <- NULL

  ymax <- max(vapply(fits, .cif_max_est, numeric(1)))
  ystep <- .cif_nice_step(max(ymax, 0.05), 4)
  ytop <- min(1, ceiling(max(ymax, 0.05) / ystep - 1e-9) * ystep)
  pal <- pub_palette(length(lv), spec$appearance$palette)
  xlab <- spec$appearance$x_label %||% tcol
  width <- .cif_width(length(causes), lv)
  ncol_leg <- .cif_legend_ncol(lv, width)
  panels <- lapply(seq_along(fits), .cif_panel, fits = fits, d = d, brks = ax$breaks,
                   xmax = tmax, ytop = ytop, ystep = ystep, pal = pal, spec = spec,
                   ecol = ecol, xlab = xlab, legend_title = if (has_group) gcol else NULL,
                   ncol_leg = ncol_leg)
  cols <- lapply(panels, .cif_column, ng = length(lv))
  body <- ggpubr::ggarrange(plotlist = cols, nrow = 1, ncol = length(cols))
  legend_rows <- if (has_group) ceiling(length(lv) / ncol_leg) else 0
  dims <- .cif_dims(width, lv, legend_rows)
  plot <- body
  if (has_group) {
    leg <- ggpubr::as_ggplot(ggpubr::get_legend(panels[[1]]$curve))
    plot <- ggpubr::ggarrange(leg, body, ncol = 1, nrow = 2,
                              heights = c(dims$legend_h, dims$height - dims$legend_h))
  }
  if (!is.null(spec$appearance$title)) {
    plot <- ggpubr::annotate_figure(plot, top = ggpubr::text_grob(spec$appearance$title,
                                                                  face = "bold", size = 14))
  }

  test_meta <- list(
    name = if (length(lv) > 1) "Cumulative incidence (Aalen-Johansen) with Gray's test"
           else "Cumulative incidence (Aalen-Johansen)",
    cause_of_interest = cause, causes = as.list(causes), groups = as.list(lv),
    time_points = as.list(tps),
    gray = stats::setNames(lapply(fits, .cif_meta_test), causes))
  resolved <- list(family = "competing_risks", label = test_meta$name)
  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .cif_methods(fits, d, causes, tps, tcol, ecol, gcol, has_group),
       build_script = list(fn = .cif_emit_script,
                           args = list(spec = spec, has_group = has_group, tps = tps,
                                       causes = causes, dims = dims, lv = lv)),
       label = slugify(ecol), df_used = d, clean_steps = clean_steps,
       curves = do.call(rbind, lapply(fits, .cif_get_curves)),
       width = dims$width, height = dims$height)
}

.cif_count_code <- function(code, d) sum(d$status == code)
.cif_max_est <- function(f) max(f$curves$est)
.cif_get_curves <- function(f) f$curves
.cif_meta_test <- function(f) f$test

# Figure width: one panel of 4.2 in per cause, wider when the risk table names
# the groups in full.
.cif_width <- function(n_causes, lv) {
  chars <- max(nchar(lv))
  extra <- if (chars > .cif_key_chars) 0 else max(0, 0.075 * chars - 0.45)
  n_causes * (4.2 + extra) + 0.2
}

# Figure height: taller with more groups and more legend rows.
.cif_dims <- function(width, lv, legend_rows) {
  legend_h <- if (legend_rows > 0) 0.45 + 0.27 * legend_rows else 0
  list(width = width,
       height = legend_h + 4.4 + 0.24 * length(lv),
       legend_h = legend_h)
}

.cif_methods <- function(fits, d, causes, tps, tcol, ecol, gcol, has_group) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions(c("cmprsk", "ggplot2"))
  grp_txt <- if (has_group) {
    sprintf(" in each group of %s (%s)", gcol, paste(sprintf("%s, n = %d", levels(d$group),
                                                              as.vector(table(d$group))),
                                                      collapse = "; "))
  } else ""
  tests <- if (has_group && length(levels(d$group)) > 1) {
    paste0(" The groups were compared for each event type with Gray's k-sample test (rho = 0): ",
           paste(vapply(fits, .cif_test_text, character(1), ecol = ecol), collapse = "; "), ".")
  } else " With one group, no test was applied."
  paste0(
    sprintf("The cumulative incidence of each event type (%s = %s; 0 = censored) was estimated with the Aalen-Johansen estimator%s, which treats the other event types as competing events rather than as censoring.",
            ecol, paste(causes, collapse = ", "), grp_txt),
    sprintf(" The event type of interest is %s = %s.", ecol, causes[1]),
    tests,
    sprintf(" The cumulative incidence at %s = %s is given in the stats table, and the number at risk (the subjects with a follow-up time of t or more) is shown below each panel.",
            tcol, .cif_join(vapply(tps, .cif_fmt_time, character(1)))),
    sprintf(" Computed in R %s with cmprsk %s and drawn with ggplot2 %s.", rver, v$cmprsk, v$ggplot2))
}

.cif_test_text <- function(f, ecol) {
  if (!is.finite(f$test$p_value)) return(sprintf("%s = %s, not estimable", ecol, f$cause))
  sprintf("%s = %s, chi-square(%d) = %.2f, P %s", ecol, f$cause, as.integer(f$test$df),
          f$test$statistic,
          if (f$test$p_value < 0.001) "< 0.001" else sprintf("= %.2g", f$test$p_value))
}

# Tick labels with the decimals of the step, so 0.0, 0.2, 0.4 and 0, 20, 40.
.cif_tick_labels <- function(values, step) {
  txt <- sprintf("%g", step)
  dec <- if (grepl(".", txt, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", txt)) else 0L
  formatC(values, format = "f", digits = dec)
}

.cif_join <- function(x) {
  if (length(x) <= 1) return(paste(x))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

.cif_emit_script <- function(spec, has_group, tps, causes, dims, lv, in_name, fig_stub) {
  tcol <- spec$data$time; ecol <- spec$data$event; gcol <- spec$data$x
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(cmprsk); library(ggplot2) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE, stringsAsFactors = FALSE)', in_name),
    sprintf('tm <- as.numeric(raw[["%s"]])', tcol),
    sprintf('status <- trimws(as.character(raw[["%s"]]))', ecol),
    if (has_group) sprintf('grp <- trimws(as.character(raw[["%s"]]))', gcol) else 'grp <- rep("all", nrow(raw))',
    "ok <- is.finite(tm) & tm >= 0 & !is.na(status) & !is.na(grp)",
    "tm <- tm[ok]; status <- status[ok]",
    sprintf("grp <- factor(grp[ok], levels = c(%s))", paste(sprintf("\"%s\"", lv), collapse = ", ")),
    'fit <- cuminc(tm, status, grp, cencode = "0")',
    "if (!is.null(fit$Tests)) print(fit$Tests)",
    sprintf("print(timepoints(fit, c(%s)))", paste(sprintf("%g", tps), collapse = ", ")),
    "curves <- NULL",
    "for (nm in setdiff(names(fit), \"Tests\")) {",
    "  curves <- rbind(curves, data.frame(curve = nm, time = fit[[nm]]$time, est = fit[[nm]]$est))",
    "}",
    "p <- ggplot(curves, aes(time, est, colour = curve)) + geom_path() +",
    sprintf('  labs(x = "%s", y = "Cumulative incidence", colour = NULL) + theme_classic(base_size = 13)',
            spec$appearance$x_label %||% tcol),
    sprintf('ggsave("%s.pdf", p, width = %.2f, height = %.2f)', fig_stub, dims$width, dims$height),
    'cat("Reproduced figure.\\n")'
  )
}
