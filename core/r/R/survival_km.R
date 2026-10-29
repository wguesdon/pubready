# Recipe: Kaplan-Meier survival curves compared between groups.
# Event column is coded 1 = event, 0 = censored. Draws step curves with
# censoring ticks, the log-rank p-value on the plot, and a number-at-risk table
# below. Tidy input: time, event, and a grouping column (--x).

recipe_survival_km <- function(df, spec) {
  time_col  <- spec$data$time
  event_col <- spec$data$event
  group_col <- spec$data$x
  if (is.null(time_col) || is.null(event_col)) {
    stop("survival_km needs --time and --event columns")
  }
  if (is.null(group_col) || !nzchar(group_col)) {
    stop("survival_km needs a grouping column via --x")
  }
  for (c in c(time_col, event_col, group_col)) {
    if (!c %in% names(df)) stop(sprintf("column '%s' not found", c))
  }

  df[[time_col]] <- suppressWarnings(as.numeric(df[[time_col]]))
  ev <- suppressWarnings(as.integer(df[[event_col]]))
  if (!all(ev %in% c(0L, 1L) | is.na(ev))) {
    stop(sprintf("event column '%s' must be coded 1 = event, 0 = censored", event_col))
  }
  df[[event_col]] <- ev

  n0   <- nrow(df)
  keep <- stats::complete.cases(df[, c(time_col, event_col, group_col)]) &
          is.finite(df[[time_col]]) & df[[time_col]] >= 0
  df <- df[keep, , drop = FALSE]
  clean_steps <- if (nrow(df) < n0) {
    sprintf("Dropped %d row(s) with missing/invalid time, event, or group (%d -> %d rows).",
            n0 - nrow(df), n0, nrow(df))
  } else "No cleaning applied; input used as-is."
  df[[group_col]] <- factor(df[[group_col]])
  lvls <- levels(df[[group_col]])

  # Bind Surv into the formula's environment so surv_fit/survdiff resolve it;
  # surv_fit (not survfit) stores the formula so ggsurvplot can recover the data.
  Surv <- survival::Surv
  fml <- stats::as.formula(sprintf("Surv(`%s`, `%s`) ~ `%s`",
                                   time_col, event_col, group_col))
  fit <- survminer::surv_fit(fml, data = df)
  sd  <- survival::survdiff(fml, data = df)
  lr_df <- length(sd$n) - 1
  lr_p  <- stats::pchisq(sd$chisq, lr_df, lower.tail = FALSE)

  pal  <- spec$appearance$palette %||% pubready_palette
  pal  <- rep(pal, length.out = length(lvls))
  xlab <- spec$appearance$x_label %||% "Time"
  ylab <- spec$appearance$y_label %||% "Survival probability"

  # The x axis runs to the last time, so the last event is on the figure, with
  # breaks of the 1-2-5 sequence that both engines compute in the same way. The
  # number at risk at each break is computed here (time >= t), not by the
  # plotting package, so the two engines print the same table.
  tmax  <- max(df[[time_col]])
  step  <- km_break_step(tmax)
  brks  <- km_breaks(tmax, step)
  risk  <- km_at_risk(df[[time_col]], df[[group_col]], lvls, brks)
  leg_labs <- wrap_labels(lvls, width = .km_legend_wrap())

  g <- survminer::ggsurvplot(
    fit, data = df, pval = FALSE, risk.table = FALSE, conf.int = FALSE,
    censor = TRUE, palette = pal[seq_along(lvls)],
    legend = "right", legend.title = group_col, legend.labs = leg_labs,
    xlab = xlab, ylab = ylab, xlim = c(0, tmax), break.time.by = step,
    ggtheme = pub_theme(spec$appearance$theme)
  )
  # The log-rank P goes in a subtitle, where no curve can cross it. ggsurvplot
  # draws it at a fixed spot inside the panel, on top of curves that drop early.
  # The legend sits right of the panel, where it never covers a curve.
  curve <- g$plot + ggplot2::labs(subtitle = .km_p_label(lr_p)) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"),
                   legend.position = "right", legend.justification = "top",
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 9),
                   legend.key.height = ggplot2::unit(0.9, "lines"))
  table <- .km_risk_plot(risk, lvls, pal[seq_along(lvls)], brks, tmax)
  plot <- ggpubr::ggarrange(curve, table, ncol = 1, nrow = 2,
                            heights = c(1 - km_table_share(length(lvls)), km_table_share(length(lvls))),
                            align = "v")

  # Per-group summary: n, events, median survival with 95% CI.
  tb <- summary(fit)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  stats_df <- data.frame(
    group           = lvls,
    n               = as.integer(.km_col(tb, "records")),
    events          = as.integer(.km_col(tb, "events")),
    median_survival = .km_col(tb, "median"),
    ci_low          = .km_col(tb, "0.95LCL"),
    ci_high         = .km_col(tb, "0.95UCL"),
    log_rank_p      = lr_p,
    stringsAsFactors = FALSE
  )

  test_meta <- list(
    name = "Kaplan-Meier estimate with log-rank test",
    log_rank = list(chisq = unname(sd$chisq), df = lr_df, p_value = lr_p),
    groups = lapply(seq_len(nrow(stats_df)), .km_group_meta, stats_df = stats_df),
    at_risk = list(times = brks,
                   counts = stats::setNames(lapply(seq_along(lvls), .km_risk_row, risk = risk),
                                            lvls))
  )

  resolved <- list(family = "survival", label = "Kaplan-Meier / log-rank",
                   group = group_col)
  methods <- .km_methods(stats_df, lr_p, sd, lr_df, spec, group_col)
  width <- km_width(leg_labs, group_col)
  build_script <- list(fn = .km_emit_script,
                       args = list(spec = spec, time_col = time_col, event_col = event_col,
                                   group_col = group_col, pal = pal[seq_along(lvls)],
                                   lvls = lvls, xlab = xlab, ylab = ylab, step = step,
                                   width = width))

  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = group_col,
       df_used = df, clean_steps = clean_steps,
       width = width, height = 5.8)
}

.km_col <- function(tb, nm) if (nm %in% colnames(tb)) tb[, nm] else rep(NA_real_, nrow(tb))

.km_group_meta <- function(i, stats_df) {
  list(group = stats_df$group[i], n = stats_df$n[i], events = stats_df$events[i],
       median_survival = stats_df$median_survival[i])
}

.km_risk_row <- function(i, risk) as.integer(risk[i, ])

# The legend labels wrap at this many characters.
.km_legend_wrap <- function() 24

# The step between the x breaks: the smallest value of the 1-2-5 sequence that
# cuts the follow-up into at most 6 intervals. Python has the same rule.
km_break_step <- function(tmax) {
  if (!is.finite(tmax) || tmax <= 0) return(1)
  cand <- sort(as.vector(outer(c(1, 2, 5), 10^(-3:6))))
  cand[which(tmax / cand <= 6 + 1e-9)[1]]
}

# The breaks from 0 to the last time, by `step`.
km_breaks <- function(tmax, step) {
  round(step * (0:floor(tmax / step + 1e-9)), 10)
}

# Subjects at risk at each break t: those whose time is t or later, one row per
# group. This is the definition of survival::summary.survfit(times = t).
km_at_risk <- function(time, group, lvls, breaks) {
  grp <- factor(as.character(group), levels = lvls)
  out <- vapply(breaks, .km_count_at, numeric(length(lvls)), time = time, grp = grp)
  matrix(out, nrow = length(lvls), dimnames = list(lvls, NULL))
}

.km_count_at <- function(b, time, grp) {
  n <- tapply(time >= b, grp, sum)
  n[is.na(n)] <- 0
  as.numeric(n)
}

# The share of the figure height that the risk table takes: more rows, more room.
# Python uses the same formula.
km_table_share <- function(k) 0.08 + 0.045 * k

# The figure width: 4.6 in for the panel and the axis, plus the legend, which
# grows with its longest line. Python uses the same formula.
km_width <- function(leg_labs, title) {
  chars <- max(vapply(leg_labs, longest_line, integer(1)), nchar(title))
  4.6 + 0.55 + 0.075 * chars
}

# The number-at-risk table, drawn on the x scale of the curves. A short line in
# the colour of each group starts its row; the legend names the groups.
.km_risk_plot <- function(risk, lvls, pal, brks, tmax) {
  k <- length(lvls)
  tab <- data.frame(time = rep(brks, each = k), row = rep(rev(seq_len(k)), times = length(brks)),
                    n = as.vector(risk), stringsAsFactors = FALSE)
  keys <- data.frame(row = rev(seq_len(k)), group = factor(lvls, levels = lvls))
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$time, y = .data$row)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), size = 3.2) +
    ggplot2::geom_segment(data = keys,
                          ggplot2::aes(x = -0.13 * tmax, xend = -0.08 * tmax,
                                       y = .data$row, yend = .data$row, colour = .data$group),
                          linewidth = 1.6) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, lvls), guide = "none") +
    ggplot2::scale_x_continuous(breaks = brks) +
    ggplot2::scale_y_continuous(limits = c(0.5, k + 0.5)) +
    ggplot2::coord_cartesian(xlim = c(0, tmax), clip = "off") +
    ggplot2::labs(title = "Number at risk") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0,
                                                      margin = ggplot2::margin(b = 4)),
                   plot.margin = ggplot2::margin(4, 14, 6, 10))
}

.km_p_label <- function(p) {
  if (p < 0.001) "Log-rank P < 0.001" else sprintf("Log-rank P = %.2g", p)
}

.km_median_text <- function(m) {
  if (is.na(m) || !is.finite(m)) "not reached" else trimws(formatC(m, format = "g", digits = 3))
}

.km_methods <- function(stats_df, lr_p, sd, lr_df, spec, group_col) {
  v    <- pkg_versions(c("survival", "survminer"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  med_str <- vapply(stats_df$median_survival, .km_median_text, character(1))
  med <- paste(sprintf("%s (median %s)", stats_df$group, med_str), collapse = ", ")
  paste0(
    "Survival was estimated by the Kaplan-Meier method and compared between ",
    sprintf("groups of %s with the log-rank test ", group_col),
    sprintf("(chi-square(%d) = %.2f, %s).", lr_df, sd$chisq, p_phrase(lr_p)),
    sprintf(" Median survival by group: %s.", med),
    " Censoring is marked with ticks and the number at risk is tabulated below the curves;",
    " it counts the subjects whose time is at or after each time point.",
    sprintf(" Analyses were performed in R %s with survival %s and survminer %s.",
            rver, v$survival, v$survminer)
  )
}

.km_emit_script <- function(spec, time_col, event_col, group_col, pal, lvls, xlab, ylab,
                            step, width, in_name, fig_stub) {
  pal_str  <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  lvls_str <- paste0("c(", paste(sprintf('"%s"', lvls), collapse = ", "), ")")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(survival); library(survminer); library(ggpubr); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf('df[["%s"]] <- factor(df[["%s"]])', group_col, group_col),
    sprintf('fit <- surv_fit(Surv(`%s`, `%s`) ~ `%s`, data = df)', time_col, event_col, group_col),
    "",
    sprintf('sd <- survdiff(Surv(`%s`, `%s`) ~ `%s`, data = df)', time_col, event_col, group_col),
    "p_lr <- pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)",
    sprintf('tmax <- max(df[["%s"]])', time_col),
    "# The risk table of survminer counts the subjects with a time at or after each break.",
    'g <- ggsurvplot(fit, data = df, pval = FALSE, risk.table = TRUE, conf.int = FALSE,',
    sprintf('                censor = TRUE, palette = %s, legend = "right", legend.title = "%s",',
            pal_str, group_col),
    sprintf('                legend.labs = %s, xlab = "%s", ylab = "%s",', lvls_str, xlab, ylab),
    sprintf('                xlim = c(0, tmax), break.time.by = %s, risk.table.y.text = FALSE,', format(step)),
    '                ggtheme = theme_classic(base_size = 13),',
    '                tables.theme = theme_cleantable())',
    'g$plot <- g$plot + labs(subtitle = if (p_lr < 0.001) "Log-rank P < 0.001" else sprintf("Log-rank P = %.2g", p_lr))',
    'p <- ggarrange(g$plot, g$table, ncol = 1, nrow = 2, heights = c(0.74, 0.26), align = "v")',
    "",
    sprintf('ggsave("%s.pdf", p, width = %.2f, height = 5.8)', fig_stub, width),
    'cat("Reproduced figure.\\n")'
  )
}
