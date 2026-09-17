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

  g <- survminer::ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    censor = TRUE, palette = pal[seq_along(lvls)],
    legend.title = group_col, legend.labs = lvls,
    xlab = xlab, ylab = ylab,
    ggtheme = pub_theme(spec$appearance$theme),
    risk.table.fontsize = 3.2,
    tables.theme = survminer::theme_cleantable()
  )
  plot <- ggpubr::ggarrange(g$plot, g$table, ncol = 1, nrow = 2,
                            heights = c(0.74, 0.26), align = "v")

  # Per-group summary: n, events, median survival with 95% CI.
  tb <- summary(fit)$table
  if (is.null(dim(tb))) tb <- t(as.matrix(tb))
  col <- function(nm) if (nm %in% colnames(tb)) tb[, nm] else rep(NA_real_, nrow(tb))
  stats_df <- data.frame(
    group           = lvls,
    n               = as.integer(col("records")),
    events          = as.integer(col("events")),
    median_survival = col("median"),
    ci_low          = col("0.95LCL"),
    ci_high         = col("0.95UCL"),
    log_rank_p      = lr_p,
    stringsAsFactors = FALSE
  )

  test_meta <- list(
    name = "Kaplan-Meier estimate with log-rank test",
    log_rank = list(chisq = unname(sd$chisq), df = lr_df, p_value = lr_p),
    groups = lapply(seq_len(nrow(stats_df)), function(i) list(
      group = stats_df$group[i], n = stats_df$n[i], events = stats_df$events[i],
      median_survival = stats_df$median_survival[i]
    ))
  )

  resolved <- list(family = "survival", label = "Kaplan-Meier / log-rank",
                   group = group_col)
  methods <- .km_methods(stats_df, lr_p, sd, lr_df, spec, group_col)
  build_script <- function(in_name, fig_stub) {
    .km_emit_script(spec, time_col, event_col, group_col, pal[seq_along(lvls)],
                    lvls, xlab, ylab, in_name, fig_stub)
  }

  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = group_col,
       df_used = df, clean_steps = clean_steps,
       width = 5.6, height = 5.8)
}

.km_methods <- function(stats_df, lr_p, sd, lr_df, spec, group_col) {
  v    <- pkg_versions(c("survival", "survminer"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  pfmt <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else formatC(p, format = "g", digits = 2)
  med_str <- vapply(stats_df$median_survival, function(m) {
    if (is.na(m)) "not reached" else formatC(m, format = "g", digits = 3)
  }, character(1))
  med <- paste(sprintf("%s (median %s)", stats_df$group, med_str), collapse = ", ")
  paste0(
    "Survival was estimated by the Kaplan-Meier method and compared between ",
    sprintf("groups of %s with the log-rank test ", group_col),
    sprintf("(chi-square(%d) = %.2f, P = %s).", lr_df, sd$chisq, pfmt(lr_p)),
    sprintf(" Median survival by group: %s.", med),
    " Censoring is marked with ticks and the number at risk is tabulated below the curves.",
    sprintf(" Analyses were performed in R %s with survival %s and survminer %s.",
            rver, v$survival, v$survminer)
  )
}

.km_emit_script <- function(spec, time_col, event_col, group_col, pal, lvls, xlab, ylab, in_name, fig_stub) {
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
    sprintf('g <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,'),
    sprintf('                censor = TRUE, palette = %s, legend.title = "%s",', pal_str, group_col),
    sprintf('                legend.labs = %s, xlab = "%s", ylab = "%s",', lvls_str, xlab, ylab),
    '                ggtheme = theme_classic(base_size = 13),',
    '                tables.theme = theme_cleantable())',
    'p <- ggarrange(g$plot, g$table, ncol = 1, nrow = 2, heights = c(0.74, 0.26), align = "v")',
    "",
    sprintf('ggsave("%s.pdf", p, width = 5.6, height = 5.8)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
