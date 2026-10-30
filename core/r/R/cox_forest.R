# Recipe: Cox proportional-hazards model with a forest plot of hazard ratios.
# Tidy input: time, event (1 = event, 0 = censored), and covariate columns.
# The forest plot shows each covariate's hazard ratio with its 95% CI and
# p-value; the proportional-hazards assumption is checked and reported.
#
# The figure is drawn with ggplot2 in the house style and matches the Python
# engine row for row: the label column on the left, one row per model term
# (the first level of each factor gets its own "reference" row), squares with
# 95% CI whiskers on a log axis, and the "HR (95% CI)" and "P" columns on the
# right.

recipe_cox_forest <- function(df, spec) {
  time_col  <- spec$data$time
  event_col <- spec$data$event
  covs      <- spec$data$covariates
  if (is.null(time_col) || is.null(event_col)) {
    stop("cox_forest needs --time and --event columns")
  }
  if (is.null(covs) || length(covs) == 0) {
    covs <- setdiff(names(df), c(time_col, event_col))   # default: everything else
  }
  for (c in c(time_col, event_col, covs)) {
    if (!c %in% names(df)) stop(sprintf("column '%s' not found", c))
  }

  df[[time_col]] <- suppressWarnings(as.numeric(df[[time_col]]))
  df[[event_col]] <- event_codes(df[[event_col]], event_col)
  for (c in covs) if (is.character(df[[c]])) df[[c]] <- factor(df[[c]])

  n0   <- nrow(df)
  keep <- stats::complete.cases(df[, c(time_col, event_col, covs)]) &
          is.finite(df[[time_col]]) & df[[time_col]] >= 0
  df <- df[keep, , drop = FALSE]
  clean_steps <- if (nrow(df) < n0) {
    sprintf("Dropped %d row(s) with missing/invalid values (%d -> %d rows).",
            n0 - nrow(df), n0, nrow(df))
  } else "No cleaning applied; input used as-is."

  Surv <- survival::Surv
  rhs <- paste(sprintf("`%s`", covs), collapse = " + ")
  fml <- stats::as.formula(sprintf("Surv(`%s`, `%s`) ~ %s", time_col, event_col, rhs))
  cox <- survival::coxph(fml, data = df)

  zph <- tryCatch(survival::cox.zph(cox), error = .cox_null_on_error)
  ph_p <- if (!is.null(zph)) unname(zph$table["GLOBAL", "p"]) else NA_real_

  tid <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE)
  stats_df <- data.frame(
    term    = tid$term,
    HR      = tid$estimate,
    ci_low  = tid$conf.low,
    ci_high = tid$conf.high,
    z       = tid$statistic,
    p_value = tid$p.value,
    stringsAsFactors = FALSE
  )

  n_subj <- nrow(df)
  n_events <- sum(df[[event_col]])
  rows <- .cox_forest_rows(cox, stats_df, covs)
  fp <- .cox_forest_plot(rows, spec, n_subj, n_events, ph_p)

  test_meta <- list(
    name = "Cox proportional-hazards model",
    n = n_subj, events = n_events,
    covariates = as.list(covs),
    proportional_hazards = list(test = "cox.zph global", p_value = ph_p),
    hazard_ratios = lapply(seq_len(nrow(stats_df)), .cox_hr_entry, stats_df = stats_df)
  )

  resolved <- list(family = "cox", label = "Cox proportional hazards", covariates = covs)
  methods  <- .cox_methods(stats_df, ph_p, covs, n_subj, n_events)
  build_script <- list(fn = .cox_emit_script,
                       args = list(spec = spec, time_col = time_col, event_col = event_col,
                                   covs = covs, ticks = fp$ticks, label_in = fp$label_in,
                                   width = fp$width, height = fp$height))

  list(plot = fp$plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = "hazard_ratios",
       df_used = df, clean_steps = clean_steps,
       width = fp$width, height = fp$height)
}

.cox_null_on_error <- function(e) NULL

.cox_hr_entry <- function(i, stats_df) {
  list(term = stats_df$term[i], HR = stats_df$HR[i],
       ci_low = stats_df$ci_low[i], ci_high = stats_df$ci_high[i],
       p_value = stats_df$p_value[i])
}

# The significance stars of the forest plot: at most three, as
# survminer::ggforest graded them. The Python engine uses the same rule.
.cox_forest_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}

.cox_pfmt <- function(p) {
  if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else formatC(p, format = "g", digits = 2)
}

# The rows of one covariate. A factor gives a "reference" row for its first
# level and one row per other level; a numeric covariate gives one row. k is
# the row of the term in the stats table, NA for a reference row.
.cox_cov_rows <- function(i, covs, asg, xlev) {
  cv <- covs[i]
  idx <- asg[[i]]
  lev <- xlev[[names(asg)[i]]]
  if (is.null(lev)) {
    return(data.frame(label = cv, ref = FALSE, k = idx[1], stringsAsFactors = FALSE))
  }
  data.frame(label = paste0(cv, ": ", lev[seq_len(length(idx) + 1)]),
             ref = c(TRUE, rep(FALSE, length(idx))),
             k = c(NA_integer_, idx), stringsAsFactors = FALSE)
}

# One row per displayed line of the forest plot, in model order, with the
# hazard ratio, its CI and P copied from the stats table.
.cox_forest_rows <- function(cox, stats_df, covs) {
  asg <- cox$assign
  rows <- do.call(rbind, lapply(seq_along(asg), .cox_cov_rows, covs = covs, asg = asg,
                                xlev = cox$xlevels %||% list()))
  rows$HR      <- stats_df$HR[rows$k]
  rows$ci_low  <- stats_df$ci_low[rows$k]
  rows$ci_high <- stats_df$ci_high[rows$k]
  rows$p_value <- stats_df$p_value[rows$k]
  rows$hr_txt <- ifelse(rows$ref, "reference",
                        sprintf("%.2f (%.2f\u2013%.2f)", rows$HR, rows$ci_low, rows$ci_high))
  p_txt <- vapply(rows$p_value, .cox_pfmt, character(1))
  rows$p_txt <- ifelse(rows$ref, "", trimws(paste(p_txt, .cox_forest_stars(rows$p_value))))
  rows$y <- rev(seq_len(nrow(rows)))
  rows
}

# The values of one tick sequence that fall in [lo, hi], with the outward
# neighbour added on a side of 1 that the CIs reach but no tick covers.
.cox_ticks_for <- function(mant, lo, hi) {
  seq_all <- sort(as.vector(outer(mant, 10^(-4:4))))
  ticks <- seq_all[seq_all >= lo & seq_all <= hi]
  if (lo < 1 && !any(ticks < 1) && any(seq_all < lo)) ticks <- c(max(seq_all[seq_all < lo]), ticks)
  if (hi > 1 && !any(ticks > 1) && any(seq_all > hi)) ticks <- c(ticks, min(seq_all[seq_all > hi]))
  ticks
}

# Ticks of the hazard-ratio axis. The 1-2-5 sequence comes first; it is kept
# when it covers the CIs on both sides of 1 without widening the axis by more
# than 60%. A narrow range falls back to finer steps, and a range over many
# decades to the powers of ten. The Python engine uses the same rule.
.cox_log_ticks <- function(lo, hi) {
  lo <- min(lo, 1); hi <- max(hi, 1)
  span <- log(hi / lo)
  sets <- list(c(1, 2, 5), c(1, 1.5, 2, 3, 5, 7), c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8))
  if (length(.cox_ticks_for(c(1, 2, 5), lo, hi)) > 8) sets <- list(1)
  ticks <- 1
  for (mant in sets) {
    ticks <- .cox_ticks_for(mant, lo, hi)
    total <- log(max(hi, ticks) / min(lo, ticks))
    if (total <= 1.6 * span) return(ticks)
  }
  ticks
}

.cox_tick_labels <- function(ticks) trimws(formatC(ticks, format = "fg", digits = 3))

# The width of the label column in inches, from the longest label.
.cox_label_inches <- function(labels) 0.3 + 0.085 * max(nchar(labels), 4)

.cox_forest_plot <- function(rows, spec, n, events, ph_p) {
  ap <- spec$appearance
  k <- nrow(rows)
  est <- rows[!rows$ref & is.finite(rows$HR), , drop = FALSE]
  fin_lo <- est$ci_low[is.finite(est$ci_low) & est$ci_low > 0]
  fin_hi <- est$ci_high[is.finite(est$ci_high)]
  lo <- if (length(fin_lo)) min(fin_lo) else 1
  hi <- if (length(fin_hi)) max(fin_hi) else 1
  ticks <- .cox_log_ticks(lo, hi)
  xlim <- c(min(lo, ticks), max(hi, ticks))
  ylim <- c(0.4, k + 1.6)
  pal <- ap$palette %||% pubready_palette
  label_in <- .cox_label_inches(rows$label)
  width <- label_in + 3.0 + 2.5
  height <- max(3.0, 0.3 * (k + 1) + 1.6 + if (!is.null(ap$title)) 0.3 else 0)

  forest <- ggplot2::ggplot(est) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(y = .data$y, xmin = .data$ci_low, xmax = .data$ci_high),
                            height = 0, linewidth = 0.5, colour = "grey15", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(x = .data$HR, y = .data$y), shape = 15, size = 3,
                        colour = pal[1], na.rm = TRUE) +
    ggplot2::scale_x_log10(breaks = ticks, labels = .cox_tick_labels(ticks), limits = xlim,
                           oob = scales::oob_keep) +
    ggplot2::scale_y_continuous(breaks = rows$y, labels = rows$label, limits = ylim,
                                expand = ggplot2::expansion(0, 0)) +
    ggplot2::labs(x = ap$x_label %||% "Hazard ratio (95% CI, log scale)", y = NULL,
                  title = ap$title, subtitle = .cox_subtitle(n, events, ph_p)) +
    pub_theme(ap$theme) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(hjust = 0, size = 10.5),
                   plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey30"),
                   plot.margin = ggplot2::margin(10, 4, 10, 10))

  table <- ggplot2::ggplot(rows, ggplot2::aes(y = .data$y)) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = .data$hr_txt), hjust = 0, size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(x = 0.66, label = .data$p_txt), hjust = 0, size = 3.5) +
    ggplot2::annotate("text", x = c(0, 0.66), y = k + 1, label = c("HR (95% CI)", "P"),
                      hjust = 0, fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(0, 0)) +
    ggplot2::labs(title = if (!is.null(ap$title)) " ", subtitle = " ", x = " ") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 13),
                   plot.subtitle = ggplot2::element_text(size = 9.5),
                   axis.title.x = ggplot2::element_text(size = 13, colour = "white"),
                   axis.text.x = ggplot2::element_text(size = 10.4, colour = "white"),
                   plot.margin = ggplot2::margin(10, 10, 10, 6))

  plot <- ggpubr::ggarrange(forest, table, ncol = 2, align = "h",
                            widths = c(label_in + 3.0, 2.5))
  list(plot = plot, ticks = ticks, label_in = label_in, width = width, height = height)
}

.cox_methods <- function(stats_df, ph_p, covs, n, events) {
  v    <- pkg_versions(c("survival", "ggplot2", "broom"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  sig <- stats_df[which(stats_df$p_value < 0.05), , drop = FALSE]
  sig_txt <- if (nrow(sig) > 0) {
    paste0(" Significant predictors: ",
           paste(sprintf("%s (HR = %.2f, 95%% CI %.2f-%.2f, %s)",
                         sig$term, sig$HR, sig$ci_low, sig$ci_high,
                         vapply(sig$p_value, p_phrase, character(1))), collapse = "; "),
           ".")
  } else " No individual predictor reached significance."
  ph_txt <- if (is.na(ph_p)) "" else {
    sprintf(" The proportional-hazards assumption was %s (global cox.zph %s).",
            if (ph_p < 0.05) "violated" else "not violated", p_phrase(ph_p))
  }

  paste0(
    sprintf("Hazard ratios were estimated with a Cox proportional-hazards model (%d subjects, %d events) including %s.",
            n, events, paste(covs, collapse = ", ")),
    sig_txt, ph_txt,
    sprintf(" Analyses were performed in R %s with survival %s, and the forest plot was drawn with ggplot2 %s (tidied with broom %s).",
            rver, v$survival, v$ggplot2, v$broom)
  )
}

.cox_emit_script <- function(spec, time_col, event_col, covs, ticks, label_in, width, height,
                             in_name, fig_stub) {
  rhs <- paste(sprintf("`%s`", covs), collapse = " + ")
  pal <- (spec$appearance$palette %||% pubready_palette)[1]
  xlab <- spec$appearance$x_label %||% "Hazard ratio (95% CI, log scale)"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(survival); library(ggplot2); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf("covs <- c(%s)", paste(sprintf('"%s"', covs), collapse = ", ")),
    "for (c in covs) if (is.character(df[[c]])) df[[c]] <- factor(df[[c]])",
    "",
    sprintf("cox <- coxph(Surv(`%s`, `%s`) ~ %s, data = df)", time_col, event_col, rhs),
    'ci <- exp(confint(cox)); pv <- summary(cox)$coefficients[, "Pr(>|z|)"]',
    "",
    "# One row per model term; the first level of each factor is the reference.",
    "rows <- NULL",
    "for (i in seq_along(cox$assign)) {",
    "  idx <- cox$assign[[i]]; lev <- cox$xlevels[[names(cox$assign)[i]]]",
    "  if (is.null(lev)) {",
    "    rows <- rbind(rows, data.frame(label = covs[i], k = idx[1]))",
    "  } else {",
    '    rows <- rbind(rows, data.frame(label = paste0(covs[i], ": ", lev[seq_len(length(idx) + 1)]),',
    "                                   k = c(NA, idx)))",
    "  }",
    "}",
    "rows$HR <- exp(coef(cox))[rows$k]; rows$lo <- ci[rows$k, 1]; rows$hi <- ci[rows$k, 2]",
    "rows$p <- pv[rows$k]; rows$y <- rev(seq_len(nrow(rows)))",
    'stars <- ifelse(is.na(rows$p), "", ifelse(rows$p < 0.001, "***",',
    '                ifelse(rows$p < 0.01, "**", ifelse(rows$p < 0.05, "*", ""))))',
    'p_txt <- ifelse(rows$p < 0.001, "< 0.001", formatC(rows$p, format = "g", digits = 2))',
    'rows$hr_txt <- ifelse(is.na(rows$k), "reference",',
    '                      sprintf("%.2f (%.2f\\u2013%.2f)", rows$HR, rows$lo, rows$hi))',
    'rows$p_txt <- ifelse(is.na(rows$k), "", trimws(paste(p_txt, stars)))',
    "est <- rows[!is.na(rows$k), ]",
    sprintf("ticks <- c(%s)", paste(format(ticks, digits = 15, trim = TRUE), collapse = ", ")),
    "ylim <- c(0.4, nrow(rows) + 1.6)",
    "",
    "forest <- ggplot(est) +",
    '  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40") +',
    '  geom_errorbarh(aes(y = y, xmin = lo, xmax = hi), height = 0, colour = "grey15") +',
    sprintf('  geom_point(aes(x = HR, y = y), shape = 15, size = 3, colour = "%s") +', pal),
    '  scale_x_log10(breaks = ticks, labels = trimws(formatC(ticks, format = "fg", digits = 3)),',
    "                limits = range(c(est$lo, est$hi, ticks), finite = TRUE)) +",
    "  scale_y_continuous(breaks = rows$y, labels = rows$label, limits = ylim, expand = c(0, 0)) +",
    sprintf("  labs(x = %s, y = NULL,", deparse(xlab)),
    '       subtitle = sprintf("%d subjects, %d events", cox$n, cox$nevent)) +',
    "  theme_classic(base_size = 13) +",
    "  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),",
    '        axis.text = element_text(colour = "black"),',
    "        axis.text.y = element_text(hjust = 0, size = 10.5))",
    "table <- ggplot(rows, aes(y = y)) +",
    "  geom_text(aes(x = 0, label = hr_txt), hjust = 0, size = 3.5) +",
    "  geom_text(aes(x = 0.66, label = p_txt), hjust = 0, size = 3.5) +",
    '  annotate("text", x = c(0, 0.66), y = nrow(rows) + 1, label = c("HR (95% CI)", "P"),',
    '           hjust = 0, fontface = "bold", size = 3.5) +',
    "  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +",
    "  scale_y_continuous(limits = ylim, expand = c(0, 0)) +",
    '  labs(subtitle = " ", x = " ") + theme_void() +',
    '  theme(axis.title.x = element_text(size = 13, colour = "white"),',
    '        axis.text.x = element_text(size = 10.4, colour = "white"))',
    sprintf('p <- ggpubr::ggarrange(forest, table, ncol = 2, align = "h", widths = c(%.2f, 2.5))',
            label_in + 3.0),
    "",
    sprintf('ggsave("%s.pdf", p, width = %.2f, height = %.2f, device = cairo_pdf)',
            fig_stub, width, height),
    'cat("Reproduced figure.\\n")'
  )
}

# The subtitle gives the counts and, when the proportional-hazards check fails,
# a warning, so the figure does not show hazard ratios of a broken model as valid.
.cox_subtitle <- function(n, events, ph_p) {
  txt <- sprintf("%d subjects, %d events", n, events)
  if (!is.na(ph_p) && ph_p < 0.05) txt <- sprintf("%s; PH assumption violated (%s)", txt, p_phrase(ph_p))
  txt
}
