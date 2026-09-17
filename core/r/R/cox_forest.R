# Recipe: Cox proportional-hazards model with a forest plot of hazard ratios.
# Tidy input: time, event (1 = event, 0 = censored), and covariate columns.
# The forest plot shows each covariate's hazard ratio with its 95% CI and
# p-value; the proportional-hazards assumption is checked and reported.

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
  ev <- suppressWarnings(as.integer(df[[event_col]]))
  if (!all(ev %in% c(0L, 1L) | is.na(ev))) {
    stop(sprintf("event column '%s' must be coded 1 = event, 0 = censored", event_col))
  }
  df[[event_col]] <- ev
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

  zph <- tryCatch(survival::cox.zph(cox), error = function(e) NULL)
  ph_p <- if (!is.null(zph)) unname(zph$table["GLOBAL", "p"]) else NA_real_

  plot <- survminer::ggforest(cox, data = df)

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

  test_meta <- list(
    name = "Cox proportional-hazards model",
    n = nrow(df), events = sum(df[[event_col]]),
    covariates = as.list(covs),
    proportional_hazards = list(test = "cox.zph global", p_value = ph_p),
    hazard_ratios = lapply(seq_len(nrow(stats_df)), function(i) list(
      term = stats_df$term[i], HR = stats_df$HR[i],
      ci_low = stats_df$ci_low[i], ci_high = stats_df$ci_high[i],
      p_value = stats_df$p_value[i]
    ))
  )

  resolved <- list(family = "cox", label = "Cox proportional hazards", covariates = covs)
  methods  <- .cox_methods(stats_df, ph_p, covs, nrow(df), sum(df[[event_col]]))
  build_script <- function(in_name, fig_stub) {
    .cox_emit_script(spec, time_col, event_col, covs, in_name, fig_stub)
  }

  n_display <- length(covs) + nrow(stats_df) + 3
  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script, label = "hazard_ratios",
       df_used = df, clean_steps = clean_steps,
       width = 7.4, height = max(3.8, 0.42 * n_display + 1.2))
}

.cox_methods <- function(stats_df, ph_p, covs, n, events) {
  v    <- pkg_versions(c("survival", "survminer", "broom"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  pfmt <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else formatC(p, format = "g", digits = 2)

  sig <- stats_df[which(stats_df$p_value < 0.05), , drop = FALSE]
  sig_txt <- if (nrow(sig) > 0) {
    paste0(" Significant predictors: ",
           paste(sprintf("%s (HR = %.2f, 95%% CI %.2f-%.2f, P = %s)",
                         sig$term, sig$HR, sig$ci_low, sig$ci_high,
                         vapply(sig$p_value, pfmt, character(1))), collapse = "; "),
           ".")
  } else " No individual predictor reached significance."
  ph_txt <- if (is.na(ph_p)) "" else {
    sprintf(" The proportional-hazards assumption was %s (global cox.zph P = %s).",
            if (ph_p < 0.05) "violated" else "not violated", pfmt(ph_p))
  }

  paste0(
    sprintf("Hazard ratios were estimated with a Cox proportional-hazards model (%d subjects, %d events) including %s.",
            n, events, paste(covs, collapse = ", ")),
    sig_txt, ph_txt,
    sprintf(" Analyses were performed in R %s with survival %s, and the forest plot was drawn with survminer %s (tidied with broom %s).",
            rver, v$survival, v$survminer, v$broom)
  )
}

.cox_emit_script <- function(spec, time_col, event_col, covs, in_name, fig_stub) {
  rhs <- paste(sprintf("`%s`", covs), collapse = " + ")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(survival); library(survminer); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    paste0(sprintf('for (c in c(%s)) if (is.character(df[[c]])) df[[c]] <- factor(df[[c]])',
                   paste(sprintf('"%s"', covs), collapse = ", "))),
    "",
    sprintf('cox <- coxph(Surv(`%s`, `%s`) ~ %s, data = df)', time_col, event_col, rhs),
    "p <- ggforest(cox, data = df)",
    "",
    sprintf('ggsave("%s.pdf", p, width = 7.4, height = 5.0)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
