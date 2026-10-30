# Recipe: a repeated-measures time course. Each subject (--id) is measured at
# every time point (--x) and may belong to a group (--group). The figure shows
# the group mean +/- SEM at each time with thin lines for the subjects. The test
# is a mixed ANOVA (group between subjects, time within subjects), or a one-way
# repeated-measures ANOVA over time when there is no group.
#
# The sums of squares are the classical split-plot decomposition, which equals
# stats::aov with an Error(id/time) stratum. The Greenhouse-Geisser epsilon and
# Mauchly's test come from the pooled within-group covariance of the time points,
# as car::Anova computes them. With --test auto the within-subject effects take
# the Greenhouse-Geisser correction when Mauchly's test rejects sphericity
# (P < 0.05); --test gg always applies it and --test uncorrected never does.
# The ANOVA needs every subject at every time, so a subject with a missing visit
# is dropped and the drop is recorded.

.tc_p_label <- function(p) if (p < 0.001) "P < 0.001" else paste("P =", sig_text(p, 2))

.tc_group_cells <- function(gl, Y, g) colMeans(Y[g == gl, , drop = FALSE])

.tc_centred_ssp <- function(gl, Y, g) {
  Yg <- Y[g == gl, , drop = FALSE]
  Z <- sweep(Yg, 2, colMeans(Yg))
  crossprod(Z)
}

# Orthonormal contrasts among the time points. Epsilon and Mauchly's W do not
# depend on which orthonormal basis is used.
.tc_contrasts <- function(n_time) stats::contr.poly(n_time)

# The mixed ANOVA on a subjects-by-times matrix Y with the group of each row.
.tc_mixed_anova <- function(Y, g) {
  lv <- levels(g); G <- length(lv); N <- nrow(Y); Tn <- ncol(Y)
  grand <- mean(Y)
  n_g <- as.numeric(table(g)[lv])
  cells <- do.call(rbind, lapply(lv, .tc_group_cells, Y = Y, g = g))
  grp_mean <- rowMeans(cells)
  ss_subj <- Tn * sum((rowMeans(Y) - grand)^2)
  ss_group <- Tn * sum(n_g * (grp_mean - grand)^2)
  ss_time <- N * sum((colMeans(Y) - grand)^2)
  ss_cells <- sum(n_g * rowSums((cells - grand)^2))
  ss_inter <- ss_cells - ss_group - ss_time
  ss_err_b <- ss_subj - ss_group
  ss_err_w <- sum((Y - grand)^2) - ss_subj - ss_time - ss_inter

  df_err_b <- N - G; df_err_w <- (N - G) * (Tn - 1)
  ssp <- Reduce(`+`, lapply(lv, .tc_centred_ssp, Y = Y, g = g))
  C <- .tc_contrasts(Tn)
  M <- t(C) %*% ssp %*% C
  p <- Tn - 1
  eps <- if (p == 1) 1 else sum(diag(M))^2 / (p * sum(diag(M %*% M)))
  mauchly <- .tc_mauchly(M, p, df_err_b)

  rows <- list(
    list(term = "time", ss = ss_time, df1 = Tn - 1, ss_err = ss_err_w, df2 = df_err_w, within = TRUE),
    list(term = "interaction", ss = ss_inter, df1 = (G - 1) * (Tn - 1), ss_err = ss_err_w, df2 = df_err_w, within = TRUE))
  if (G >= 2) {
    rows <- c(list(list(term = "group", ss = ss_group, df1 = G - 1, ss_err = ss_err_b, df2 = df_err_b,
                        within = FALSE)), rows)
  } else {
    rows <- rows[1]
  }
  list(rows = rows, eps = eps, mauchly = mauchly, cells = cells, grand = grand)
}

# Mauchly's sphericity test with the second-order correction of stats::mauchly.test
# and car::Anova. p is the number of contrasts (time points - 1). The correction
# term of stats::mauchly.test reads 3 * (p + 1), the number of time points, at one
# position, and it is kept here so that the P value matches R and car exactly.
.tc_mauchly <- function(M, p, df) {
  if (p < 2) return(list(W = 1, statistic = 0, df = 0, p_value = 1, available = TRUE))
  # The test needs as many error degrees of freedom (subjects - groups) as
  # contrasts; with fewer, M is singular and W is 0, so the test is not available.
  na <- list(W = NA_real_, statistic = NA_real_, df = NA_real_, p_value = NA_real_, available = FALSE)
  if (df < p) return(na)
  logW <- suppressWarnings(log(det(M))) - p * log(sum(diag(M)) / p)
  if (!is.finite(logW)) return(na)
  rho <- 1 - (2 * p^2 + p + 2) / (6 * p * df)
  w2 <- (p + 2) * (p - 1) * (p - 2) * (2 * p^3 + 6 * p^2 + 3 * (p + 1) + 2) /
    (288 * (df * p * rho)^2)
  z <- -df * rho * logW
  f <- p * (p + 1) / 2 - 1
  pr1 <- stats::pchisq(z, f, lower.tail = FALSE)
  pr2 <- stats::pchisq(z, f + 4, lower.tail = FALSE)
  list(W = exp(logW), statistic = z, df = f, p_value = pr1 + w2 * (pr2 - pr1), available = TRUE)
}

.tc_effect_row <- function(r, eps, use_gg, labels) {
  fval <- (r$ss / r$df1) / (r$ss_err / r$df2)
  p <- stats::pf(fval, r$df1, r$df2, lower.tail = FALSE)
  p_gg <- if (r$within) stats::pf(fval, eps * r$df1, eps * r$df2, lower.tail = FALSE) else NA_real_
  corrected <- r$within && use_gg
  data.frame(effect = labels[[r$term]], ss = r$ss, df1 = r$df1, df2 = r$df2, F = fval,
             p_value = p, gg_epsilon = if (r$within) eps else NA_real_, p_gg = p_gg,
             correction = if (corrected) "Greenhouse-Geisser" else "none",
             p_reported = if (corrected) p_gg else p,
             partial_eta2 = r$ss / (r$ss + r$ss_err), stringsAsFactors = FALSE)
}

.tc_n_distinct <- function(v) length(unique(v))

recipe_timecourse <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id; grp <- spec$data$group
  if (is.null(x) || is.null(y)) stop("timecourse needs --x (the time) and --y (the measured value)")
  if (is.null(idc) || !nzchar(idc)) stop("timecourse needs --id (the subject column)")
  has_group <- !is.null(grp) && nzchar(grp)
  for (cn in c(idc, x, y, if (has_group) grp)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(id = as.character(df[[idc]]), time = df[[x]],
                  y = suppressWarnings(as.numeric(df[[y]])),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  n_missing <- n0 - nrow(d)
  if (anyDuplicated(paste(d$id, d$time, sep = "\r"))) {
    stop("timecourse needs one value per subject and time; average any technical replicates first")
  }
  if (any(tapply(d$group, d$id, .tc_n_distinct) > 1)) {
    stop("each subject must belong to one group only")
  }
  # The time points follow --levels, else numeric order, else the number that each
  # label holds ("day 10"), else first appearance (time_level_order, in sankey.R).
  # A numeric time keeps its spacing on the axis unless --levels sets the order.
  d$time_key <- level_labels(d$time)
  t_levels <- unlist(spec$data$levels)
  tord <- time_level_order(d$time_key, t_levels, "time point", "--levels")
  t_key <- tord$levels
  time_numeric <- is.numeric(d$time) && length(t_levels) == 0
  t_lv <- if (time_numeric) as.numeric(t_key) else t_key
  counts <- table(d$id)
  complete_ids <- names(counts)[counts == length(t_key)]
  incomplete <- setdiff(unique(d$id), complete_ids)
  d <- d[d$id %in% complete_ids, , drop = FALSE]
  g_lv <- unique(d$group)
  clean_steps <- c(
    if (n_missing > 0) sprintf("Dropped %d row(s) with a missing subject, time, value or group.", n_missing),
    if (length(incomplete) > 0) sprintf("Dropped %d subject(s) without a value at every one of the %d time points, because the repeated-measures ANOVA needs complete subjects.",
                                        length(incomplete), length(t_key)),
    sprintf("Kept %d subject(s) at %d time points (%d values).", length(complete_ids), length(t_key), nrow(d)),
    sprintf("The time points of '%s' are in %s: %s.", x, tord$rule, paste(t_key, collapse = ", ")))
  if (length(t_key) < 2) stop("timecourse needs at least 2 time points")
  if (length(complete_ids) - length(g_lv) < 1) stop("timecourse needs more complete subjects than groups")
  if (length(unique(d$y)) < 2) {
    stop(sprintf("the outcome '%s' has zero variance: every value of the %d complete subject(s) is %s, so the ANOVA has no variation to test",
                 y, length(complete_ids), level_labels(d$y[1])))
  }

  ids <- sort(complete_ids, method = "radix")
  Y <- matrix(NA_real_, length(ids), length(t_key))
  Y[cbind(match(d$id, ids), match(d$time_key, t_key))] <- d$y
  g <- factor(d$group[match(ids, d$id)], levels = g_lv)

  res <- .tc_mixed_anova(Y, g)
  method <- tolower(spec$test$method %||% "auto")
  if (!method %in% c("auto", "gg", "uncorrected")) {
    stop(sprintf("unknown timecourse test '%s' (use auto, gg or uncorrected)", method))
  }
  use_gg <- switch(method, gg = TRUE, uncorrected = FALSE,
                   auto = length(t_key) > 2 &&
                     (!res$mauchly$available || res$mauchly$p_value < 0.05))
  labels <- list(group = grp %||% "group", time = x,
                 interaction = sprintf("%s:%s", grp %||% "group", x))
  stats_df <- do.call(rbind, lapply(res$rows, .tc_effect_row, eps = res$eps, use_gg = use_gg,
                                    labels = labels))
  rownames(stats_df) <- NULL

  summ <- .tc_summary(Y, g, t_lv)
  p <- .tc_plot(summ, Y, g, ids, t_lv, time_numeric, stats_df, has_group, spec, x, y, grp)
  resid <- Y - res$cells[as.integer(g), , drop = FALSE]
  test_meta <- list(
    name = if (has_group) "mixed ANOVA (group between, time within)" else "one-way repeated-measures ANOVA",
    n_subjects = length(ids), n_times = length(t_key), n_dropped_subjects = length(incomplete),
    sphericity = list(test = "Mauchly", W = res$mauchly$W, statistic = res$mauchly$statistic,
                      df = res$mauchly$df, p_value = res$mauchly$p_value,
                      available = res$mauchly$available),
    gg_epsilon = res$eps, correction_applied = use_gg,
    effects = stats::setNames(lapply(seq_len(nrow(stats_df)), .tc_meta_row, rows = stats_df),
                              stats_df$effect))
  resolved <- list(family = "timecourse", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .tc_methods(stats_df, res, use_gg, method, length(ids), length(t_key),
                             length(incomplete), has_group, spec, y),
       build_script = list(fn = .tc_emit_script, args = list(spec = spec, has_group = has_group)),
       label = slugify(y), df_used = d[, c("id", "group", "time", "y")], clean_steps = clean_steps,
       qc = list(quantity = "residuals from the group-by-time means",
                 panels = list(residuals = as.numeric(resid))),
       width = if (has_group) 6.0 else 5.0, height = 4.4)
}

.tc_meta_row <- function(i, rows) {
  r <- rows[i, ]
  list(F = r$F, df1 = r$df1, df2 = r$df2, p_value = r$p_value, p_gg = r$p_gg,
       p_reported = r$p_reported, partial_eta2 = r$partial_eta2)
}

.tc_summary <- function(Y, g, t_lv) {
  out <- NULL
  for (gl in levels(g)) {
    Yg <- Y[g == gl, , drop = FALSE]
    n <- nrow(Yg)
    sem <- if (n > 1) apply(Yg, 2, stats::sd) / sqrt(n) else rep(NA_real_, ncol(Yg))
    out <- rbind(out, data.frame(group = gl, time = seq_along(t_lv), mean = colMeans(Yg),
                                 sem = sem, stringsAsFactors = FALSE))
  }
  out$group <- factor(out$group, levels = levels(g))
  out
}

.tc_plot <- function(summ, Y, g, ids, t_lv, time_numeric, stats_df, has_group, spec, x, y, grp) {
  pos <- if (time_numeric) as.numeric(t_lv) else seq_along(t_lv)
  summ$x <- pos[summ$time]
  subj <- data.frame(id = rep(ids, times = ncol(Y)), group = rep(g, times = ncol(Y)),
                     x = rep(pos, each = nrow(Y)), y = as.vector(Y))
  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = nlevels(g))
  sub <- paste(sprintf("%s: %s%s", stats_df$effect, vapply(stats_df$p_reported, .tc_p_label, character(1)),
                       ifelse(stats_df$correction == "none", "", " (GG)")), collapse = "; ")
  p <- ggplot2::ggplot(summ, ggplot2::aes(x = x, y = mean, colour = group))
  if (isTRUE(spec$appearance$show_points)) {
    p <- p + ggplot2::geom_line(data = subj, ggplot2::aes(x = x, y = y, group = id, colour = group),
                                alpha = 0.18, linewidth = 0.3)
  }
  dodge <- if (nlevels(g) > 1) 0.012 * diff(range(pos)) else 0
  p <- p +
    ggplot2::geom_line(linewidth = 0.9, position = ggplot2::position_dodge(width = dodge)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sem, ymax = mean + sem), width = 0,
                           linewidth = 0.6, position = ggplot2::position_dodge(width = dodge)) +
    ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(width = dodge)) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = spec$appearance$x_label %||% x, y = spec$appearance$y_label %||% y,
                  colour = if (has_group) grp else NULL, title = spec$appearance$title,
                  subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = if (has_group) "right" else "none",
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
  if (!time_numeric) {
    p <- p + ggplot2::scale_x_continuous(breaks = pos, labels = as.character(t_lv))
  } else if (length(pos) <= 8) {
    p <- p + ggplot2::scale_x_continuous(breaks = pos)
  }
  p
}

.tc_methods <- function(rows, res, use_gg, method, n_subj, n_time, n_drop, has_group, spec, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  y_lab <- spec$appearance$y_label %||% y
  design <- if (has_group) {
    "a mixed ANOVA with the group as a between-subjects factor and time as a within-subjects factor"
  } else "a one-way repeated-measures ANOVA over time"
  drop_txt <- if (n_drop > 0) sprintf(" %d subject(s) without a value at every time point were excluded.", n_drop) else ""
  n_grp <- nrow(res$cells)
  corr_txt <- if (method == "auto" && !res$mauchly$available) {
    sprintf(" Mauchly's test was not available, because it needs at least %d complete subjects (%d subject(s) in %d group(s) at %d time points); the within-subject effects carry the Greenhouse-Geisser correction (epsilon = %.3f).",
            n_grp + n_time - 1, n_subj, n_grp, n_time, res$eps)
  } else if (method == "auto") {
    sprintf(" Sphericity was assessed with Mauchly's test (W = %.3f, P = %s); %s",
            res$mauchly$W, sig_text(res$mauchly$p_value, 2),
            if (use_gg) sprintf("because it was rejected, the within-subject effects carry the Greenhouse-Geisser correction (epsilon = %.3f).", res$eps)
            else "it was not rejected, so no correction was applied.")
  } else if (use_gg) {
    sprintf(" The within-subject effects carry the Greenhouse-Geisser correction (epsilon = %.3f).", res$eps)
  } else " No sphericity correction was applied."
  eff <- paste(sprintf("%s, F(%s, %s) = %.2f, %s, partial eta squared = %.3f",
                       rows$effect, .tc_df_label(rows$df1, rows$gg_epsilon, rows$correction),
                       .tc_df_label(rows$df2, rows$gg_epsilon, rows$correction), rows$F,
                       vapply(rows$p_reported, .tc_p_label, character(1)), rows$partial_eta2),
               collapse = "; ")
  paste0(
    sprintf("%s was analysed with %s (n = %d subjects, %d time points).", cap_first(y_lab), design, n_subj, n_time),
    drop_txt, corr_txt, " ", eff, ". The figure shows the mean +/- SEM at each time point, with one thin line per subject.",
    sprintf(" Computed in R %s and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.tc_df_label <- function(df, eps, correction) {
  ifelse(correction == "none", sprintf("%g", df), sprintf("%.2f", df * eps))
}

.tc_emit_script <- function(spec, has_group, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id; grp <- spec$data$group
  fml <- if (has_group) {
    sprintf("y ~ group * time_f + Error(id / time_f)")
  } else "y ~ time_f + Error(id / time_f)"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder. The ANOVA uses complete subjects.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(id = factor(raw[["%s"]]), time = raw[["%s"]], y = as.numeric(raw[["%s"]]),', idc, x, y),
    sprintf('                group = %s)', if (has_group) sprintf('factor(raw[["%s"]])', grp) else 'factor("all")'),
    "d <- d[complete.cases(d), ]",
    "n_times <- length(unique(d$time))",
    "d <- d[d$id %in% names(which(table(d$id) == n_times)), ]",
    "d$id <- droplevels(d$id); d$time_f <- factor(d$time)",
    sprintf("print(summary(aov(%s, data = d)))", fml),
    "m <- aggregate(y ~ group + time, data = d, FUN = mean)",
    "p <- ggplot(m, aes(time, y, colour = group)) + geom_line() + geom_point() +",
    sprintf('  labs(x = "%s", y = "%s") + theme_classic()', spec$appearance$x_label %||% x,
            spec$appearance$y_label %||% y),
    sprintf('ggsave("%s.pdf", p, width = 6.0, height = 4.4)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
