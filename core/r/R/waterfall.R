# Recipe: waterfall plot of the best response of each patient. One bar per
# patient, sorted from the largest growth on the left to the deepest shrinkage
# on the right, with dashed lines at the RECIST 1.1 thresholds (+20% and -30% by
# default, --pd_threshold and --pr_threshold).
#
# Two inputs are accepted. With --x (the time), the data is one row per visit, as
# for spider_response, and each patient's best change is the smallest percentage
# change from their own baseline over the visits after it. The category then
# follows spider_response: a partial response when the best change reaches the
# PR line, else progressive disease when the worst change reaches the PD line,
# else stable disease. Without --x, the data is one row per patient and --y is
# already the best percentage change, so the category reads that value alone.
#
# The plot is descriptive, so there is no hypothesis test. The stats table holds
# one row per patient, and the metadata counts the categories per group.

# The legend labels wrap at this many characters, as in the Python engine.
.wf_wrap <- 24

# The figure width: 0.16 in per bar between 5 and 9 in, wider when the legend
# needs more than 2 in, so a long label widens the figure instead of squeezing
# the panel. The Python engine uses the same formula.
waterfall_width <- function(n, legend_labels, title) {
  max(5.0, min(9.0, 3.2 + 0.16 * n)) +
    max(0, legend_width(wrap_labels(legend_labels, .wf_wrap), title) - 2.0)
}

# The figure height: 4.4 in, taller when the patient labels under the bars
# (drawn up to 30 patients) need more than 0.5 in at 0.055 in per character.
waterfall_height <- function(patients) {
  chars <- if (length(patients) <= 30) max(nchar(patients), 0L) else 0L
  4.4 + max(0, 0.055 * chars - 0.5)
}

.wf_first <- function(v) v[1]
.wf_last <- function(v) v[length(v)]

.wf_category <- function(best, worst, pr_cut, pd_cut) {
  ifelse(best <= pr_cut, "partial response",
         ifelse(worst >= pd_cut, "progressive disease", "stable disease"))
}

# Best and worst change from baseline per patient from one row per visit.
.wf_from_visits <- function(d) {
  d <- d[order(d$id, d$time, method = "radix"), , drop = FALSE]
  base <- tapply(d$y, d$id, .wf_first)
  usable <- names(base)[is.finite(base) & base != 0]
  dropped_zero <- length(base) - length(usable)
  d <- d[d$id %in% usable, , drop = FALSE]
  d$change <- 100 * (d$y - base[d$id]) / base[d$id]
  post <- d[duplicated(d$id), , drop = FALSE]
  no_follow <- setdiff(unique(d$id), unique(post$id))
  per <- data.frame(patient = names(tapply(post$change, post$id, min)),
                    best_change_pct = as.numeric(tapply(post$change, post$id, min)),
                    worst_change_pct = as.numeric(tapply(post$change, post$id, max)),
                    stringsAsFactors = FALSE)
  per$group <- as.character(tapply(d$group, d$id, .wf_first)[per$patient])
  list(per = per, dropped_zero = dropped_zero, no_follow = length(no_follow))
}

recipe_waterfall <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id; grp <- spec$data$group
  if (is.null(y)) stop("waterfall needs --y (the measured value, or the best percentage change)")
  if (is.null(idc) || !nzchar(idc)) stop("waterfall needs --id (the patient column)")
  has_time <- !is.null(x) && nzchar(x)
  has_group <- !is.null(grp) && nzchar(grp)
  for (cn in c(idc, y, if (has_time) x, if (has_group) grp)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }
  pd_cut <- spec$appearance$pd_threshold %||% 20
  pr_cut <- spec$appearance$pr_threshold %||% -30

  if (has_time) check_numeric_time(df[[x]], x)
  n0 <- nrow(df)
  d <- data.frame(id = as.character(df[[idc]]),
                  time = if (has_time) suppressWarnings(as.numeric(df[[x]])) else 0,
                  y = suppressWarnings(as.numeric(df[[y]])),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  n_missing <- n0 - nrow(d)
  if (has_time) {
    check_one_row_per_time(d$id, d$time, "waterfall", x)
    late_step <- late_baseline_step(d$id, d$time, x)
    fv <- .wf_from_visits(d)
    per <- fv$per
    per$category <- .wf_category(per$best_change_pct, per$worst_change_pct, pr_cut, pd_cut)
    extra <- c(
      late_step,
      if (fv$dropped_zero > 0) sprintf("Dropped %d patient(s) whose baseline was zero, because a percentage change needs a non-zero baseline.", fv$dropped_zero),
      if (fv$no_follow > 0) sprintf("Dropped %d patient(s) with no visit after the baseline.", fv$no_follow))
  } else {
    if (anyDuplicated(d$id)) {
      stop("without --x, waterfall needs one row per patient; pass --x (the time) for one row per visit")
    }
    per <- data.frame(patient = d$id, best_change_pct = d$y, group = d$group,
                      stringsAsFactors = FALSE)
    per$category <- .wf_category(per$best_change_pct, per$best_change_pct, pr_cut, pd_cut)
    extra <- NULL
  }
  if (nrow(per) < 1) stop("waterfall needs at least one patient with a response")
  clean_steps <- c(
    if (n_missing > 0) sprintf("Dropped %d row(s) with a missing patient, time, value or group.", n_missing),
    extra, sprintf("Kept %d patient(s).", nrow(per)))

  ord <- order(-per$best_change_pct, per$patient, method = "radix")
  per <- per[ord, , drop = FALSE]
  rownames(per) <- NULL
  per$rank <- seq_len(nrow(per))
  cats <- c("partial response", "stable disease", "progressive disease")
  if (has_group) {
    fill_var <- factor(per$group, levels = unique(d$group))
    legend_title <- grp
  } else {
    fill_var <- factor(per$category, levels = cats[cats %in% per$category])
    legend_title <- "Response"
  }
  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = nlevels(fill_var))
  per_plot <- data.frame(rank = per$rank, best = per$best_change_pct, fill = fill_var)
  n_pr <- sum(per$category == "partial response")

  p <- ggplot2::ggplot(per_plot, ggplot2::aes(x = rank, y = best, fill = fill)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(pd_cut, pr_cut), linetype = "dashed", colour = "grey45",
                        linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, labels = wrap_labels(levels(fill_var), .wf_wrap)) +
    ggplot2::scale_x_continuous(breaks = if (nrow(per) <= 30) per$rank else NULL,
                                labels = if (nrow(per) <= 30) per$patient else NULL,
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(x = spec$appearance$x_label %||% "Patient",
                  y = spec$appearance$y_label %||% "Best change from baseline (%)",
                  fill = legend_title, title = spec$appearance$title,
                  subtitle = sprintf("%d patients; %d (%.0f%%) reach the %+g%% line", nrow(per), n_pr,
                                     100 * n_pr / nrow(per), pr_cut)) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  stats_df <- per[, c("patient", "group", "best_change_pct",
                      if (has_time) "worst_change_pct", "category"), drop = FALSE]
  names(stats_df)[names(stats_df) == "category"] <- "category_from_change"
  if (!has_group) stats_df$group <- NULL
  test_meta <- list(
    name = "Waterfall plot (best change from baseline)",
    n_patients = nrow(per), pd_threshold = pd_cut, pr_threshold = pr_cut,
    input = if (has_time) "one row per visit" else "one row per patient",
    counts = .wf_counts(per, cats))
  resolved <- list(family = "waterfall", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .wf_methods(per, cats, has_time, pd_cut, pr_cut, spec, y),
       build_script = list(fn = .wf_emit_script,
                           args = list(spec = spec, has_time = has_time, pd_cut = pd_cut,
                                       pr_cut = pr_cut)),
       label = slugify(y), df_used = d, clean_steps = clean_steps,
       width = waterfall_width(nrow(per), levels(fill_var), legend_title),
       height = waterfall_height(per$patient))
}

.wf_count_group <- function(gl, per, cats) {
  sel <- per$category[per$group == gl]
  as.list(stats::setNames(as.integer(table(factor(sel, levels = cats))), gsub(" ", "_", cats)))
}

.wf_counts <- function(per, cats) {
  g_lv <- unique(per$group)
  stats::setNames(lapply(g_lv, .wf_count_group, per = per, cats = cats), g_lv)
}

.wf_methods <- function(per, cats, has_time, pd_cut, pr_cut, spec, y) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  tab <- table(factor(per$category, levels = cats))
  src <- if (has_time) {
    sprintf("The best response of each patient was the smallest percentage change in %s from that patient's own baseline over the visits after it. ", y)
  } else sprintf("The best percentage change of each patient was read from %s. ", y)
  paste0(
    src,
    sprintf("The %d patients are drawn as bars sorted from the largest increase to the deepest decrease, with dashed reference lines at %+g%% and %+g%%, the RECIST 1.1 thresholds for progressive disease and for a partial response. ",
            nrow(per), pd_cut, pr_cut),
    sprintf("On the change from baseline alone, %d patient(s) reached a partial response, %d were stable and %d reached progressive disease. ",
            tab[["partial response"]], tab[["stable disease"]], tab[["progressive disease"]]),
    "A full RECIST 1.1 assessment also reads the non-target lesions and any new lesion. ",
    "The figure is descriptive, so no statistical test is applied.",
    sprintf(" Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.wf_emit_script <- function(spec, has_time, pd_cut, pr_cut, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id
  body <- if (has_time) {
    c(sprintf('raw <- raw[complete.cases(raw[, c("%s", "%s", "%s")]), ]', idc, x, y),
      sprintf('raw <- raw[order(raw[["%s"]], raw[["%s"]]), ]', idc, x),
      sprintf('base <- tapply(raw[["%s"]], raw[["%s"]], function(v) v[1])', y, idc),
      sprintf('raw$change <- 100 * (raw[["%s"]] - base[raw[["%s"]]]) / base[raw[["%s"]]]', y, idc, idc),
      sprintf('post <- raw[duplicated(raw[["%s"]]), ]', idc),
      sprintf('best <- tapply(post$change, post[["%s"]], min)', idc))
  } else {
    c(sprintf('best <- setNames(raw[["%s"]], raw[["%s"]])', y, idc))
  }
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    body,
    "d <- data.frame(patient = names(best), best = as.numeric(best))",
    "d <- d[order(-d$best), ]; d$rank <- seq_len(nrow(d))",
    "p <- ggplot(d, aes(rank, best)) + geom_col() + geom_hline(yintercept = 0) +",
    sprintf('  geom_hline(yintercept = c(%g, %g), linetype = "dashed") +', pd_cut, pr_cut),
    '  labs(x = "Patient", y = "Best change from baseline (%)") + theme_classic()',
    sprintf('ggsave("%s.pdf", p, width = 6, height = 4.4)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
