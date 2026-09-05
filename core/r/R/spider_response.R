# Recipe: spider plot of the change from baseline of each patient over time.
# One line per patient (--id), the time on --x and the measured value on --y.
# Each value becomes a percentage change from that patient's own baseline, which
# is the row with the smallest time. Two dashed reference lines carry the RECIST
# 1.1 thresholds: +20% (--pd_threshold) and -30% (--pr_threshold).
#
# The plot is descriptive, so there is no hypothesis test. The stats table holds
# one row per patient with the baseline value, the best and the worst change and
# the category those two put the patient in. The best and the worst are measured
# over the visits after the baseline, because every patient reads zero at the
# baseline itself. That categorisation reads the change from baseline alone. A
# full RECIST 1.1 assessment also reads the nadir, the non-target lesions and any
# new lesion, so the column is named for what it measures.

recipe_spider_response <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id; grp <- spec$data$group
  if (is.null(x) || is.null(y)) {
    stop("spider_response needs --x (the time) and --y (the measured value)")
  }
  if (is.null(idc) || !nzchar(idc)) {
    stop("spider_response needs --id (the patient column)")
  }
  cols <- c(idc, x, y, if (!is.null(grp) && nzchar(grp)) grp)
  for (cn in cols) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  d <- df[, cols, drop = FALSE]
  d[[x]] <- suppressWarnings(as.numeric(d[[x]]))
  d[[y]] <- suppressWarnings(as.numeric(d[[y]]))
  if (all(is.na(d[[x]]))) {
    stop(sprintf("column '%s' holds no number, and the time has to be numeric", x))
  }
  n_rows0 <- nrow(d)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d[[idc]] <- as.character(d[[idc]])
  d <- d[order(d[[idc]], d[[x]]), , drop = FALSE]

  baseline <- stats::aggregate(d[[y]], by = list(id = d[[idc]]),
                               FUN = function(v) v[1])
  names(baseline) <- c("id", "baseline")
  usable <- baseline$id[is.finite(baseline$baseline) & baseline$baseline != 0]
  dropped <- setdiff(baseline$id, usable)
  d <- d[d[[idc]] %in% usable, , drop = FALSE]
  baseline <- baseline[baseline$id %in% usable, , drop = FALSE]
  if (nrow(d) < 2) {
    stop("spider_response needs at least one patient with a non-zero baseline")
  }

  d$change <- 100 * (d[[y]] - baseline$baseline[match(d[[idc]], baseline$id)]) /
    baseline$baseline[match(d[[idc]], baseline$id)]

  pd_cut <- spec$appearance$pd_threshold %||% 20
  pr_cut <- spec$appearance$pr_threshold %||% -30

  # The response is measured over the visits after the baseline. Every patient
  # reads zero at the baseline, so a patient who only grows would otherwise
  # carry a best change of zero rather than the smallest growth they reached.
  post <- d[duplicated(d[[idc]]), , drop = FALSE]
  no_follow_up <- setdiff(unique(d[[idc]]), unique(post[[idc]]))
  if (length(no_follow_up) > 0) {
    d <- d[!d[[idc]] %in% no_follow_up, , drop = FALSE]
    post <- post[!post[[idc]] %in% no_follow_up, , drop = FALSE]
  }
  if (nrow(post) < 1) {
    stop("spider_response needs at least one patient with a visit after the baseline")
  }

  best  <- stats::aggregate(post$change, by = list(id = post[[idc]]), FUN = min)
  worst <- stats::aggregate(post$change, by = list(id = post[[idc]]), FUN = max)
  last  <- stats::aggregate(post$change, by = list(id = post[[idc]]),
                            FUN = function(v) v[length(v)])
  visits <- stats::aggregate(d$change, by = list(id = d[[idc]]), FUN = length)
  names(best) <- c("id", "best"); names(worst) <- c("id", "worst")
  names(last) <- c("id", "last"); names(visits) <- c("id", "visits")

  category <- ifelse(best$best <= pr_cut, "partial response",
                     ifelse(worst$worst[match(best$id, worst$id)] >= pd_cut,
                            "progressive disease", "stable disease"))
  category <- factor(category, levels = c("partial response", "stable disease",
                                          "progressive disease"))

  clean_steps <- c(
    sprintf("Kept %d patient(s) across %d measurement(s).",
            nrow(best), nrow(d)),
    if (n_rows0 > nrow(d)) {
      sprintf("Dropped %d row(s) with a missing time, value or patient.",
              n_rows0 - nrow(d))
    } else NULL,
    if (length(dropped) > 0) {
      sprintf("Dropped %d patient(s) whose baseline was zero or absent, because a percentage change needs a non-zero baseline.",
              length(dropped))
    } else NULL,
    if (length(no_follow_up) > 0) {
      sprintf("Dropped %d patient(s) with no visit after the baseline, because a response needs a second measurement.",
              length(no_follow_up))
    } else NULL)

  plotdf <- data.frame(id = d[[idc]], time = d[[x]], change = d$change,
                       stringsAsFactors = FALSE)
  if (!is.null(grp) && nzchar(grp)) {
    plotdf$colour <- as.character(d[[grp]])
    legend_title <- grp
  } else {
    plotdf$colour <- as.character(category[match(plotdf$id, best$id)])
    legend_title <- "Response"
  }
  # A group keeps the order it appears in, so a control arm stays left of a
  # treated arm. A response category keeps the clinical order instead, so the
  # legend reads from the best outcome to the worst whatever the data holds.
  levels_seen <- if (!is.null(grp) && nzchar(grp)) {
    unique(plotdf$colour)
  } else {
    levels(category)[levels(category) %in% plotdf$colour]
  }
  plotdf$colour <- factor(plotdf$colour, levels = levels_seen)
  pal <- spec$appearance$palette %||% pubplot_palette
  pal <- rep(pal, length.out = length(levels_seen))

  p <- ggplot2::ggplot(plotdf, ggplot2::aes(x = time, y = change,
                                            group = id, colour = colour)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = pd_cut, linetype = "dashed",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = pr_cut, linetype = "dashed",
                        colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
    ggplot2::geom_point(size = 1.8, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% "Change from baseline (%)",
                  colour = legend_title,
                  title = spec$appearance$title) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right")

  stats_df <- data.frame(
    patient = best$id,
    baseline = baseline$baseline[match(best$id, baseline$id)],
    visits = visits$visits[match(best$id, visits$id)],
    best_change_pct = best$best,
    worst_change_pct = worst$worst[match(best$id, worst$id)],
    last_change_pct = last$last[match(best$id, last$id)],
    category_from_change = as.character(category),
    stringsAsFactors = FALSE)
  if (!is.null(grp) && nzchar(grp)) {
    arm <- d[!duplicated(d[[idc]]), c(idc, grp)]
    stats_df$group <- arm[[grp]][match(stats_df$patient, arm[[idc]])]
  }
  stats_df <- stats_df[order(stats_df$best_change_pct), , drop = FALSE]

  counts <- table(category)
  test_meta <- list(
    name = "Spider plot (change from baseline)",
    n_patients = nrow(best),
    n_measurements = nrow(d),
    n_dropped_patients = length(dropped),
    n_without_follow_up = length(no_follow_up),
    pd_threshold = pd_cut,
    pr_threshold = pr_cut,
    partial_response = unname(counts[["partial response"]]),
    stable_disease = unname(counts[["stable disease"]]),
    progressive_disease = unname(counts[["progressive disease"]]))

  resolved <- list(family = "spider", label = "Spider plot")
  methods  <- .spider_methods(nrow(best), nrow(d), pd_cut, pr_cut, counts, spec)
  build_script <- function(in_name, fig_stub) {
    .spider_emit_script(spec, in_name, fig_stub, pd_cut, pr_cut)
  }

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(y), df_used = d, clean_steps = clean_steps,
       width = 5.8, height = 4.2)
}

.spider_methods <- function(n_patients, n_rows, pd_cut, pr_cut, counts, spec) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  y_lab <- spec$appearance$y_label %||% spec$data$y
  paste0(
    sprintf("%s was expressed as a percentage change from each patient's own baseline, which is the first measurement in time, and drawn as one line per patient (n = %d patients, %d measurements). ",
            cap_first(y_lab), n_patients, n_rows),
    sprintf("Dashed reference lines mark %+g%% and %+g%%, the RECIST 1.1 thresholds for progressive disease and for a partial response. ",
            pd_cut, pr_cut),
    sprintf("On the change from baseline alone, %d patient(s) reached a partial response, %d were stable and %d reached progressive disease. ",
            counts[["partial response"]], counts[["stable disease"]],
            counts[["progressive disease"]]),
    "That categorisation reads the target lesion change from baseline only. A full RECIST 1.1 assessment also reads the nadir, the non-target lesions and any new lesion. ",
    "The figure is descriptive, so no statistical test is applied.",
    sprintf(" Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.spider_emit_script <- function(spec, in_name, fig_stub, pd_cut, pr_cut) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('raw[["%s"]] <- as.numeric(raw[["%s"]]); raw[["%s"]] <- as.numeric(raw[["%s"]])',
            x, x, y, y),
    sprintf('raw <- raw[complete.cases(raw[, c("%s","%s","%s")]), ]', idc, x, y),
    sprintf('raw <- raw[order(raw[["%s"]], raw[["%s"]]), ]', idc, x),
    sprintf('base <- tapply(raw[["%s"]], raw[["%s"]], function(v) v[1])', y, idc),
    sprintf('raw$change <- 100 * (raw[["%s"]] - base[raw[["%s"]]]) / base[raw[["%s"]]]',
            y, idc, idc),
    "",
    sprintf('p <- ggplot(raw, aes(.data[["%s"]], change, group = .data[["%s"]])) +',
            x, idc),
    '  geom_hline(yintercept = 0, colour = "grey30") +',
    sprintf('  geom_hline(yintercept = %g, linetype = "dashed", colour = "grey55") +',
            pd_cut),
    sprintf('  geom_hline(yintercept = %g, linetype = "dashed", colour = "grey55") +',
            pr_cut),
    "  geom_line() + geom_point(size = 1.8) +",
    sprintf('  labs(x = "%s", y = "Change from baseline (%%)") + theme_classic()', x),
    sprintf('ggsave("%s.pdf", p, width = 5.8, height = 4.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
