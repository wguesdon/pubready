# Recipe: swimmer plot of the time each patient spends on treatment or on study.
# Input is one row per patient: the patient (--id), the duration (--time), the
# arm (--group, optional) and the status (--event, optional: 1 = the course
# ended, 0 = ongoing). One horizontal bar per patient, sorted from the longest
# course at the top, with an arrow at the end of each ongoing course. Columns
# named in --markers hold the time of an event on the course, for example the
# start of a response or a progression, and each one is drawn with its own shape.
#
# The plot is descriptive, so there is no hypothesis test. The stats table holds
# one row per patient, and the metadata gives the median duration and the count
# of ongoing courses per arm.

.sw_count_present <- function(col, d) sum(is.finite(d[[col]]))

.sw_arm_summary <- function(g, d) {
  sel <- d$group == g
  list(n = sum(sel), median_duration = stats::median(d$duration[sel]),
       ongoing = sum(d$ongoing[sel]))
}

recipe_swimmer <- function(df, spec) {
  idc <- spec$data$id; tcol <- spec$data$time; ev <- spec$data$event; grp <- spec$data$group
  markers <- spec$data$markers %||% character(0)
  if (is.null(idc) || !nzchar(idc)) stop("swimmer needs --id (the patient column)")
  if (is.null(tcol) || !nzchar(tcol)) stop("swimmer needs --time (the duration on treatment)")
  has_event <- !is.null(ev) && nzchar(ev)
  has_group <- !is.null(grp) && nzchar(grp)
  for (cn in c(idc, tcol, if (has_event) ev, if (has_group) grp, markers)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  n0 <- nrow(df)
  d <- data.frame(patient = as.character(df[[idc]]),
                  duration = suppressWarnings(as.numeric(df[[tcol]])),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  status <- if (has_event) {
    event_codes(df[[ev]], ev, "1 = the course ended, 0 = ongoing")
  } else rep(1L, n0)
  d$ongoing <- status == 0L
  for (m in markers) d[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  keep <- !is.na(d$patient) & is.finite(d$duration) & d$duration >= 0 & !is.na(d$group) &
    !is.na(d$ongoing)
  d <- d[keep, , drop = FALSE]
  if (anyDuplicated(d$patient)) stop("swimmer needs one row per patient")
  if (nrow(d) < 1) stop("swimmer needs at least one patient with a duration")
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing patient, duration, status or arm.", sum(!keep)),
    sprintf("Kept %d of %d patient(s).", nrow(d), n0))

  d <- d[order(d$duration, d$patient, method = "radix"), , drop = FALSE]
  d$rank <- seq_len(nrow(d))
  # The arms are sorted, as the other recipes sort their groups; "radix" gives
  # the C-locale order that the Python engine uses.
  g_lv <- if (has_group) sort(unique(as.character(df[[grp]][keep])), method = "radix") else "all"
  d$group <- factor(d$group, levels = g_lv)
  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = length(g_lv))
  span <- max(d$duration)
  arrow_len <- 0.04 * span

  p <- ggplot2::ggplot(d) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = duration, y = rank, yend = rank, colour = group),
                         linewidth = 3.2, lineend = "butt")
  ong <- d[d$ongoing, , drop = FALSE]
  if (nrow(ong) > 0) {
    p <- p + ggplot2::geom_segment(
      data = ong, ggplot2::aes(x = duration, xend = duration + arrow_len, y = rank, yend = rank),
      colour = "grey20", linewidth = 0.5,
      arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"))
  }
  if (length(markers) > 0) {
    mk <- do.call(rbind, lapply(markers, .sw_marker_rows, d = d))
    if (!is.null(mk) && nrow(mk) > 0) {
      mk$marker <- factor(mk$marker, levels = markers)
      p <- p + ggplot2::geom_point(data = mk, ggplot2::aes(x = time, y = rank, shape = marker),
                                   size = 2.2, fill = "white", colour = "black") +
        ggplot2::scale_shape_manual(values = c(21, 24, 23, 22, 25)[seq_along(markers)], drop = FALSE)
    }
  }
  p <- p +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_y_continuous(breaks = if (nrow(d) <= 40) d$rank else NULL,
                                labels = if (nrow(d) <= 40) d$patient else NULL,
                                expand = ggplot2::expansion(add = 0.8)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
    ggplot2::labs(x = spec$appearance$x_label %||% tcol,
                  y = spec$appearance$y_label %||% "Patient",
                  colour = if (has_group) grp else NULL, shape = NULL,
                  title = spec$appearance$title,
                  subtitle = sprintf("%d patients; %d ongoing (arrow)", nrow(d), sum(d$ongoing))) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   axis.text.y = ggplot2::element_text(size = 7),
                   legend.text = ggplot2::element_text(size = 9),
                   legend.title = ggplot2::element_text(size = 10),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))
  if (!has_group) p <- p + ggplot2::guides(colour = "none")

  stats_df <- d[rev(seq_len(nrow(d))), c("patient", "group", "duration", "ongoing", markers), drop = FALSE]
  stats_df$group <- as.character(stats_df$group)
  if (!has_group) stats_df$group <- NULL
  rownames(stats_df) <- NULL
  test_meta <- list(
    name = "Swimmer plot (duration per patient)",
    n_patients = nrow(d), n_ongoing = sum(d$ongoing),
    arms = stats::setNames(lapply(g_lv, .sw_arm_summary, d = d), g_lv),
    markers = stats::setNames(lapply(markers, .sw_count_present, d = d), markers))
  resolved <- list(family = "swimmer", label = test_meta$name)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .sw_methods(d, g_lv, markers, tcol, has_event),
       build_script = list(fn = .sw_emit_script,
                           args = list(spec = spec, markers = markers, has_event = has_event,
                                       has_group = has_group)),
       label = slugify(tcol), df_used = d, clean_steps = clean_steps,
       width = 6.4, height = max(3.6, min(10, 1.6 + 0.18 * nrow(d))))
}

.sw_marker_rows <- function(m, d) {
  sel <- is.finite(d[[m]])
  if (!any(sel)) return(NULL)
  data.frame(rank = d$rank[sel], time = d[[m]][sel], marker = m, stringsAsFactors = FALSE)
}

.sw_methods <- function(d, g_lv, markers, tcol, has_event) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  arms <- paste(vapply(g_lv, .sw_arm_text, character(1), d = d), collapse = "; ")
  mk <- if (length(markers) > 0) {
    sprintf(" Symbols mark %s.", paste(sprintf("%s (%d patient(s))", markers,
                                                vapply(markers, .sw_count_present, integer(1), d = d)),
                                        collapse = ", "))
  } else ""
  paste0(
    sprintf("Each of the %d patients is drawn as one horizontal bar whose length is the value of %s, sorted from the longest course at the top.",
            nrow(d), tcol),
    if (has_event) " An arrow marks a course that was ongoing at the data cut-off." else "",
    mk, " ", arms, ".",
    " The figure is descriptive, so no statistical test is applied.",
    sprintf(" Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.sw_arm_text <- function(g, d) {
  s <- .sw_arm_summary(g, d)
  sprintf("%s%d patient(s), median duration %s, %d ongoing",
          if (g == "all") "" else paste0(g, ": "), s$n, sprintf("%.3g", s$median_duration), s$ongoing)
}

.sw_emit_script <- function(spec, markers, has_event, has_group, in_name, fig_stub) {
  idc <- spec$data$id; tcol <- spec$data$time; ev <- spec$data$event; grp <- spec$data$group
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(patient = as.character(raw[["%s"]]), duration = as.numeric(raw[["%s"]]),', idc, tcol),
    sprintf('                group = %s, ongoing = %s)',
            if (has_group) sprintf('as.character(raw[["%s"]])', grp) else '"all"',
            if (has_event) sprintf('raw[["%s"]] == 0', ev) else "FALSE"),
    "d <- d[order(d$duration), ]; d$rank <- seq_len(nrow(d))",
    "p <- ggplot(d) + geom_segment(aes(x = 0, xend = duration, y = rank, yend = rank, colour = group),",
    "                              linewidth = 3) +",
    "  geom_point(data = d[d$ongoing, ], aes(x = duration, y = rank), shape = 62, size = 4) +",
    "  scale_y_continuous(breaks = d$rank, labels = d$patient) +",
    sprintf('  labs(x = "%s", y = "Patient") + theme_classic()', spec$appearance$x_label %||% tcol),
    sprintf('ggsave("%s.pdf", p, width = 6.4, height = %s)', fig_stub,
            "max(3.6, min(10, 1.6 + 0.18 * nrow(d)))"),
    'cat("Reproduced figure.\\n")'
  )
}
