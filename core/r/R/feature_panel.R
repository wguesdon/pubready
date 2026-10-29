# Recipe: small multiples, one panel per feature (a gene, a protein, a
# cytokine), each with a test across the groups, and the Benjamini-Hochberg
# correction across the features. Input is long (--x group --y value --facet
# feature) or wide (--x group, and every other numeric column is a feature).
#
# One test family covers the whole panel. Two groups take the Welch t-test or
# the Mann-Whitney U test; three or more take the one-way ANOVA or the
# Kruskal-Wallis test. --test auto runs the Shapiro-Wilk test on every group of
# every feature: when all are consistent with normality every feature gets the
# parametric test, otherwise every feature gets the rank test. --test welch_t,
# wilcoxon, anova or kruskal forces the test.

.fp_max_features  <- 36
.fp_max_qc_panels <- 48
.fp_min_n         <- 2
.fp_labels <- c(welch_t = "Welch t-test", wilcoxon = "Mann-Whitney U test",
                anova = "one-way ANOVA", kruskal = "Kruskal-Wallis test")
.fp_effects <- c(welch_t = "Cohen's d", wilcoxon = "rank-biserial r",
                 anova = "eta squared", kruskal = "epsilon squared")

# Significance stars with the cutpoints of rstatix::add_significance.
.fp_signif <- function(p) {
  out <- as.character(cut(p, breaks = c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                          labels = c("****", "***", "**", "*", "ns"), right = FALSE))
  out[is.na(out)] <- "NA"
  out
}

.fp_na_on_error <- function(e) NULL

# One wide column as long rows.
.fp_stack_col <- function(col, df, x) {
  data.frame(feature = col, group = df[[x]], value = suppressWarnings(as.numeric(df[[col]])),
             stringsAsFactors = FALSE)
}

# The group levels: numeric groups in numeric order, text in C-locale order.
.fp_group_levels <- function(g) {
  if (is.numeric(g)) return(as.character(sort(unique(g))))
  sort(unique(as.character(g)), method = "radix")
}

.fp_is_feature_col <- function(v) is.numeric(v)

.fp_long <- function(df, spec) {
  x <- spec$data$x; fcol <- spec$data$facet; y <- spec$data$y
  if (is.null(x) || !nzchar(x)) stop("feature_panel needs --x (the group column)")
  if (!x %in% names(df)) stop(sprintf("column '%s' not found in data", x))
  if (!is.null(fcol) && nzchar(fcol)) {
    if (is.null(y) || !nzchar(y)) stop("feature_panel in long form needs --y (the value column)")
    for (cn in c(fcol, y)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
    d <- data.frame(feature = as.character(df[[fcol]]), group = df[[x]],
                    value = suppressWarnings(as.numeric(df[[y]])), stringsAsFactors = FALSE)
    layout <- "long"
  } else {
    skip <- c(x, spec$data$id, spec$data$label)
    cols <- setdiff(names(df)[vapply(df, .fp_is_feature_col, logical(1))], skip)
    if (length(cols) == 0) stop("feature_panel found no numeric feature column; pass --facet for long data")
    d <- do.call(rbind, lapply(cols, .fp_stack_col, df = df, x = x))
    layout <- "wide"
  }
  n0 <- nrow(d)
  keep <- !is.na(d$feature) & nzchar(d$feature) & !is.na(d$group) & is.finite(d$value)
  d <- d[keep, , drop = FALSE]
  glv <- .fp_group_levels(d$group)
  d$group <- as.character(d$group)
  feats <- unique(d$feature)
  steps <- c(
    sprintf("Read %s data: %d feature(s), %d group(s).", layout, length(feats), length(glv)),
    if (sum(!keep) > 0) sprintf("Dropped %d value(s) with a missing feature, group or value (%d -> %d).",
                                sum(!keep), n0, nrow(d)))
  list(d = d, features = feats, groups = glv, steps = steps, layout = layout)
}

.fp_resolve <- function(method, n_groups, all_normal) {
  method <- tolower(method %||% "auto")
  if (method %in% c("mann_whitney", "mwu", "wilcox")) method <- "wilcoxon"
  if (method %in% c("kw", "kruskal_wallis")) method <- "kruskal"
  if (method == "auto") {
    if (isTRUE(all_normal)) return(if (n_groups == 2) "welch_t" else "anova")
    return(if (n_groups == 2) "wilcoxon" else "kruskal")
  }
  if (!method %in% names(.fp_labels)) {
    stop(sprintf("feature_panel: unknown test '%s'; use auto, welch_t, wilcoxon, anova or kruskal", method))
  }
  if (method %in% c("welch_t", "wilcoxon") && n_groups != 2) {
    stop(sprintf("feature_panel: --test %s needs exactly 2 groups; found %d. Use anova or kruskal.",
                 method, n_groups))
  }
  method
}

.fp_shapiro_normal <- function(v) isTRUE(shapiro_safe(v)$normal)

.fp_cohens_d <- function(a, b) {
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * stats::var(a) + (nb - 1) * stats::var(b)) / (na + nb - 2))
  (mean(a) - mean(b)) / sp
}

# The test of one feature: statistic, df1, df2, p and the effect size.
.fp_run_test <- function(by_group, test) {
  if (test == "welch_t") {
    a <- by_group[[1]]; b <- by_group[[2]]
    tt <- stats::t.test(a, b, var.equal = FALSE)
    return(c(unname(tt$statistic), unname(tt$parameter), NA_real_, tt$p.value, .fp_cohens_d(a, b)))
  }
  if (test == "wilcoxon") {
    a <- by_group[[1]]; b <- by_group[[2]]
    wt <- suppressWarnings(stats::wilcox.test(a, b))
    w <- unname(wt$statistic)
    return(c(w, NA_real_, NA_real_, wt$p.value, rank_biserial_unpaired(w, length(a), length(b))))
  }
  v <- unlist(by_group, use.names = FALSE)
  g <- factor(rep(names(by_group), lengths(by_group)), levels = names(by_group))
  if (test == "anova") {
    ow <- stats::oneway.test(v ~ g, var.equal = TRUE)
    ssb <- sum(lengths(by_group) * (vapply(by_group, mean, numeric(1)) - mean(v))^2)
    return(c(unname(ow$statistic), unname(ow$parameter[1]), unname(ow$parameter[2]), ow$p.value,
             ssb / sum((v - mean(v))^2)))
  }
  kt <- stats::kruskal.test(v, g)
  h <- unname(kt$statistic)
  c(h, unname(kt$parameter), NA_real_, kt$p.value, h / (length(v) - 1))
}

.fp_test_feature <- function(f, d, glv, test) {
  sub <- d[d$feature == f, , drop = FALSE]
  by_group <- split(sub$value, factor(sub$group, levels = glv))
  if (any(lengths(by_group) < .fp_min_n)) return(rep(NA_real_, 5))
  res <- tryCatch(.fp_run_test(by_group, test), error = .fp_na_on_error)
  if (is.null(res)) rep(NA_real_, 5) else res
}

.fp_count_group <- function(g, d, f) sum(d$feature == f & d$group == g)

.fp_counts <- function(f, d, glv) vapply(glv, .fp_count_group, numeric(1), d = d, f = f)

.fp_qc_panel <- function(key, d) d$value[paste(d$feature, d$group, sep = ": ") == key]

.fp_wrap <- function(s, width = 12) scales::label_wrap(width)(s)

# A feature name in a panel title is wrapped at 24 characters, and a line that
# is still longer than 26 characters is cut and ends in an ellipsis.
.fp_cut_line <- function(s) ifelse(nchar(s) > 26, paste0(substr(s, 1, 25), "\u2026"), s)

.fp_title_one <- function(s) paste(.fp_cut_line(strwrap(s, width = 24)), collapse = "\n")

.fp_title_name <- function(s) vapply(s, .fp_title_one, character(1), USE.NAMES = FALSE)

recipe_feature_panel <- function(df, spec) {
  lg <- .fp_long(df, spec)
  d <- lg$d; feats <- lg$features; glv <- lg$groups
  if (length(feats) > .fp_max_features) {
    stop(sprintf("feature_panel draws at most %d features and found %d; use the heatmap recipe for a larger panel",
                 .fp_max_features, length(feats)))
  }
  if (length(glv) < 2) stop("feature_panel needs at least 2 groups")

  keys <- as.vector(t(outer(feats, glv, paste, sep = ": ")))
  sw <- vapply(split(d$value, factor(paste(d$feature, d$group, sep = ": "), levels = keys)),
               .fp_shapiro_normal, logical(1))
  all_normal <- all(sw)
  test <- .fp_resolve(spec$test$method, length(glv), all_normal)
  forced <- !identical(tolower(spec$test$method %||% "auto"), "auto")

  res <- vapply(feats, .fp_test_feature, numeric(5), d = d, glv = glv, test = test)
  res <- matrix(res, nrow = 5)
  p <- res[4, ]
  q <- stats::p.adjust(p, method = "BH")
  counts <- matrix(vapply(feats, .fp_counts, numeric(length(glv)), d = d, glv = glv),
                   nrow = length(glv))
  stats_df <- data.frame(feature = feats, test = .fp_labels[[test]], statistic = res[1, ],
                         df1 = res[2, ], df2 = res[3, ], p_value = p, p_adj = q,
                         p_adj_signif = .fp_signif(q), effect_size = res[5, ],
                         effect_size_name = .fp_effects[[test]], stringsAsFactors = FALSE)
  for (i in seq_along(glv)) stats_df[[paste0("n_", glv[i])]] <- as.integer(counts[i, ])
  n_untested <- sum(is.na(p))

  # Panels in the input order of the features, titled with the adjusted p.
  qtxt <- ifelse(is.na(q), "q = NA", paste0("q = ", formatC(q, format = "g", digits = 2)))
  titles <- paste0(.fp_title_name(feats), "\n", qtxt)
  d$panel <- factor(d$feature, levels = feats)
  d$group <- factor(d$group, levels = glv)
  ncol <- ceiling(sqrt(length(feats)))
  nrow <- ceiling(length(feats) / ncol)
  pal <- pub_palette(length(glv), spec$appearance$palette)
  n_sig <- sum(q < 0.05, na.rm = TRUE)
  sub <- sprintf("%s per feature; BH-adjusted q across %d features; %d with q < 0.05",
                 .fp_labels[[test]], length(feats), n_sig)
  y_lab <- spec$appearance$y_label %||% if (lg$layout == "long") spec$data$y else "Value"

  plot <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), width = 0.6, outlier.shape = NA,
                          alpha = 0.85, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.3, stroke = 0.25,
                        alpha = 0.8, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 1)) +
    ggplot2::facet_wrap(~ panel, ncol = ncol, scales = "free_y", axes = "all", axis.labels = "all_y",
                        labeller = ggplot2::as_labeller(stats::setNames(titles, feats))) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = .fp_wrap) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_extended(5)) +
    ggplot2::labs(x = spec$appearance$x_label %||% spec$data$x, y = y_lab,
                  title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10, lineheight = 1.05),
                   axis.text.x = ggplot2::element_text(size = 9, lineheight = 0.95),
                   axis.text.y = ggplot2::element_text(size = 9),
                   panel.spacing = ggplot2::unit(0.9, "lines"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  qc <- if (length(keys) <= .fp_max_qc_panels) {
    list(quantity = "group values of each feature",
         panels = stats::setNames(lapply(keys, .fp_qc_panel, d = d), keys))
  } else NULL

  test_meta <- list(
    name = .fp_labels[[test]], test = test, p_adjust = "BH", n_features = length(feats),
    n_groups = length(glv), groups = glv, n_significant = n_sig, n_untested = n_untested,
    forced = forced,
    assumptions = list(normality = list(test = "shapiro-wilk", all_normal = all_normal,
                                        n_checked = length(keys), n_non_normal = sum(!sw))))
  resolved <- list(family = "feature_panel", label = .fp_labels[[test]])
  panel_w <- max(2.2, 0.55 * length(glv) + 1.0)
  list(plot = plot, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .fp_methods(test_meta, stats_df, glv, y_lab, spec$data$x),
       build_script = list(fn = .fp_emit_script,
                           args = list(spec = spec, layout = lg$layout, test = test, ncol = ncol)),
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = d[, c("feature", "group", "value")], clean_steps = lg$steps, qc = qc,
       width = 0.9 + panel_w * ncol, height = 1.2 + 2.45 * nrow)
}

.fp_methods <- function(tm, stats_df, glv, y_lab, x) {
  v <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  choice <- if (isTRUE(tm$forced)) {
    sprintf("The %s was set for every feature by the user. ", tm$name)
  } else if (isTRUE(tm$assumptions$normality$all_normal)) {
    sprintf("One test family was chosen for the whole panel: the Shapiro-Wilk test was run on each of the %d feature-by-group samples, and because none departed from normality (P > 0.05) the parametric %s was used for every feature. ",
            tm$assumptions$normality$n_checked, tm$name)
  } else {
    sprintf("One test family was chosen for the whole panel: the Shapiro-Wilk test was run on each of the %d feature-by-group samples, and because %d of them departed from normality (P <= 0.05, or fewer than 3 values) the rank-based %s was used for every feature. ",
            tm$assumptions$normality$n_checked, tm$assumptions$normality$n_non_normal, tm$name)
  }
  nt <- paste(sprintf("%s (n = %s)", glv, vapply(paste0("n_", glv), .fp_n_range, character(1),
                                                 stats_df = stats_df)), collapse = ", ")
  paste0(
    sprintf("%s was compared between the levels of %s (%s) separately for each of %d features. ",
            cap_first(y_lab), x, paste(glv, collapse = ", "), tm$n_features),
    choice,
    sprintf("P values were adjusted across the %d features with the Benjamini-Hochberg method; %d features had an adjusted P (q) below 0.05. ",
            tm$n_features, tm$n_significant),
    sprintf("The effect size is %s. ", .fp_effects[[tm$test]]),
    if (tm$n_untested > 0) sprintf("%d feature(s) had a group with fewer than 2 values and were not tested. ", tm$n_untested) else "",
    sprintf("Group sizes per feature: %s. ", nt),
    sprintf("Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.fp_n_range <- function(col, stats_df) {
  r <- range(stats_df[[col]])
  if (r[1] == r[2]) as.character(r[1]) else sprintf("%d to %d", r[1], r[2])
}

.fp_emit_script <- function(spec, layout, test, ncol, in_name, fig_stub) {
  x <- spec$data$x
  read_lines <- if (layout == "long") {
    c(sprintf('d <- data.frame(feature = as.character(raw[["%s"]]), group = as.character(raw[["%s"]]), value = as.numeric(raw[["%s"]]))',
              spec$data$facet, x, spec$data$y))
  } else {
    c(sprintf('stack_col <- function(cl) data.frame(feature = cl, group = as.character(raw[["%s"]]), value = raw[[cl]])', x),
      sprintf('cols <- setdiff(names(raw)[vapply(raw, is.numeric, logical(1))], "%s")', x),
      "d <- do.call(rbind, lapply(cols, stack_col))")
  }
  test_call <- switch(test,
    welch_t  = "t.test(value ~ group, data = d[d$feature == f, ])$p.value",
    wilcoxon = "wilcox.test(value ~ group, data = d[d$feature == f, ])$p.value",
    anova    = "oneway.test(value ~ group, data = d[d$feature == f, ], var.equal = TRUE)$p.value",
    kruskal  = "kruskal.test(value ~ group, data = d[d$feature == f, ])$p.value")
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    read_lines,
    "d <- d[is.finite(d$value) & !is.na(d$group), ]",
    "feats <- unique(d$feature)",
    sprintf("test_one <- function(f) %s", test_call),
    "p <- vapply(feats, test_one, numeric(1))",
    'q <- p.adjust(p, method = "BH")',
    'titles <- paste0(feats, "\\nq = ", formatC(q, format = "g", digits = 2))',
    "d$panel <- factor(titles[match(d$feature, feats)], levels = titles)",
    "g <- ggplot(d, aes(group, value, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.85) +",
    "  geom_jitter(shape = 21, width = 0.12, height = 0) +",
    sprintf("  facet_wrap(~ panel, ncol = %d, scales = \"free_y\") +", ncol),
    '  theme_classic(base_size = 13) + theme(legend.position = "none", strip.background = element_blank())',
    sprintf('ggsave("%s.pdf", g, width = 9, height = 8)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
