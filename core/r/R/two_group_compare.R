# Recipe: compare two groups on one continuous outcome.
# Chooses t-test or Mann-Whitney from the assumption checks (unless the caller
# forces a test), draws a box or violin plot, and puts the significance bracket
# on top. Returns the plot, a tidy stats table, and test metadata.

# The --test values that two_group_compare knows (F63).
.two_group_compare_tests <- c("auto", "t", "t_test", "student_t", "welch_t", "wilcoxon",
                              "mann_whitney", "wilcox", "mwu")

recipe_two_group_compare <- function(df, spec) {
  check_test_method("two_group_compare", spec$test$method, .two_group_compare_tests)
  cl <- clean_tidy(df, spec)
  df <- cl$df
  x <- spec$data$x
  y <- spec$data$y

  # One level rule in both engines: the --levels order, else numeric order when
  # every label is a number, else code-point order. The first level is group 1.
  lab <- level_labels(df[[x]])
  lvls <- level_order(lab, unlist(spec$data$levels), "group")
  df[[x]] <- factor(lab, levels = lvls)
  if (length(lvls) != 2) {
    stop(sprintf("two_group_compare needs exactly 2 groups in '%s'; found %d (%s)",
                 x, length(lvls), paste(lvls, collapse = ", ")))
  }
  check_spread(df[[y]], "two_group_compare", sprintf("the outcome '%s'", y))

  paired <- isTRUE(spec$test$paired)

  # Assumption checks. Equal variance uses Levene's test, as the multi-group
  # recipe and both Python recipes do, so the column compares across engines.
  # The tests run on a copy with syntactic column names (.x, .y), so a column
  # name with a space or a symbol cannot break a formula. The labels and the
  # tables keep the original names.
  sdf <- data.frame(.x = df[[x]], .y = df[[y]])
  fml <- .y ~ .x
  by_group <- split(df[[y]], df[[x]])
  norm <- lapply(by_group, shapiro_safe)
  all_normal <- all(vapply(norm, is_normal_check, logical(1)))
  # Under --paired the rows pair by their order in each group, and the test is
  # chosen from the normality of the within-pair differences (group1 minus
  # group2), the quantity that the paired tests use.
  diffs <- NULL
  norm_diff <- list(p = NA_real_)
  if (paired) {
    diffs <- .tg_paired_diffs(by_group, lvls, x)
    norm_diff <- shapiro_safe(diffs)
    all_normal <- isTRUE(norm_diff$normal)
  }
  var_p <- tryCatch(
    rstatix::levene_test(sdf, fml)$p[1],
    error = na_on_error
  )
  equal_var <- is.na(var_p) || var_p > 0.05

  rank_min_p <- .tg_rank_min_p(by_group, lvls, paired)
  resolved <- resolve_two_group_method(spec$test$method, all_normal, equal_var, paired,
                                       rank_min_p, .tg_t_can_run(by_group, lvls, paired))
  .tg_check_variance(resolved, by_group, lvls)

  if (resolved$family == "t") {
    stat <- rstatix::t_test(sdf, fml, paired = paired, var.equal = isTRUE(resolved$var_equal))
    # Cohen's d over the pooled SD (paired: dz), as the Python engine and
    # feature_panel compute it; rstatix::cohens_d uses the average variance.
    eff_name <- if (paired) "Cohen's dz" else "Cohen's d"
    eff_val <- cohens_d_signed(by_group[[lvls[1]]], by_group[[lvls[2]]], paired)
  } else {
    stat <- rstatix::wilcox_test(sdf, fml, paired = paired)
    eff_name <- "rank-biserial r"
    eff_val <- if (paired) {
      rank_biserial_paired(by_group[[lvls[1]]] - by_group[[lvls[2]]])
    } else {
      rank_biserial_unpaired(stat$statistic[1],
                             length(by_group[[lvls[1]]]),
                             length(by_group[[lvls[2]]]))
    }
  }

  # rstatix rounds the p-value to three significant digits. The table keeps the
  # full-precision value, so it matches the Python engine digit for digit. The
  # vectors are passed rather than the formula, because the formula method of
  # stats::t.test and stats::wilcox.test rejects paired = TRUE.
  p_full <- tryCatch(
    {
      g1v <- by_group[[lvls[1]]]
      g2v <- by_group[[lvls[2]]]
      if (resolved$family == "t") {
        stats::t.test(g1v, g2v, paired = paired,
                      var.equal = isTRUE(resolved$var_equal))$p.value
      } else {
        suppressWarnings(stats::wilcox.test(g1v, g2v, paired = paired)$p.value)
      }
    },
    error = na_on_error
  )
  if (!is.na(p_full)) stat$p <- p_full

  stat <- rstatix::add_significance(stat)
  xy <- tryCatch(
    rstatix::add_xy_position(stat, x = ".x"),
    error = null_on_error
  )
  if (!is.null(xy)) stat <- xy
  if (is.null(stat$y.position)) stat$y.position <- max(df[[y]], na.rm = TRUE) * 1.08

  bracket_label <- spec$appearance$bracket$label %||% "p.signif"
  if (bracket_label == "p.format") {
    stat$plabel <- paste0("p = ", formatC(stat$p, format = "g", digits = 2))
    label_col <- "plabel"
  } else {
    label_col <- "p.signif"
  }

  pal  <- spec$appearance$palette %||% pubready_palette[seq_len(2)]
  geom <- spec$appearance$geom %||% "box"

  # A violin with trim = FALSE reaches past the data, so its bracket starts
  # above the drawn violin, not above the largest value.
  if (identical(geom, "violin")) {
    tops <- vapply(by_group, geom_top, numeric(1), geom = "violin")
    span <- max(tops) - min(df[[y]], na.rm = TRUE)
    stat$y.position <- max(tops) + 0.06 * span
  }

  # fill is mapped inside each geom, not globally: the significance-bracket
  # layer inherits global aesthetics, and a global fill would break it because
  # its data has no grouping column.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  p <- add_dist_geom(p, x, geom, spec$appearance$show_points)
  p <- p +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = stats::setNames(wrap_labels(lvls), lvls)) +
    ggplot2::scale_y_continuous(labels = plain_number_labels) +
    ggpubr::stat_pvalue_manual(stat, label = label_col, tip.length = 0.01,
                               bracket.size = 0.4) +
    ggplot2::labs(
      x     = spec$appearance$x_label %||% x,
      y     = spec$appearance$y_label %||% y,
      title = spec$appearance$title
    ) +
    pub_theme(spec$appearance$theme)

  ylim <- spec$appearance$y_limits
  ymin <- if (length(ylim) >= 1) ylim[[1]] else NULL
  ymax <- if (length(ylim) >= 2) ylim[[2]] else NULL
  if (!is.null(ymin) || !is.null(ymax)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax))
  }

  g1 <- as.character(stat$group1[1]); g2 <- as.character(stat$group2[1])
  stats_df <- data.frame(
    outcome          = y,
    group1           = g1,
    group2           = g2,
    contrast         = paste(g1, "-", g2),
    n1               = if (!is.null(stat$n1)) stat$n1[1] else length(by_group[[g1]]),
    n2               = if (!is.null(stat$n2)) stat$n2[1] else length(by_group[[g2]]),
    test             = resolved$label,
    statistic        = if (!is.null(stat$statistic)) unname(stat$statistic[1]) else NA_real_,
    df               = if ("df" %in% names(stat)) unname(stat$df[1]) else NA_real_,
    p_value          = stat$p[1],
    p_adjust_method  = spec$test$p_adjust,
    effect_size_name = eff_name,
    effect_size      = eff_val,
    shapiro_p_group1 = norm[[g1]]$p,
    shapiro_p_group2 = norm[[g2]]$p,
    var_equal_p      = var_p,
    stringsAsFactors = FALSE
  )
  stats_df$shapiro_p_diff <- norm_diff$p

  test_meta <- list(
    name        = resolved$label,
    statistic   = if (!is.null(stat$statistic)) unname(stat$statistic[1]) else NULL,
    p_value     = stat$p[1],
    adjustment  = spec$test$p_adjust,
    effect_size = list(name = eff_name, value = eff_val),
    assumptions = list(
      normality      = list(test = "shapiro-wilk",
                            p_group1 = norm[[g1]]$p, p_group2 = norm[[g2]]$p),
      equal_variance = list(test = "Levene", p = var_p)
    )
  )
  test_meta$rank_test_attainable_p <- rank_min_p
  if (paired) {
    test_meta$assumptions$normality_of_differences <- list(test = "shapiro-wilk", p = norm_diff$p)
  }

  methods <- .tg_methods(resolved, stats_df, test_meta, spec)
  build_script <- list(fn = .tg_emit_script,
                       args = list(spec = spec, resolved = resolved, lvls = lvls))

  qc <- if (paired) {
    list(quantity = "within-pair differences", panels = list(differences = diffs))
  } else {
    list(quantity = "group values", panels = by_group)
  }
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       df_used = df, clean_steps = cl$steps, qc = qc)
}

# The sentence that defines a Cohen's d of the methods text. The Python engine
# writes the same words.
.tg_effect_definition <- function(name) {
  if (identical(name, "Cohen's d")) {
    " Cohen's d is the mean difference over the pooled standard deviation."
  } else if (identical(name, "Cohen's dz")) {
    " Cohen's dz is the mean of the within-pair differences over their standard deviation."
  } else ""
}

# The within-pair differences, group1 minus group2, with the rows paired by
# their order in each group. Groups of unequal size cannot pair.
.tg_paired_diffs <- function(by_group, lvls, x) {
  a <- by_group[[lvls[1]]]; b <- by_group[[lvls[2]]]
  if (length(a) != length(b)) {
    stop(sprintf("two_group_compare --paired pairs the rows by their order in each group of '%s', so the groups need the same size; '%s' has %d and '%s' has %d",
                 x, lvls[1], length(a), lvls[2], length(b)))
  }
  a - b
}

# A t-test divides by a variance: the two-sample test by the variance of both
# groups, the paired test by the variance of the differences. When that variance
# is zero the statistic is infinite, so the run stops and names the samples. The
# Python engine has the same rule.
.tg_check_variance <- function(resolved, by_group, lvls) {
  paired <- isTRUE(resolved$paired)
  if (resolved$family != "t" || .tg_t_can_run(by_group, lvls, paired)) return(invisible(NULL))
  subject <- if (paired) "the within-pair differences have" else zero_variance_subject(lvls)
  stop(zero_variance_message(resolved$label, subject, "use --test wilcoxon"))
}

# Whether the t-test can run: the variance that it divides by is not zero.
.tg_t_can_run <- function(by_group, lvls, paired) {
  a <- by_group[[lvls[1]]]; b <- by_group[[lvls[2]]]
  if (paired) return(!(length(a) == length(b) && is_zero_variance(a - b)))
  length(zero_variance_groups(by_group)) < 2
}

# The smallest P that the rank test can reach at the group sizes (F61).
.tg_rank_min_p <- function(by_group, lvls, paired) {
  n <- lengths(by_group[lvls])
  if (paired) rank_test_min_p(min(n), "signed_rank") else rank_test_min_p(n, "mann_whitney")
}

# The methods paragraph. The Python engine writes the same sentences
# (two_group_compare._methods).
.tg_methods <- function(resolved, stats_df, test_meta, spec) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v <- pkg_versions(c("rstatix", "ggpubr", "ggprism"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  norm_txt <- if (isTRUE(resolved$paired)) {
    sprintf("The rows were paired by their order in each group, and the normality of the within-pair differences (%s minus %s) was assessed with the Shapiro-Wilk test.",
            stats_df$group1[1], stats_df$group2[1])
  } else "Normality was assessed with the Shapiro-Wilk test."
  # A paired t-test works on the within-pair differences, so no variance
  # sentence applies to it.
  var_txt <- if (resolved$family == "t") {
    if (isTRUE(resolved$paired)) ""
    else if (isTRUE(resolved$var_equal)) " Equal variances were assumed."
    else " The Welch t-test does not assume equal variances; Levene's test is reported in the stats table."
  } else " A non-parametric test was used because the normality assumption was not met."
  if (tolower(spec$test$method %||% "auto") != "auto") {
    var_txt <- paste0(if (resolved$family == "t") var_txt else "", " The test was set by the user.")
  }
  if (isTRUE(resolved$kept_parametric) || resolved$family == "wilcoxon") {
    floor_txt <- rank_floor_text(resolved$rank_min_p, resolved$kept_parametric)
    if (nzchar(floor_txt)) var_txt <- paste0(var_txt, " ", floor_txt)
  }
  eff_txt <- if (!is.null(test_meta$effect_size)) {
    # The value is given, not only the name, so the paragraph can be read on its
    # own. The Python engine writes the same sentence.
    eff_v <- test_meta$effect_size$value
    if (is.null(eff_v) || !is.finite(eff_v)) {
      sprintf(" Effect size is reported as %s (%s minus %s).", test_meta$effect_size$name,
              stats_df$group1[1], stats_df$group2[1])
    } else {
      sprintf(" Effect size is reported as %s = %s (%s minus %s).", test_meta$effect_size$name,
              formatC(eff_v, format = "f", digits = 2), stats_df$group1[1], stats_df$group2[1])
    }
  } else ""
  eff_txt <- paste0(eff_txt, .tg_effect_definition(test_meta$effect_size$name))
  adj_txt <- if (!identical(spec$test$p_adjust, "none")) {
    sprintf(" P-values were adjusted for multiple comparisons using the %s method.", spec$test$p_adjust)
  } else ""
  n_txt <- if (all(c("n1", "n2") %in% names(stats_df))) {
    sprintf(" Group sizes were n = %d and n = %d.", stats_df$n1[1], stats_df$n2[1])
  } else ""

  paste0(
    sprintf("%s was compared between the two groups with a %s.",
            cap_first(y_label), resolved$label),
    " ", norm_txt, var_txt, eff_txt, adj_txt,
    " Significance was set at P < 0.05.",
    sprintf(" Analyses were performed in R %s with rstatix %s; figures were produced with ggpubr %s and ggprism %s.",
            rver, v$rstatix, v$ggpubr, v$ggprism),
    n_txt
  )
}

# A standalone, readable script that reproduces the figure from the bundled
# input using the concrete resolved parameters.
.tg_emit_script <- function(spec, resolved, lvls, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  pal <- spec$appearance$palette %||% pubready_palette[1:2]
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  paired <- if (isTRUE(resolved$paired)) "TRUE" else "FALSE"
  geom <- spec$appearance$geom %||% "box"
  bracket_label <- spec$appearance$bracket$label %||% "p.signif"

  if (identical(bracket_label, "p.format")) {
    plabel_line <- 'stat$plabel <- paste0("p = ", formatC(stat$p, format = "g", digits = 2))'
    label_arg <- "plabel"
  } else {
    plabel_line <- character(0)
    label_arg <- "p.signif"
  }

  test_line <- if (resolved$family == "t") {
    sprintf('stat <- rstatix::t_test(sdf, .y ~ .x, paired = %s, var.equal = %s)',
            paired, if (isTRUE(resolved$var_equal)) "TRUE" else "FALSE")
  } else {
    sprintf('stat <- rstatix::wilcox_test(sdf, .y ~ .x, paired = %s)', paired)
  }
  gl <- script_geom_lines(x, geom, spec$appearance$show_points %||% TRUE)

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction of this figure. Run inside the pinned pubready",
    "# container (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(ggpubr); library(rstatix); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf('df[["%s"]] <- factor(df[["%s"]], levels = c(%s))', x, x,
            paste(sprintf('"%s"', lvls), collapse = ", ")),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "# The test runs on a copy with syntactic column names.",
    sprintf('sdf <- data.frame(.x = df[["%s"]], .y = df[["%s"]])', x, y),
    "",
    test_line,
    "stat <- rstatix::add_significance(stat)",
    'stat <- rstatix::add_xy_position(stat, x = ".x")',
    plabel_line,
    "",
    gl$prep,
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]])) +', x, y),
    gl$layers,
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  stat_pvalue_manual(stat, label = "%s", tip.length = 0.01, bracket.size = 0.4) +',
            label_arg),
    sprintf('  labs(x = "%s", y = "%s") +',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    "  theme_classic(base_size = 13) +",
    '  theme(legend.position = "none")',
    "",
    sprintf('ggsave("%s.pdf", p, width = 3.8, height = 4.0)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
