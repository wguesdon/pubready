# Recipe: compare two groups on one continuous outcome.
# Chooses t-test or Mann-Whitney from the assumption checks (unless the caller
# forces a test), draws a box or violin plot, and puts the significance bracket
# on top. Returns the plot, a tidy stats table, and test metadata.

recipe_two_group_compare <- function(df, spec) {
  x <- spec$data$x
  y <- spec$data$y

  df[[x]] <- factor(df[[x]])
  lvls <- levels(df[[x]])
  if (length(lvls) != 2) {
    stop(sprintf("two_group_compare needs exactly 2 groups in '%s'; found %d (%s)",
                 x, length(lvls), paste(lvls, collapse = ", ")))
  }

  paired <- isTRUE(spec$test$paired)

  # Assumption checks.
  by_group <- split(df[[y]], df[[x]])
  norm <- lapply(by_group, shapiro_safe)
  all_normal <- all(vapply(norm, function(n) isTRUE(n$normal), logical(1)))
  var_p <- tryCatch(
    stats::var.test(df[[y]] ~ df[[x]])$p.value,
    error = function(e) NA_real_
  )
  equal_var <- is.na(var_p) || var_p > 0.05

  resolved <- resolve_two_group_method(spec$test$method, all_normal, equal_var, paired)

  fml <- stats::as.formula(sprintf("`%s` ~ `%s`", y, x))
  if (resolved$family == "t") {
    stat <- rstatix::t_test(df, fml, paired = paired, var.equal = isTRUE(resolved$var_equal))
    eff  <- tryCatch(
      rstatix::cohens_d(df, fml, paired = paired),
      error = function(e) NULL
    )
    eff_name <- "Cohen's d"
  } else {
    stat <- rstatix::wilcox_test(df, fml, paired = paired)
    eff  <- tryCatch(rstatix::wilcox_effsize(df, fml, paired = paired),
                     error = function(e) NULL)
    eff_name <- "r (rank-biserial)"
  }
  eff_val <- if (!is.null(eff) && "effsize" %in% names(eff)) eff$effsize[1] else NA_real_

  stat <- rstatix::add_significance(stat)
  stat <- tryCatch(
    rstatix::add_xy_position(stat, x = x),
    error = function(e) { stat$y.position <- max(df[[y]], na.rm = TRUE) * 1.08; stat }
  )
  if (is.null(stat$y.position)) stat$y.position <- max(df[[y]], na.rm = TRUE) * 1.08

  bracket_label <- spec$appearance$bracket$label %||% "p.signif"
  if (bracket_label == "p.format") {
    stat$plabel <- paste0("p = ", formatC(stat$p, format = "g", digits = 2))
    label_col <- "plabel"
  } else {
    label_col <- "p.signif"
  }

  pal  <- spec$appearance$palette %||% pubplot_palette[seq_len(2)]
  geom <- spec$appearance$geom %||% "box"

  # fill is mapped inside each geom, not globally: the significance-bracket
  # layer inherits global aesthetics, and a global fill would break it because
  # its data has no grouping column.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  p <- if (geom == "violin") {
    p + ggplot2::geom_violin(ggplot2::aes(fill = .data[[x]]),
                             trim = FALSE, width = 0.7, alpha = 0.9)
  } else {
    p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[x]]),
                              width = 0.6, outlier.shape = NA, alpha = 0.9)
  }
  if (isTRUE(spec$appearance$show_points)) {
    p <- p + ggplot2::geom_jitter(ggplot2::aes(fill = .data[[x]]),
                                  width = 0.12, size = 1.6, alpha = 0.75,
                                  shape = 21, stroke = 0.3)
  }
  p <- p +
    ggplot2::scale_fill_manual(values = pal) +
    ggpubr::stat_pvalue_manual(stat, label = label_col, tip.length = 0.01,
                               bracket.size = 0.4) +
    ggplot2::labs(
      x     = spec$appearance$x_label %||% x,
      y     = spec$appearance$y_label %||% y,
      title = spec$appearance$title
    ) +
    theme_pubplot()

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

  test_meta <- list(
    name        = resolved$label,
    statistic   = if (!is.null(stat$statistic)) unname(stat$statistic[1]) else NULL,
    p_value     = stat$p[1],
    adjustment  = spec$test$p_adjust,
    effect_size = list(name = eff_name, value = eff_val),
    assumptions = list(
      normality      = list(test = "shapiro-wilk",
                            p_group1 = norm[[g1]]$p, p_group2 = norm[[g2]]$p),
      equal_variance = list(test = "F test", p = var_p)
    )
  )

  methods <- methods_paragraph(resolved, stats_df, test_meta, spec)
  build_script <- function(in_name, fig_stub) emit_script(spec, resolved, in_name, fig_stub)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script)
}
