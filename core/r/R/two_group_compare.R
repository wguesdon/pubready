# Recipe: compare two groups on one continuous outcome.
# Chooses t-test or Mann-Whitney from the assumption checks (unless the caller
# forces a test), draws a box or violin plot, and puts the significance bracket
# on top. Returns the plot, a tidy stats table, and test metadata.

recipe_two_group_compare <- function(df, spec) {
  cl <- clean_tidy(df, spec)
  df <- cl$df
  x <- spec$data$x
  y <- spec$data$y

  df[[x]] <- factor(df[[x]])
  lvls <- levels(df[[x]])
  if (length(lvls) != 2) {
    stop(sprintf("two_group_compare needs exactly 2 groups in '%s'; found %d (%s)",
                 x, length(lvls), paste(lvls, collapse = ", ")))
  }

  paired <- isTRUE(spec$test$paired)

  # Assumption checks. Equal variance uses Levene's test, as the multi-group
  # recipe and both Python recipes do, so the column compares across engines.
  fml <- stats::as.formula(sprintf("`%s` ~ `%s`", y, x))
  by_group <- split(df[[y]], df[[x]])
  norm <- lapply(by_group, shapiro_safe)
  all_normal <- all(vapply(norm, is_normal_check, logical(1)))
  var_p <- tryCatch(
    rstatix::levene_test(df, fml)$p[1],
    error = na_on_error
  )
  equal_var <- is.na(var_p) || var_p > 0.05

  resolved <- resolve_two_group_method(spec$test$method, all_normal, equal_var, paired)
  .tg_check_variance(resolved, by_group, lvls)

  if (resolved$family == "t") {
    stat <- rstatix::t_test(df, fml, paired = paired, var.equal = isTRUE(resolved$var_equal))
    eff  <- tryCatch(
      rstatix::cohens_d(df, fml, paired = paired),
      error = null_on_error
    )
    eff_name <- "Cohen's d"
    eff_val <- if (!is.null(eff) && "effsize" %in% names(eff)) eff$effsize[1] else NA_real_
  } else {
    stat <- rstatix::wilcox_test(df, fml, paired = paired)
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
    rstatix::add_xy_position(stat, x = x),
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

  methods <- methods_paragraph(resolved, stats_df, test_meta, spec)
  build_script <- list(fn = emit_script, args = list(spec = spec, resolved = resolved))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       df_used = df, clean_steps = cl$steps,
       qc = list(quantity = "group values", panels = by_group))
}

# A t-test divides by a variance: the two-sample test by the variance of both
# groups, the paired test by the variance of the differences. When that variance
# is zero the statistic is infinite, so the run stops and names the samples. The
# Python engine has the same rule.
.tg_check_variance <- function(resolved, by_group, lvls) {
  if (resolved$family != "t") return(invisible(NULL))
  a <- by_group[[lvls[1]]]; b <- by_group[[lvls[2]]]
  if (isTRUE(resolved$paired)) {
    if (length(a) == length(b) && is_zero_variance(a - b)) {
      stop(zero_variance_message(resolved$label, "the within-pair differences have",
                                 "use --test wilcoxon"))
    }
  } else if (length(zero_variance_groups(by_group)) == 2) {
    stop(zero_variance_message(resolved$label, zero_variance_subject(lvls), "use --test wilcoxon"))
  }
  invisible(NULL)
}
