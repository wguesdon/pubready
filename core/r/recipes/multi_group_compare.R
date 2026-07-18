# Recipe: compare a continuous outcome across 3+ groups (also works for 2).
# Picks one-way ANOVA, Welch's ANOVA, or Kruskal-Wallis from the assumption
# checks, runs the matching post-hoc test, and brackets only the comparisons
# that reach significance (adjusted P < 0.05). The full pairwise table is saved.

resolve_multi_group_method <- function(method, all_normal, equal_var) {
  method <- tolower(method %||% "auto")
  fam <- if (method == "auto") {
    if (all_normal && equal_var) "anova"
    else if (all_normal && !equal_var) "welch"
    else "kruskal"
  } else if (method %in% c("anova", "aov")) "anova"
  else if (method %in% c("welch", "welch_anova")) "welch"
  else if (method %in% c("kruskal", "kruskal_wallis", "kw")) "kruskal"
  else stop(sprintf("unknown ANOVA method '%s'", method))

  list(
    family        = fam,
    label         = switch(fam, anova = "one-way ANOVA",
                                 welch = "Welch's one-way ANOVA",
                                 kruskal = "Kruskal-Wallis test"),
    article       = switch(fam, anova = "a ", welch = "", kruskal = "a "),
    posthoc       = switch(fam, anova = "tukey", welch = "games_howell", kruskal = "dunn"),
    posthoc_label = switch(fam, anova = "Tukey's HSD",
                                 welch = "the Games-Howell test",
                                 kruskal = "Dunn's test with Benjamini-Hochberg adjustment")
  )
}

recipe_multi_group_compare <- function(df, spec) {
  x <- spec$data$x
  y <- spec$data$y

  df[[x]] <- factor(df[[x]])
  lvls <- levels(df[[x]])
  if (length(lvls) < 2) {
    stop(sprintf("multi_group_compare needs 2 or more groups in '%s'; found %d", x, length(lvls)))
  }

  fml <- stats::as.formula(sprintf("`%s` ~ `%s`", y, x))
  by_group <- split(df[[y]], df[[x]])
  norm <- lapply(by_group, shapiro_safe)
  all_normal <- all(vapply(norm, function(n) isTRUE(n$normal), logical(1)))
  lev_p <- tryCatch(rstatix::levene_test(df, fml)$p[1], error = function(e) NA_real_)
  equal_var <- is.na(lev_p) || lev_p > 0.05

  resolved <- resolve_multi_group_method(spec$test$method, all_normal, equal_var)

  # Omnibus test and matching post-hoc.
  if (resolved$family == "anova") {
    at <- as.data.frame(rstatix::anova_test(df, fml))
    omnibus <- list(statistic = at$F[1], df1 = at$DFn[1], df2 = at$DFd[1], p = at$p[1],
                    effect_size = list(name = "generalized eta-squared", value = at$ges[1]))
    ph <- as.data.frame(rstatix::tukey_hsd(df, fml))
  } else if (resolved$family == "welch") {
    wt <- as.data.frame(rstatix::welch_anova_test(df, fml))
    omnibus <- list(statistic = wt$statistic[1], df1 = wt$DFn[1], df2 = wt$DFd[1], p = wt$p[1],
                    effect_size = NULL)
    ph <- as.data.frame(rstatix::games_howell_test(df, fml))
  } else {
    kt <- as.data.frame(rstatix::kruskal_test(df, fml))
    ke <- tryCatch(as.data.frame(rstatix::kruskal_effsize(df, fml)), error = function(e) NULL)
    omnibus <- list(statistic = kt$statistic[1], df1 = kt$df[1], df2 = NA_real_, p = kt$p[1],
                    effect_size = if (!is.null(ke)) list(name = "eta-squared (H)", value = ke$effsize[1]) else NULL)
    ph <- as.data.frame(rstatix::dunn_test(df, fml, p.adjust.method = "BH"))
  }

  ph$group1 <- as.character(ph$group1)
  ph$group2 <- as.character(ph$group2)
  estimate <- if ("estimate" %in% names(ph)) ph$estimate else NA_real_

  stats_df <- data.frame(
    outcome      = y,
    group1       = ph$group1,
    group2       = ph$group2,
    n1           = vapply(ph$group1, function(g) length(by_group[[g]]), integer(1)),
    n2           = vapply(ph$group2, function(g) length(by_group[[g]]), integer(1)),
    post_hoc     = resolved$posthoc_label,
    estimate     = estimate,
    p_adj        = ph$p.adj,
    p_adj_signif = ph$p.adj.signif,
    displayed    = ph$p.adj < 0.05,
    stringsAsFactors = FALSE
  )

  # Significant comparisons only, stacked above the data.
  sig <- ph[which(ph$p.adj < 0.05), , drop = FALSE]
  bracket_label <- spec$appearance$bracket$label %||% "p.signif"
  if (nrow(sig) > 0) {
    ymax <- max(df[[y]], na.rm = TRUE)
    rng  <- diff(range(df[[y]], na.rm = TRUE))
    if (!is.finite(rng) || rng == 0) rng <- max(abs(ymax), 1)
    sig$y.position <- ymax + rng * (0.07 + 0.09 * (seq_len(nrow(sig)) - 1))
    sig$lab <- if (identical(bracket_label, "p.format")) {
      paste0("p = ", formatC(sig$p.adj, format = "g", digits = 2))
    } else {
      sig$p.adj.signif
    }
  }

  pal  <- spec$appearance$palette %||% pubplot_palette
  pal  <- rep(pal, length.out = length(lvls))
  geom <- spec$appearance$geom %||% "box"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  p <- if (geom == "violin") {
    p + ggplot2::geom_violin(ggplot2::aes(fill = .data[[x]]), trim = FALSE, width = 0.7, alpha = 0.9)
  } else {
    p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[x]]), width = 0.6, outlier.shape = NA, alpha = 0.9)
  }
  if (isTRUE(spec$appearance$show_points)) {
    p <- p + ggplot2::geom_jitter(ggplot2::aes(fill = .data[[x]]),
                                  width = 0.12, size = 1.4, alpha = 0.7, shape = 21, stroke = 0.3)
  }
  omni_lab <- sprintf("%s, P = %s", cap_first(resolved$label),
                      formatC(omnibus$p, format = "g", digits = 2))
  p <- p +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% y,
                  title = spec$appearance$title,
                  subtitle = omni_lab) +
    theme_pubplot() +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))
  if (nrow(sig) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(sig, label = "lab", tip.length = 0.01, bracket.size = 0.4)
  }

  ylim <- spec$appearance$y_limits
  ymin <- if (length(ylim) >= 1) ylim[[1]] else NULL
  ymax2 <- if (length(ylim) >= 2) ylim[[2]] else NULL
  if (!is.null(ymin) || !is.null(ymax2)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax2))
  }

  test_meta <- list(
    name = resolved$label,
    omnibus = list(statistic = omnibus$statistic, df1 = omnibus$df1, df2 = omnibus$df2,
                   p_value = omnibus$p, effect_size = omnibus$effect_size),
    post_hoc = resolved$posthoc_label,
    n_comparisons = nrow(ph),
    n_significant = sum(ph$p.adj < 0.05, na.rm = TRUE),
    assumptions = list(
      normality      = list(test = "shapiro-wilk", all_normal = all_normal),
      equal_variance = list(test = "Levene", p = lev_p)
    )
  )

  n_by_group <- vapply(by_group, length, integer(1))
  methods <- .mg_methods(resolved, omnibus, spec, n_by_group)
  build_script <- function(in_name, fig_stub) {
    .mg_emit_script(spec, resolved, pal[seq_len(length(lvls))], in_name, fig_stub)
  }

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script)
}

.mg_methods <- function(resolved, omnibus, spec, n_by_group) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v    <- pkg_versions(c("rstatix", "ggpubr", "ggprism"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  pfmt <- function(p) formatC(p, format = "g", digits = 2)

  omni <- if (resolved$family == "kruskal") {
    sprintf("a Kruskal-Wallis test (H(%s) = %.2f, P = %s)",
            omnibus$df1, omnibus$statistic, pfmt(omnibus$p))
  } else {
    sprintf("%s%s (F(%s, %s) = %.2f, P = %s)",
            resolved$article, resolved$label,
            omnibus$df1, omnibus$df2, omnibus$statistic, pfmt(omnibus$p))
  }
  eff <- if (!is.null(omnibus$effect_size)) {
    sprintf(" The effect size was %s = %.2f.", omnibus$effect_size$name, omnibus$effect_size$value)
  } else ""
  ph_txt <- sprintf(" Pairwise comparisons used %s; only comparisons reaching significance (adjusted P < 0.05) are annotated on the figure.",
                    resolved$posthoc_label)
  ns <- paste(sprintf("%s (n = %d)", names(n_by_group), n_by_group), collapse = ", ")

  paste0(
    sprintf("%s was compared across %d groups with %s.",
            cap_first(y_label), length(n_by_group), omni),
    eff, ph_txt,
    sprintf(" Analyses were performed in R %s with rstatix %s; figures were produced with ggpubr %s and ggprism %s.",
            rver, v$rstatix, v$ggpubr, v$ggprism),
    sprintf(" Group sizes: %s.", ns)
  )
}

.mg_emit_script <- function(spec, resolved, pal, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  geom <- spec$appearance$geom %||% "box"
  geom_line <- if (geom == "violin") {
    sprintf('geom_violin(aes(fill = .data[["%s"]]), trim = FALSE, width = 0.7, alpha = 0.9)', x)
  } else {
    sprintf('geom_boxplot(aes(fill = .data[["%s"]]), width = 0.6, outlier.shape = NA, alpha = 0.9)', x)
  }
  posthoc_line <- switch(resolved$posthoc,
    tukey        = sprintf('ph <- rstatix::tukey_hsd(df, %s ~ %s)', y, x),
    games_howell = sprintf('ph <- rstatix::games_howell_test(df, %s ~ %s)', y, x),
    dunn         = sprintf('ph <- rstatix::dunn_test(df, %s ~ %s, p.adjust.method = "BH")', y, x)
  )
  omni_line <- switch(resolved$family,
    anova   = sprintf('omni <- as.data.frame(rstatix::anova_test(df, %s ~ %s))', y, x),
    welch   = sprintf('omni <- as.data.frame(rstatix::welch_anova_test(df, %s ~ %s))', y, x),
    kruskal = sprintf('omni <- as.data.frame(rstatix::kruskal_test(df, %s ~ %s))', y, x)
  )
  omni_name <- cap_first(resolved$label)
  lab_line <- if (identical(spec$appearance$bracket$label %||% "p.signif", "p.format")) {
    'sig$lab <- paste0("p = ", formatC(sig$p.adj, format = "g", digits = 2))'
  } else {
    "sig$lab <- sig$p.adj.signif"
  }

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(ggpubr); library(rstatix); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf('df[["%s"]] <- factor(df[["%s"]])', x, x),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "",
    posthoc_line,
    "ph$group1 <- as.character(ph$group1); ph$group2 <- as.character(ph$group2)",
    "sig <- ph[which(ph$p.adj < 0.05), , drop = FALSE]",
    "if (nrow(sig) > 0) {",
    sprintf('  ymax <- max(df[["%s"]], na.rm = TRUE); rng <- diff(range(df[["%s"]], na.rm = TRUE))', y, y),
    "  if (!is.finite(rng) || rng == 0) rng <- max(abs(ymax), 1)",
    "  sig$y.position <- ymax + rng * (0.07 + 0.09 * (seq_len(nrow(sig)) - 1))",
    sprintf("  %s", lab_line),
    "}",
    "",
    omni_line,
    sprintf('subtitle_txt <- paste0("%s, P = ", formatC(omni$p[1], format = "g", digits = 2))', omni_name),
    "",
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]])) +', x, y),
    sprintf("  %s +", geom_line),
    sprintf('  geom_jitter(aes(fill = .data[["%s"]]), width = 0.12, size = 1.4, alpha = 0.7, shape = 21, stroke = 0.3) +', x),
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  labs(x = "%s", y = "%s", subtitle = subtitle_txt) +',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    "  theme_classic(base_size = 13) +",
    '  theme(legend.position = "none",',
    '        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"))',
    'if (nrow(sig) > 0) p <- p + stat_pvalue_manual(sig, label = "lab", tip.length = 0.01, bracket.size = 0.4)',
    "",
    sprintf('ggsave("%s.pdf", p, width = 4.2, height = 4.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
