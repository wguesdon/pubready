# Recipe: compare a continuous outcome across 3+ groups (also works for 2).
# Picks one-way ANOVA, Welch's ANOVA, or Kruskal-Wallis from the assumption
# checks, runs the matching post-hoc test, and brackets only the comparisons
# that reach significance (adjusted P < 0.05). The full pairwise table is saved.

# When the normality check fails but the Kruskal-Wallis test cannot reach
# P < 0.05 at the group sizes (rank_min_p), --test auto keeps the parametric
# test that Levene's test picks, when that test can run (can_run, a named
# logical with the entries anova and welch).
resolve_multi_group_method <- function(method, all_normal, equal_var, rank_min_p = NA_real_,
                                       can_run = c(anova = TRUE, welch = TRUE)) {
  method <- tolower(method %||% "auto")
  para <- if (equal_var) "anova" else "welch"
  kept <- method == "auto" && !all_normal && rank_test_cannot_reach(rank_min_p) &&
    isTRUE(can_run[[para]])
  fam <- if (method == "auto") {
    if (all_normal || kept) para else "kruskal"
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
                                 kruskal = "Dunn's test with Benjamini-Hochberg adjustment"),
    kept_parametric = kept,
    rank_min_p    = rank_min_p
  )
}

# Significance stars with the cutpoints of rstatix::add_significance.
.mg_signif <- function(p) {
  as.character(cut(p, breaks = c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                   labels = c("****", "***", "**", "*", "ns"), right = FALSE))
}

# Tukey's HSD from stats::TukeyHSD, in the pair order of combn(levels, 2). The
# estimate is the mean of group2 minus the mean of group1, as rstatix reports it.
.mg_tukey <- function(df, fml, lvls) {
  tk <- stats::TukeyHSD(stats::aov(fml, data = df))[[1]]
  pairs <- utils::combn(lvls, 2)
  key <- paste(pairs[2, ], pairs[1, ], sep = "-")
  p <- unname(tk[key, "p adj"])
  data.frame(group1 = pairs[1, ], group2 = pairs[2, ], estimate = unname(tk[key, "diff"]),
             p.adj = p, p.adj.signif = .mg_signif(p), stringsAsFactors = FALSE)
}

# One Games-Howell comparison: a Welch t on the pair, referred to the studentized
# range distribution with the Welch degrees of freedom, as in rstatix and pingouin.
.mg_gh_pair <- function(k, pairs, by_group, n_groups) {
  a <- by_group[[pairs[1, k]]]; b <- by_group[[pairs[2, k]]]
  va <- stats::var(a) / length(a); vb <- stats::var(b) / length(b)
  diff <- mean(b) - mean(a)
  tval <- abs(diff) / sqrt(va + vb)
  dfw <- (va + vb)^2 / (va^2 / (length(a) - 1) + vb^2 / (length(b) - 1))
  c(diff, stats::ptukey(tval * sqrt(2), n_groups, dfw, lower.tail = FALSE))
}

.mg_games_howell <- function(by_group, lvls) {
  pairs <- utils::combn(lvls, 2)
  res <- vapply(seq_len(ncol(pairs)), .mg_gh_pair, numeric(2), pairs = pairs,
                by_group = by_group, n_groups = length(lvls))
  data.frame(group1 = pairs[1, ], group2 = pairs[2, ], estimate = res[1, ],
             p.adj = res[2, ], p.adj.signif = .mg_signif(res[2, ]), stringsAsFactors = FALSE)
}

# The --test values that multi_group_compare knows (F63).
.multi_group_compare_tests <- c("auto", "anova", "aov", "welch", "welch_anova", "kruskal",
                                "kruskal_wallis", "kw")

recipe_multi_group_compare <- function(df, spec) {
  check_test_method("multi_group_compare", spec$test$method, .multi_group_compare_tests)
  cl <- clean_tidy(df, spec)
  df <- cl$df
  x <- spec$data$x
  y <- spec$data$y

  # One level rule in both engines: the --levels order, else numeric order when
  # every label is a number, else code-point order.
  lab <- level_labels(df[[x]])
  lvls <- level_order(lab, unlist(spec$data$levels), "group")
  df[[x]] <- factor(lab, levels = lvls)
  if (length(lvls) < 2) {
    stop(sprintf("multi_group_compare needs 2 or more groups in '%s'; found %d", x, length(lvls)))
  }

  # The tests run on a copy with syntactic column names (.x, .y), so a column
  # name with a space or a symbol cannot break a formula. The labels and the
  # tables keep the original names.
  sdf <- data.frame(.x = df[[x]], .y = df[[y]])
  fml <- .y ~ .x
  by_group <- split(df[[y]], df[[x]])
  norm <- lapply(by_group, shapiro_safe)
  all_normal <- all(vapply(norm, is_normal_check, logical(1)))
  lev_p <- tryCatch(rstatix::levene_test(sdf, fml)$p[1], error = na_on_error)
  equal_var <- is.na(lev_p) || lev_p > 0.05

  rank_min_p <- rank_test_min_p(lengths(by_group), "kruskal")
  zv <- zero_variance_groups(by_group)
  resolved <- resolve_multi_group_method(spec$test$method, all_normal, equal_var, rank_min_p,
                                         c(anova = length(zv) < length(by_group),
                                           welch = length(zv) == 0))
  .mg_check_variance(resolved, by_group)

  # Omnibus test and matching post-hoc. The base R functions give full precision;
  # rstatix rounds the Tukey and Games-Howell P values to three decimals.
  if (resolved$family == "anova") {
    at <- stats::anova(stats::lm(fml, data = sdf))
    omnibus <- list(statistic = at[["F value"]][1], df1 = at$Df[1], df2 = at$Df[2],
                    p = at[["Pr(>F)"]][1],
                    effect_size = list(name = "partial eta-squared",
                                       value = at[["Sum Sq"]][1] / sum(at[["Sum Sq"]])))
    ph <- .mg_tukey(sdf, fml, lvls)
  } else if (resolved$family == "welch") {
    wt <- stats::oneway.test(fml, data = sdf, var.equal = FALSE)
    omnibus <- list(statistic = unname(wt$statistic), df1 = unname(wt$parameter[1]),
                    df2 = unname(wt$parameter[2]), p = wt$p.value, effect_size = NULL)
    ph <- .mg_games_howell(by_group, lvls)
  } else {
    kt <- stats::kruskal.test(fml, data = sdf)
    ke <- tryCatch(as.data.frame(rstatix::kruskal_effsize(sdf, fml)), error = null_on_error)
    omnibus <- list(statistic = unname(kt$statistic), df1 = unname(kt$parameter), df2 = NA_real_,
                    p = kt$p.value,
                    effect_size = if (!is.null(ke)) list(name = "eta-squared (H)", value = ke$effsize[1]) else NULL)
    ph <- as.data.frame(rstatix::dunn_test(sdf, fml, p.adjust.method = "BH"))
  }

  ph$group1 <- as.character(ph$group1)
  ph$group2 <- as.character(ph$group2)
  estimate <- if ("estimate" %in% names(ph)) ph$estimate else NA_real_

  stats_df <- data.frame(
    outcome      = y,
    group1       = ph$group1,
    group2       = ph$group2,
    contrast     = paste(ph$group2, "-", ph$group1),
    n1           = vapply(ph$group1, .mg_group_n, integer(1), by_group = by_group),
    n2           = vapply(ph$group2, .mg_group_n, integer(1), by_group = by_group),
    post_hoc     = resolved$posthoc_label,
    estimate     = estimate,
    p_adj        = ph$p.adj,
    p_adj_signif = ph$p.adj.signif,
    displayed    = ph$p.adj < 0.05,
    omnibus_test = resolved$label,
    omnibus_statistic = omnibus$statistic,
    omnibus_df1  = omnibus$df1,
    omnibus_df2  = omnibus$df2,
    omnibus_p    = omnibus$p,
    omnibus_effect_size_name = if (is.null(omnibus$effect_size)) NA_character_ else omnibus$effect_size$name,
    omnibus_effect_size = if (is.null(omnibus$effect_size)) NA_real_ else omnibus$effect_size$value,
    stringsAsFactors = FALSE
  )

  # Significant comparisons only. Up to .mg_max_brackets() of them get a
  # bracket each, on its own level; more than that would stack too many lines,
  # so the figure shows a compact letter display instead.
  sig <- ph[which(ph$p.adj < 0.05), , drop = FALSE]
  display <- if (nrow(sig) > .mg_max_brackets()) "compact letter display" else "brackets"
  bracket_label <- spec$appearance$bracket$label %||% "p.signif"

  pal  <- spec$appearance$palette %||% pubready_palette
  pal  <- rep(pal, length.out = length(lvls))
  geom <- spec$appearance$geom %||% "box"

  # The brackets and the letters start above the drawn geom. The spacing is a
  # share of the span of the drawn axis (from 0 for a bar), so a label never
  # reaches the next line.
  tops <- vapply(lvls, .mg_group_top, numeric(1), by_group = by_group, geom = geom)
  ylo <- min(df[[y]], na.rm = TRUE)
  if (identical(geom, "bar")) ylo <- min(ylo, 0)
  top_all <- max(tops)
  span <- top_all - ylo
  if (!is.finite(span) || span <= 0) span <- max(abs(top_all), 1)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  p <- add_dist_geom(p, x, geom, spec$appearance$show_points)
  omni_lab <- sprintf("%s, P = %s", cap_first(resolved$label),
                      formatC(omnibus$p, format = "g", digits = 2))
  p <- p +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = stats::setNames(wrap_labels(lvls), lvls)) +
    ggplot2::scale_y_continuous(labels = plain_number_labels) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% y,
                  title = spec$appearance$title,
                  subtitle = omni_lab) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))
  cld <- NULL
  if (display == "brackets" && nrow(sig) > 0) {
    sig <- .mg_bracket_rows(sig, lvls, top_all, span, bracket_label)
    p <- p +
      ggpubr::stat_pvalue_manual(sig, label = "lab", tip.length = 0.01, bracket.size = 0.4,
                                 label.size = 3.6, vjust = 0) +
      ggplot2::expand_limits(y = max(sig$y.position) + 0.10 * span)
  } else if (display == "compact letter display") {
    cld <- cld_letters(lvls, ph$group1, ph$group2, ph$p.adj)
    lab_df <- data.frame(xg = factor(lvls, levels = lvls), yv = tops + 0.04 * span,
                         lab = unname(cld), stringsAsFactors = FALSE)
    p <- p +
      ggplot2::geom_text(data = lab_df, ggplot2::aes(x = .data$xg, y = .data$yv, label = .data$lab),
                         inherit.aes = FALSE, vjust = 0, size = 3.6) +
      ggplot2::expand_limits(y = top_all + 0.14 * span)
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
    display = display,
    letters = if (!is.null(cld)) as.list(cld) else NULL,
    assumptions = list(
      normality      = list(test = "shapiro-wilk", all_normal = all_normal),
      equal_variance = list(test = "Levene", p = lev_p)
    ),
    rank_test_attainable_p = rank_min_p
  )

  n_by_group <- vapply(by_group, length, integer(1))
  methods <- .mg_methods(resolved, omnibus, spec, n_by_group, display, nrow(sig))
  build_script <- list(fn = .mg_emit_script,
                       args = list(spec = spec, resolved = resolved, lvls = lvls,
                                   pal = pal[seq_len(length(lvls))], cld = cld,
                                   width = category_width(length(lvls), wrap_labels(lvls))))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       df_used = df, clean_steps = cl$steps,
       qc = list(quantity = "group values", panels = by_group),
       width = category_width(length(lvls), wrap_labels(lvls)), height = 4.2)
}

# Welch's ANOVA weights each group by the inverse of its variance, and the
# one-way ANOVA divides by the pooled variance. A variance of zero there makes
# the statistic infinite, so the run stops and names the groups. The Python
# engine has the same rule.
.mg_check_variance <- function(resolved, by_group) {
  zv <- zero_variance_groups(by_group)
  if (resolved$family == "welch" && length(zv) > 0) {
    stop(zero_variance_message(resolved$label, zero_variance_subject(zv),
                               "use --test anova or --test kruskal"))
  }
  if (resolved$family == "anova" && length(zv) == length(by_group)) {
    stop(zero_variance_message(resolved$label, zero_variance_subject(zv), "use --test kruskal"))
  }
  invisible(NULL)
}

# More significant pairs than this switch the figure to a compact letter display.
.mg_max_brackets <- function() 6

.mg_group_top <- function(g, by_group, geom) geom_top(by_group[[g]], geom)

.mg_group_n <- function(g, by_group) length(by_group[[g]])

# One level per bracket, the narrowest pairs lowest, from 8% of the span above
# the tallest geom, 13% of the span apart. The Python engine uses the same rule.
.mg_bracket_rows <- function(sig, lvls, top_all, span, bracket_label) {
  i1 <- match(sig$group1, lvls); i2 <- match(sig$group2, lvls)
  sig <- sig[order(abs(i2 - i1), pmin(i1, i2)), , drop = FALSE]
  sig$y.position <- top_all + span * (0.08 + 0.13 * (seq_len(nrow(sig)) - 1))
  sig$lab <- if (identical(bracket_label, "p.format")) {
    paste0("p = ", formatC(sig$p.adj, format = "g", digits = 2))
  } else {
    p_stars(sig$p.adj)
  }
  sig
}

.mg_p_text <- function(p) formatC(p, format = "g", digits = 2)

.mg_methods <- function(resolved, omnibus, spec, n_by_group, display, n_sig) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v    <- pkg_versions(c("rstatix", "ggpubr", "ggprism"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  omni <- if (resolved$family == "kruskal") {
    sprintf("a Kruskal-Wallis test (H(%s) = %.2f, P = %s)",
            omnibus$df1, omnibus$statistic, .mg_p_text(omnibus$p))
  } else {
    sprintf("%s%s (F(%s, %s) = %.2f, P = %s)",
            resolved$article, resolved$label,
            df_text(omnibus$df1), df_text(omnibus$df2), omnibus$statistic, .mg_p_text(omnibus$p))
  }
  eff <- if (!is.null(omnibus$effect_size)) {
    sprintf(" The effect size was %s = %.2f.", omnibus$effect_size$name, omnibus$effect_size$value)
  } else ""
  if (isTRUE(resolved$kept_parametric) || resolved$family == "kruskal") {
    floor_txt <- rank_floor_text(resolved$rank_min_p, resolved$kept_parametric)
    if (nzchar(floor_txt)) eff <- paste0(" ", floor_txt, eff)
  }
  ph_txt <- if (identical(display, "compact letter display")) {
    sprintf(" Pairwise comparisons used %s. Because %d pairs differed (adjusted P < 0.05), the figure shows a compact letter display (the insert-and-absorb algorithm of Piepho, 2004) instead of brackets: groups that share a letter do not differ significantly.",
            resolved$posthoc_label, n_sig)
  } else {
    sprintf(" Pairwise comparisons used %s; only comparisons reaching significance (adjusted P < 0.05) are annotated on the figure with brackets.",
            resolved$posthoc_label)
  }
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

.mg_emit_script <- function(spec, resolved, lvls, pal, cld, width, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  gl <- script_geom_lines(x, spec$appearance$geom %||% "box",
                          spec$appearance$show_points %||% TRUE)
  posthoc_line <- switch(resolved$posthoc,
    tukey        = 'ph <- rstatix::tukey_hsd(sdf, .y ~ .x)',
    games_howell = 'ph <- rstatix::games_howell_test(sdf, .y ~ .x)',
    dunn         = 'ph <- rstatix::dunn_test(sdf, .y ~ .x, p.adjust.method = "BH")'
  )
  omni_line <- switch(resolved$family,
    anova   = 'omni <- as.data.frame(rstatix::anova_test(sdf, .y ~ .x))',
    welch   = 'omni <- as.data.frame(rstatix::welch_anova_test(sdf, .y ~ .x))',
    kruskal = 'omni <- as.data.frame(rstatix::kruskal_test(sdf, .y ~ .x))'
  )
  omni_name <- cap_first(resolved$label)
  lab_line <- if (identical(spec$appearance$bracket$label %||% "p.signif", "p.format")) {
    'sig$lab <- paste0("p = ", formatC(sig$p.adj, format = "g", digits = 2))'
  } else {
    "sig$lab <- sig$p.adj.signif"
  }
  # With more than six significant pairs the figure shows letters, not brackets.
  annot_lines <- if (!is.null(cld)) {
    c("# Compact letter display (Piepho 2004), as computed by pubready.",
      sprintf("cld <- c(%s)", paste(sprintf('"%s" = "%s"', names(cld), unname(cld)), collapse = ", ")),
      sprintf('tops <- tapply(df[["%s"]], df[["%s"]], max)', y, x),
      sprintf('lab <- data.frame(g = factor(names(cld), levels = names(cld)), lab = unname(cld), y = tops[names(cld)] + 0.04 * diff(range(df[["%s"]])))', y),
      'p <- p + geom_text(data = lab, aes(x = g, y = y, label = lab), inherit.aes = FALSE, vjust = 0)')
  } else {
    c("if (nrow(sig) > 0) {",
      sprintf('  ymax <- max(df[["%s"]], na.rm = TRUE); rng <- diff(range(df[["%s"]], na.rm = TRUE))', y, y),
      "  if (!is.finite(rng) || rng == 0) rng <- max(abs(ymax), 1)",
      sprintf('  lv <- levels(df[["%s"]]); i1 <- match(sig$group1, lv); i2 <- match(sig$group2, lv)', x),
      "  sig <- sig[order(abs(i2 - i1), pmin(i1, i2)), , drop = FALSE]",
      "  sig$y.position <- ymax + rng * (0.08 + 0.13 * (seq_len(nrow(sig)) - 1))",
      sprintf("  %s", lab_line),
      '  p <- p + stat_pvalue_manual(sig, label = "lab", tip.length = 0.01, bracket.size = 0.4)',
      "}")
  }

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(ggpubr); library(rstatix); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf('df[["%s"]] <- factor(df[["%s"]], levels = c(%s))', x, x,
            paste(sprintf('"%s"', lvls), collapse = ", ")),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "# The tests run on a copy with syntactic column names.",
    sprintf('sdf <- data.frame(.x = df[["%s"]], .y = df[["%s"]])', x, y),
    "",
    posthoc_line,
    "ph$group1 <- as.character(ph$group1); ph$group2 <- as.character(ph$group2)",
    "sig <- ph[which(ph$p.adj < 0.05), , drop = FALSE]",
    "",
    omni_line,
    sprintf('subtitle_txt <- paste0("%s, P = ", formatC(omni$p[1], format = "g", digits = 2))', omni_name),
    "",
    gl$prep,
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]])) +', x, y),
    gl$layers,
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  labs(x = "%s", y = "%s", subtitle = subtitle_txt) +',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    "  theme_classic(base_size = 13) +",
    '  theme(legend.position = "none",',
    '        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"))',
    annot_lines,
    "",
    sprintf('ggsave("%s.pdf", p, width = %.1f, height = 4.2)', fig_stub, width),
    'cat("Reproduced figure.\\n")'
  )
}
