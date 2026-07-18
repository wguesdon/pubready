# Recipe: two-way and three-way factorial ANOVA on one continuous outcome.
# Parametric (Type II ANOVA) by default; falls back to the aligned rank
# transform (ARTool) when the residuals are non-normal, or when forced with
# --test art. The grouped figure shows every effect's p-value on top; the full
# effects table is saved.
#
# Factors: --x is the x-axis factor, --fill the grouping factor (second factor),
# --facet the optional third factor.

.fa_plab <- function(p) {
  if (is.na(p)) "P = NA"
  else if (p < 0.001) "P < 0.001"
  else paste0("P = ", formatC(p, format = "g", digits = 2))
}

recipe_factorial_anova <- function(df, spec) {
  y  <- spec$data$y
  f1 <- spec$data$x
  f2 <- spec$data$fill
  f3 <- spec$data$facet
  if (is.null(f2) || !nzchar(f2)) {
    stop("factorial_anova needs at least two factors: pass --x and --fill (and --facet for three-way)")
  }
  factors <- c(f1, f2, if (!is.null(f3) && nzchar(f3)) f3)
  for (f in factors) if (!f %in% names(df)) stop(sprintf("factor column '%s' not found", f))
  if (!y %in% names(df)) stop(sprintf("outcome column '%s' not found", y))
  for (f in factors) df[[f]] <- factor(df[[f]])
  df[[y]] <- as.numeric(df[[y]])
  n0 <- nrow(df)
  df <- df[stats::complete.cases(df[, c(y, factors)]), , drop = FALSE]
  clean_steps <- if (nrow(df) < n0) {
    sprintf("Dropped %d row(s) with missing values in %s (%d -> %d rows).",
            n0 - nrow(df), paste(c(y, factors), collapse = ", "), n0, nrow(df))
  } else "No cleaning applied; input used as-is."

  rhs <- paste(sprintf("`%s`", factors), collapse = " * ")
  fml <- stats::as.formula(sprintf("`%s` ~ %s", y, rhs))

  method <- tolower(spec$test$method %||% "auto")
  res_p  <- tryCatch(shapiro_safe(stats::residuals(stats::aov(fml, data = df)))$p,
                     error = function(e) NA_real_)
  use_art <- if (method %in% c("anova", "aov")) FALSE
             else if (method %in% c("art", "aligned_rank")) TRUE
             else (!is.na(res_p) && res_p < 0.05)   # auto

  if (use_art) {
    m  <- ARTool::art(fml, data = df)
    at <- as.data.frame(stats::anova(m))
    eff <- data.frame(
      effect      = at$Term,
      df1         = at$Df,
      df2         = at$Df.res,
      statistic   = at$F,
      p_value     = at[["Pr(>F)"]],
      effect_size = NA_real_,
      stringsAsFactors = FALSE
    )
    test_name <- sprintf("aligned rank transform ANOVA (%d-way)", length(factors))
  } else {
    at <- as.data.frame(rstatix::anova_test(df, fml, type = 2))
    eff <- data.frame(
      effect      = at$Effect,
      df1         = at$DFn,
      df2         = at$DFd,
      statistic   = at$F,
      p_value     = at$p,
      effect_size = at$ges,
      stringsAsFactors = FALSE
    )
    test_name <- sprintf("%d-way ANOVA", length(factors))
  }
  eff$significant <- eff$p_value < 0.05

  # Subtitle: significant effects only, interaction ":" shown as a cross.
  sig_eff <- eff[which(eff$significant), , drop = FALSE]
  subtitle <- if (nrow(sig_eff) > 0) {
    paste(sprintf("%s: %s", gsub(":", " × ", sig_eff$effect),
                  vapply(sig_eff$p_value, .fa_plab, character(1))),
          collapse = "\n")
  } else {
    "No significant effects (P ≥ 0.05)"
  }

  pal  <- spec$appearance$palette %||% pubplot_palette
  pal  <- rep(pal, length.out = nlevels(df[[f2]]))
  geom <- spec$appearance$geom %||% "box"
  pd   <- ggplot2::position_dodge(width = 0.8)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[f1]], y = .data[[y]], fill = .data[[f2]]))
  if (geom == "bar") {
    p <- p +
      ggplot2::stat_summary(fun = mean, geom = "bar", position = pd, width = 0.7,
                            alpha = 0.9, colour = "black", linewidth = 0.3) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                            position = pd, width = 0.2)
  } else {
    p <- p + ggplot2::geom_boxplot(position = pd, width = 0.7, outlier.shape = NA, alpha = 0.9)
  }
  if (isTRUE(spec$appearance$show_points)) {
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = 1.3, alpha = 0.7, shape = 21, stroke = 0.3)
  }
  p <- p +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = spec$appearance$x_label %||% f1,
                  y = spec$appearance$y_label %||% y,
                  fill = f2, title = spec$appearance$title, subtitle = subtitle) +
    theme_pubplot() +
    ggplot2::theme(legend.position = "right",
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
  if (!is.null(f3) && nzchar(f3)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[f3]]))
  }

  stats_df <- eff[, c("effect", "df1", "df2", "statistic", "p_value", "effect_size", "significant")]

  test_meta <- list(
    name    = test_name,
    factors = as.list(factors),
    type    = if (use_art) "aligned_rank_transform" else "type_II_anova",
    residual_normality = list(test = "shapiro-wilk", p = res_p),
    effects = lapply(seq_len(nrow(eff)), function(i) list(
      effect = eff$effect[i], df1 = eff$df1[i], df2 = eff$df2[i],
      statistic = eff$statistic[i], p_value = eff$p_value[i],
      effect_size = if (is.na(eff$effect_size[i])) NULL else eff$effect_size[i]
    ))
  )

  resolved <- list(family = if (use_art) "art" else "anova", label = test_name,
                   factors = factors, use_art = use_art)

  methods <- .fa_methods(eff, spec, factors, use_art)
  build_script <- function(in_name, fig_stub) {
    .fa_emit_script(spec, factors, use_art, pal[seq_len(nlevels(df[[f2]]))], in_name, fig_stub)
  }

  has_facet <- !is.null(f3) && nzchar(f3)
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       df_used = df, clean_steps = clean_steps,
       width = if (has_facet) 6.8 else 5.2, height = 4.4)
}

.fa_methods <- function(eff, spec, factors, use_art) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v    <- pkg_versions(c("rstatix", "ARTool", "ggplot2"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  intro <- if (use_art) {
    sprintf("%s was analyzed by an aligned rank transform ANOVA, a non-parametric factorial method, with factors %s.",
            cap_first(y_label), paste(factors, collapse = ", "))
  } else {
    sprintf("%s was analyzed by a %d-way ANOVA (Type II sums of squares) with factors %s.",
            cap_first(y_label), length(factors), paste(factors, collapse = ", "))
  }
  parts <- vapply(seq_len(nrow(eff)), function(i) {
    e <- eff[i, ]
    sprintf("%s (F(%s, %s) = %.2f, %s)",
            gsub(":", " x ", e$effect), e$df1, e$df2, e$statistic, .fa_plab(e$p_value))
  }, character(1))
  soft <- if (use_art) sprintf("ARTool %s", v$ARTool) else sprintf("rstatix %s", v$rstatix)

  paste0(
    intro, " Effects tested: ", paste(parts, collapse = "; "), ".",
    " Significance was set at P < 0.05.",
    sprintf(" Analyses were performed in R %s with %s; figures were produced with ggplot2 %s.",
            rver, soft, v$ggplot2)
  )
}

.fa_emit_script <- function(spec, factors, use_art, pal, in_name, fig_stub) {
  y  <- spec$data$y
  f1 <- factors[1]; f2 <- factors[2]; f3 <- if (length(factors) >= 3) factors[3] else NULL
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  rhs <- paste(sprintf("`%s`", factors), collapse = " * ")
  geom <- spec$appearance$geom %||% "box"

  omni_lines <- if (use_art) c(
    sprintf('m  <- ARTool::art(%s ~ %s, data = df)', y, rhs),
    'at <- as.data.frame(anova(m))',
    'eff <- data.frame(effect = at$Term, p = at[["Pr(>F)"]])'
  ) else c(
    sprintf('at <- as.data.frame(rstatix::anova_test(df, %s ~ %s, type = 2))', y, rhs),
    'eff <- data.frame(effect = at$Effect, p = at$p)'
  )
  geom_lines <- if (geom == "bar") c(
    '  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.8), width = 0.7, alpha = 0.9, colour = "black", linewidth = 0.3) +',
    '  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.8), width = 0.2) +'
  ) else
    '  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.9) +'
  facet_line <- if (!is.null(f3)) sprintf('p <- p + facet_wrap(vars(.data[["%s"]]))', f3) else character(0)

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubplot container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(rstatix); library(readr); library(ARTool)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    paste0(sprintf('df[["%s"]] <- factor(df[["%s"]]); ', factors, factors), collapse = ""),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "",
    omni_lines,
    'plab <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))',
    'sig <- eff[which(eff$p < 0.05), , drop = FALSE]',
    'subtitle_txt <- if (nrow(sig) > 0) paste(sprintf("%s: %s", gsub(":", " × ", sig$effect), vapply(sig$p, plab, character(1))), collapse = "\\n") else "No significant effects (P ≥ 0.05)"',
    "",
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]], fill = .data[["%s"]])) +', f1, y, f2),
    geom_lines,
    '  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.3, alpha = 0.7, shape = 21, stroke = 0.3) +',
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  labs(x = "%s", y = "%s", fill = "%s", subtitle = subtitle_txt) +',
            spec$appearance$x_label %||% f1, spec$appearance$y_label %||% y, f2),
    "  theme_classic(base_size = 13) +",
    '  theme(plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey30"))',
    facet_line,
    "",
    sprintf('ggsave("%s.pdf", p, width = %s, height = 4.4)', fig_stub, if (!is.null(f3)) "6.8" else "5.2"),
    'cat("Reproduced figure.\\n")'
  )
}
