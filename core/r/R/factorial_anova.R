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

# The --test values that factorial_anova knows (F63).
.factorial_anova_tests <- c("auto", "anova", "aov", "art", "aligned_rank")

recipe_factorial_anova <- function(df, spec) {
  check_test_method("factorial_anova", spec$test$method, .factorial_anova_tests)
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
  # Each reason for a dropped row is counted on its own: a text value of the
  # outcome that becomes NA, and a missing value in each column.
  y_raw <- df[[y]]
  df[[y]] <- suppressWarnings(as.numeric(df[[y]]))
  coerced <- sum(is.na(df[[y]]) & !is.na(y_raw))
  n_missing <- vapply(c(y, factors), .fa_n_missing, integer(1), df = df)
  n0 <- nrow(df)
  df <- df[stats::complete.cases(df[, c(y, factors)]), , drop = FALSE]
  # The factors are made after the rows are dropped, so a level whose every
  # outcome is missing does not stay as an empty slot of the axis. --levels
  # orders the --x levels; --fill and --facet take the rule without a list:
  # numeric order when every label is a number, else code-point order.
  df[[f1]] <- .fa_factor(df[[f1]], unlist(spec$data$levels))
  for (f in factors[-1]) df[[f]] <- .fa_factor(df[[f]])
  clean_steps <- .fa_clean_steps(y, coerced, n_missing, n0, nrow(df))
  check_spread(df[[y]], "factorial_anova", sprintf("the outcome '%s'", y))
  .fa_check_design(df, factors)

  # The models run on a copy with syntactic column names (.y, .f1, .f2, .f3),
  # so a column name with a space or a symbol cannot break a formula, also in
  # ARTool. The effects take the original names back.
  safe <- sprintf(".f%d", seq_along(factors))
  sdf <- stats::setNames(data.frame(df[[y]], df[factors]), c(".y", safe))
  fml <- stats::as.formula(sprintf(".y ~ %s", paste(safe, collapse = " * ")))

  method <- tolower(spec$test$method %||% "auto")
  resid_vec <- tryCatch(as.numeric(stats::residuals(stats::aov(fml, data = sdf))),
                        error = .fa_null_on_error)
  res_p  <- if (!is.null(resid_vec)) shapiro_safe(resid_vec)$p else NA_real_
  use_art <- if (method %in% c("anova", "aov")) FALSE
             else if (method %in% c("art", "aligned_rank")) TRUE
             else (!is.na(res_p) && res_p < 0.05)   # auto

  if (use_art) {
    m  <- ARTool::art(fml, data = sdf)
    at <- as.data.frame(stats::anova(m))
    eff <- data.frame(
      effect      = vapply(as.character(at$Term), .fa_term_name, character(1),
                           safe = safe, factors = factors, USE.NAMES = FALSE),
      df1         = at$Df,
      df2         = at$Df.res,
      statistic   = at$F,
      p_value     = at[["Pr(>F)"]],
      effect_size = NA_real_,
      stringsAsFactors = FALSE
    )
    test_name <- sprintf("aligned rank transform ANOVA (%d-way)", length(factors))
  } else {
    # car::Anova gives the Type II table at full precision; rstatix::anova_test
    # rounds F to three decimals and P to three significant digits. The effect
    # size is the partial eta squared, which the Python engine also reports.
    at <- as.data.frame(car::Anova(stats::lm(fml, data = sdf), type = 2))
    terms <- setdiff(rownames(at), "Residuals")
    ss_res <- at["Residuals", "Sum Sq"]
    eff <- data.frame(
      effect      = vapply(terms, .fa_term_name, character(1), safe = safe,
                           factors = factors, USE.NAMES = FALSE),
      df1         = at[terms, "Df"],
      df2         = at["Residuals", "Df"],
      statistic   = at[terms, "F value"],
      p_value     = at[terms, "Pr(>F)"],
      effect_size = at[terms, "Sum Sq"] / (at[terms, "Sum Sq"] + ss_res),
      stringsAsFactors = FALSE
    )
    test_name <- sprintf("%d-way ANOVA", length(factors))
  }
  eff$effect_size_name <- if (use_art) NA_character_ else "partial eta squared"
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

  pal  <- spec$appearance$palette %||% pubready_palette
  pal  <- rep(pal, length.out = nlevels(df[[f2]]))
  geom <- spec$appearance$geom %||% "box"
  if (!geom %in% c("box", "bar", "interaction")) geom <- "box"
  pd   <- ggplot2::position_dodge(width = 0.8)

  if (geom == "interaction") {
    # The interaction plot: the mean +/- SEM of each cell, one line per --fill
    # level across the --x levels, slightly dodged so the error bars do not
    # overlap. The statistics are those of the box geom.
    summ <- .fa_cell_summary(df, y, factors)
    pdi <- ggplot2::position_dodge(width = .fa_int_dodge)
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data[[f1]], y = mean, colour = .data[[f2]],
                                            group = .data[[f2]])) +
      ggplot2::geom_line(position = pdi, linewidth = 0.8) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sem, ymax = mean + sem), width = 0.12,
                             position = pdi, linewidth = 0.6, na.rm = TRUE) +
      ggplot2::geom_point(position = pdi, size = 2.4) +
      ggplot2::scale_colour_manual(values = pal) +
      ggplot2::scale_y_continuous(breaks = scales::breaks_extended(n = 6))
  } else {
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
      p <- p + ggplot2::geom_point(data = .fa_point_x(df, f1, f2, f3),
                                   ggplot2::aes(x = .data$.point_x),
                                   size = 1.3, alpha = 0.7, shape = 21, stroke = 0.3)
    }
    p <- p + ggplot2::scale_fill_manual(values = pal)
  }
  # Long --x levels wrap at 12 characters, so neighbouring labels do not collide.
  if (.fa_needs_wrap(levels(df[[f1]]))) {
    p <- p + ggplot2::scale_x_discrete(labels = .fa_wrap_labels)
  }
  p <- p +
    ggplot2::labs(x = spec$appearance$x_label %||% f1,
                  y = spec$appearance$y_label %||% y,
                  fill = f2, colour = f2, title = spec$appearance$title, subtitle = subtitle) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
  if (!is.null(f3) && nzchar(f3)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[f3]]))
  }

  stats_df <- eff[, c("effect", "df1", "df2", "statistic", "p_value", "effect_size",
                      "effect_size_name", "significant")]

  test_meta <- list(
    name    = test_name,
    factors = as.list(factors),
    type    = if (use_art) "aligned_rank_transform" else "type_II_anova",
    residual_normality = list(test = "shapiro-wilk", p = res_p),
    effects = lapply(seq_len(nrow(eff)), .fa_meta_effect, eff = eff)
  )

  resolved <- list(family = if (use_art) "art" else "anova", label = test_name,
                   factors = factors, use_art = use_art)

  methods <- .fa_methods(eff, spec, factors, use_art)
  has_facet <- !is.null(f3) && nzchar(f3)
  size <- .fa_fig_size(levels(df[[f1]]), levels(df[[f2]]), f2,
                       if (has_facet) nlevels(df[[f3]]) else 1, has_facet)
  build_script <- list(fn = .fa_emit_script,
                       args = list(spec = spec, factors = factors, use_art = use_art,
                                   pal = pal[seq_len(nlevels(df[[f2]]))],
                                   width = size$width, height = size$height,
                                   wrap_x = .fa_needs_wrap(levels(df[[f1]])),
                                   fac_levels = lapply(df[factors], levels)))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       df_used = df, clean_steps = clean_steps,
       qc = if (!is.null(resid_vec))
              list(quantity = "model residuals", panels = list(residuals = resid_vec))
            else NULL,
       width = size$width, height = size$height)
}

.fa_null_on_error <- function(e) NULL

# The original name of a model term: ".f1:.f2" becomes "genotype:treatment".
.fa_term_name <- function(term, safe, factors) {
  paste(factors[match(strsplit(term, ":", fixed = TRUE)[[1]], safe)], collapse = ":")
}

# The design check before the fit, in the same words as the Python engine. Each
# factor needs 2 or more levels (F92). Every cell of the design needs a value,
# because an empty cell makes some terms not estimable and the two engines then
# test different hypotheses (F53). Some cell needs 2 or more values, else no
# residual degree of freedom is left (F92).
.fa_check_design <- function(df, factors) {
  for (f in factors) {
    lv <- levels(df[[f]])
    if (length(lv) < 2) {
      stop(sprintf("factorial_anova needs 2 or more levels in each factor; '%s' has %d (%s)",
                   f, length(lv), quote_join(lv)))
    }
  }
  cells <- as.data.frame(table(df[factors]), stringsAsFactors = FALSE)
  names(cells) <- c(factors, "Freq")
  empty <- cells[cells$Freq == 0, factors, drop = FALSE]
  if (nrow(empty) > 0) {
    txt <- vapply(seq_len(nrow(empty)), .fa_cell_text, character(1), cells = empty,
                  factors = factors)
    more <- if (length(txt) > 5) sprintf("; and %d more", length(txt) - 5) else ""
    stop(sprintf("factorial_anova needs a value in every cell of the design; no row holds %s%s",
                 paste(utils::head(txt, 5), collapse = "; "), more))
  }
  if (nrow(df) <= nrow(cells)) {
    stop(sprintf("factorial_anova needs a cell with 2 or more values to estimate the residual variance; the %d cells hold one value each",
                 nrow(cells)))
  }
  invisible(NULL)
}

# One cell of the design as text: genotype = 'ko', treatment = 'drug'.
.fa_cell_text <- function(i, cells, factors) {
  paste(sprintf("%s = '%s'", factors, unlist(cells[i, factors], use.names = FALSE)),
        collapse = ", ")
}

# The count of missing values of one column.
.fa_n_missing <- function(col, df) sum(is.na(df[[col]]))

# The data log of the cleaning: the coerced outcome values, the rows with no
# value in each column (a row with two missing values counts in both), and the
# rows dropped in all. The Python engine writes the same steps.
.fa_clean_steps <- function(y, coerced, n_missing, n0, n1) {
  steps <- character(0)
  if (coerced > 0) {
    steps <- c(steps, sprintf("Coerced column '%s' to numeric; %d value(s) became NA.", y, coerced))
  }
  for (col in names(n_missing)[n_missing > 0]) {
    steps <- c(steps, sprintf("%d row(s) have no value in '%s'.", n_missing[[col]], col))
  }
  if (n1 < n0) {
    steps <- c(steps, sprintf("Dropped %d row(s) with a missing value (%d -> %d rows).",
                              n0 - n1, n0, n1))
  }
  if (length(steps) == 0) "No cleaning applied; input used as-is." else steps
}

# A factor column with the level rule of both engines (util level_order).
.fa_factor <- function(v, levels = NULL) {
  lab <- level_labels(v)
  factor(lab, levels = level_order(lab, levels, "level"))
}

# The line of an emitted script that makes one factor with the levels of the
# figure.
.fa_factor_line <- function(i, factors, fac_levels) {
  sprintf('df[["%s"]] <- factor(df[["%s"]], levels = c(%s))', factors[i], factors[i],
          paste(sprintf('"%s"', fac_levels[[i]]), collapse = ", "))
}

# The dodge width of the interaction plot; the Python engine offsets each line
# by the same amount.
.fa_int_dodge <- 0.25

.fa_meta_effect <- function(i, eff) {
  list(effect = eff$effect[i], df1 = eff$df1[i], df2 = eff$df2[i],
       statistic = eff$statistic[i], p_value = eff$p_value[i],
       effect_size = if (is.na(eff$effect_size[i])) NULL else eff$effect_size[i])
}

# The mean, SEM and n of each cell of the design, for the interaction geom.
.fa_cell_summary <- function(df, y, factors) {
  by <- df[factors]
  out <- stats::aggregate(df[y], by = by, FUN = mean)
  names(out)[ncol(out)] <- "mean"
  sds <- stats::aggregate(df[y], by = by, FUN = stats::sd)[[y]]
  ns <- stats::aggregate(df[y], by = by, FUN = length)[[y]]
  out$sem <- sds / sqrt(ns)
  out$n <- ns
  out
}

# Wrap the --x labels on spaces at 12 characters. strwrap keeps a line shorter
# than its width, so the width is 13; the Python engine uses
# textwrap.fill(12, break_on_hyphens = False) and gives the same lines.
.fa_wrap_one <- function(s) paste(strwrap(s, width = 13), collapse = "\n")
.fa_wrap_labels <- function(x) vapply(as.character(x), .fa_wrap_one, character(1), USE.NAMES = FALSE)
.fa_needs_wrap <- function(lv) any(nchar(as.character(lv)) > 12)

# The figure size: wide enough for every --x level times every --fill level in
# each panel, plus the legend; never smaller than the original 5.2 (6.8 with a
# facet) by 4.4 inches. Wrapped labels add a line of height each.
.fa_fig_size <- function(x_lv, fill_lv, fill_name, n_facet, has_facet) {
  slot <- 0.3 * length(fill_lv) + 0.35
  legend_w <- 0.5 + 0.08 * max(nchar(c(as.character(fill_lv), fill_name)))
  width <- max(if (has_facet) 6.8 else 5.2, 0.9 + n_facet * length(x_lv) * slot + legend_w)
  lines <- max(lengths(strsplit(.fa_wrap_labels(x_lv), "\n", fixed = TRUE)))
  list(width = width, height = 4.4 + if (.fa_needs_wrap(x_lv)) 0.18 * (lines - 1) else 0)
}

.fa_methods <- function(eff, spec, factors, use_art) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v    <- pkg_versions(c("car", "ARTool", "ggplot2"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  intro <- if (use_art) {
    sprintf("%s was analyzed by an aligned rank transform ANOVA, a non-parametric factorial method, with factors %s.",
            cap_first(y_label), paste(factors, collapse = ", "))
  } else {
    sprintf("%s was analyzed by a %d-way ANOVA (Type II sums of squares) with factors %s.",
            cap_first(y_label), length(factors), paste(factors, collapse = ", "))
  }
  parts <- vapply(seq_len(nrow(eff)), .fa_effect_part, character(1), eff = eff)
  soft <- if (use_art) sprintf("ARTool %s", v$ARTool) else sprintf("car %s", v$car)
  es_txt <- if (use_art) "" else " Effect sizes are partial eta squared."

  paste0(
    intro, " Effects tested: ", paste(parts, collapse = "; "), ".", es_txt,
    " Significance was set at P < 0.05.",
    sprintf(" Analyses were performed in R %s with %s; figures were produced with ggplot2 %s.",
            rver, soft, v$ggplot2)
  )
}

.fa_effect_part <- function(i, eff) {
  e <- eff[i, ]
  sprintf("%s (F(%s, %s) = %.2f, %s)",
          gsub(":", " x ", e$effect), e$df1, e$df2, e$statistic, .fa_plab(e$p_value))
}

# The x position of each point of the box and bar geoms: the centre of its cell,
# where position_dodge(0.8) puts the box or the bar, plus a fixed golden-ratio
# offset by its row in the cell. The offsets span +/- 0.15 / (m + 2) for m --fill
# levels, the jitter of position_jitterdodge(jitter.width = 0.15), so the points
# keep their places in every run. The Python engine uses the same offsets.
.fa_point_x <- function(df, f1, f2, f3 = NULL) {
  m <- nlevels(df[[f2]])
  a <- 0.15 / (m + 2)
  key <- paste(df[[f1]], df[[f2]], if (is.null(f3) || !nzchar(f3)) "" else df[[f3]], sep = "\r")
  k <- stats::ave(seq_len(nrow(df)), key, FUN = seq_along)
  df$.point_x <- as.integer(df[[f1]]) - 0.4 + 0.8 * (as.integer(df[[f2]]) - 0.5) / m +
    rain_offsets(max(k, 0), -a, a)[k]
  df
}

.fa_emit_script <- function(spec, factors, use_art, pal, width, height, wrap_x, fac_levels,
                            in_name, fig_stub) {
  y  <- spec$data$y
  f1 <- factors[1]; f2 <- factors[2]; f3 <- if (length(factors) >= 3) factors[3] else NULL
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  geom <- spec$appearance$geom %||% "box"
  if (geom == "interaction") {
    return(.fa_emit_interaction_script(spec, factors, use_art, pal_str, width, height, wrap_x,
                                       fac_levels, in_name, fig_stub))
  }

  omni_lines <- .fa_script_model_lines(y, factors, use_art)
  geom_lines <- if (geom == "bar") c(
    '  stat_summary(fun = mean, geom = "bar", position = position_dodge(0.8), width = 0.7, alpha = 0.9, colour = "black", linewidth = 0.3) +',
    '  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(0.8), width = 0.2) +'
  ) else
    '  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.9) +'
  facet_line <- if (!is.null(f3)) sprintf('p <- p + facet_wrap(vars(.data[["%s"]]))', f3) else character(0)
  show_points <- isTRUE(spec$appearance$show_points %||% TRUE)
  point_prep <- if (show_points) c(
    "# Each point sits at its dodged cell centre plus a fixed golden-ratio offset",
    "# by its row in the cell, as in the figure.",
    sprintf('point_m <- nlevels(df[["%s"]]); point_a <- 0.15 / (point_m + 2)', f2),
    sprintf('point_k <- ave(seq_len(nrow(df)), paste(df[["%s"]], df[["%s"]]%s, sep = "|"), FUN = seq_along)',
            f1, f2, if (!is.null(f3)) sprintf(', df[["%s"]]', f3) else ""),
    sprintf('df$.point_x <- as.integer(df[["%s"]]) - 0.4 + 0.8 * (as.integer(df[["%s"]]) - 0.5) / point_m - point_a + 2 * point_a * ((point_k * 0.6180339887498949) %%%% 1)',
            f1, f2)
  ) else character(0)
  point_layer <- if (show_points) {
    "  geom_point(aes(x = .point_x), size = 1.3, alpha = 0.7, shape = 21, stroke = 0.3) +"
  } else character(0)

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(rstatix); library(readr); library(ARTool)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    vapply(seq_along(factors), .fa_factor_line, character(1), factors = factors,
           fac_levels = fac_levels),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "",
    omni_lines,
    'plab <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))',
    'sig <- eff[which(eff$p < 0.05), , drop = FALSE]',
    'subtitle_txt <- if (nrow(sig) > 0) paste(sprintf("%s: %s", gsub(":", " × ", sig$effect), vapply(sig$p, plab, character(1))), collapse = "\\n") else "No significant effects (P ≥ 0.05)"',
    "",
    point_prep,
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]], fill = .data[["%s"]])) +', f1, y, f2),
    geom_lines,
    point_layer,
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  labs(x = "%s", y = "%s", fill = "%s", subtitle = subtitle_txt) +',
            spec$appearance$x_label %||% f1, spec$appearance$y_label %||% y, f2),
    "  theme_classic(base_size = 13) +",
    '  theme(plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey30"))',
    facet_line,
    if (isTRUE(wrap_x)) "p <- p + scale_x_discrete(labels = scales::label_wrap(13))" else character(0),
    "",
    sprintf('ggsave("%s.pdf", p, width = %s, height = %s)', fig_stub,
            sprintf("%.4g", width), sprintf("%.4g", height)),
    'cat("Reproduced figure.\\n")'
  )
}

# The lines of an emitted script that fit the model of the recipe on a copy with
# syntactic column names and give the effects their column names back.
.fa_script_model_lines <- function(y, factors, use_art) {
  safe <- sprintf(".f%d", seq_along(factors))
  map <- paste(sprintf('"%s" = "%s"', safe, factors), collapse = ", ")
  cols <- paste(sprintf('"%s"', factors), collapse = ", ")
  fml <- sprintf(".y ~ %s", paste(safe, collapse = " * "))
  model <- if (use_art) c(
    sprintf("at <- as.data.frame(anova(ARTool::art(%s, data = sdf)))", fml),
    'eff <- data.frame(effect = as.character(at$Term), p = at[["Pr(>F)"]])'
  ) else c(
    sprintf("at <- as.data.frame(car::Anova(lm(%s, data = sdf), type = 2))", fml),
    'at <- at[rownames(at) != "Residuals", ]',
    'eff <- data.frame(effect = rownames(at), p = at[["Pr(>F)"]])'
  )
  c("# The models run on a copy with syntactic column names; the effects take the",
    "# names of the columns back.",
    sprintf("safe <- c(%s)", map),
    sprintf('sdf <- setNames(data.frame(df[["%s"]], df[c(%s)]), c(".y", names(safe)))', y, cols),
    model,
    'term_name <- function(t) paste(safe[strsplit(t, ":", fixed = TRUE)[[1]]], collapse = ":")',
    "eff$effect <- vapply(eff$effect, term_name, character(1), USE.NAMES = FALSE)")
}

# The standalone script of the interaction geom: the same ANOVA, then the mean
# +/- SEM of each cell joined by one line per --fill level.
.fa_emit_interaction_script <- function(spec, factors, use_art, pal_str, width, height, wrap_x,
                                        fac_levels, in_name, fig_stub) {
  y  <- spec$data$y
  f1 <- factors[1]; f2 <- factors[2]; f3 <- if (length(factors) >= 3) factors[3] else NULL
  omni_lines <- .fa_script_model_lines(y, factors, use_art)
  by_line <- sprintf('by <- df[c(%s)]', paste(sprintf('"%s"', factors), collapse = ", "))
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(readr); library(ARTool)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    vapply(seq_along(factors), .fa_factor_line, character(1), factors = factors,
           fac_levels = fac_levels),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "",
    omni_lines,
    'plab <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))',
    'sig <- eff[which(eff$p < 0.05), , drop = FALSE]',
    'subtitle_txt <- if (nrow(sig) > 0) paste(sprintf("%s: %s", gsub(":", " × ", sig$effect), vapply(sig$p, plab, character(1))), collapse = "\\n") else "No significant effects (P ≥ 0.05)"',
    "",
    "# The mean and SEM of each cell.",
    by_line,
    sprintf('summ <- aggregate(df["%s"], by = by, FUN = mean); names(summ)[ncol(summ)] <- "mean"', y),
    sprintf('summ$sem <- aggregate(df["%s"], by = by, FUN = sd)[["%s"]] / sqrt(aggregate(df["%s"], by = by, FUN = length)[["%s"]])', y, y, y, y),
    "",
    "pd <- position_dodge(width = 0.25)",
    sprintf('p <- ggplot(summ, aes(x = .data[["%s"]], y = mean, colour = .data[["%s"]], group = .data[["%s"]])) +', f1, f2, f2),
    "  geom_line(position = pd, linewidth = 0.8) +",
    "  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.12, position = pd, linewidth = 0.6) +",
    "  geom_point(position = pd, size = 2.4) +",
    sprintf("  scale_colour_manual(values = %s) +", pal_str),
    sprintf('  labs(x = "%s", y = "%s", colour = "%s", subtitle = subtitle_txt) +',
            spec$appearance$x_label %||% f1, spec$appearance$y_label %||% y, f2),
    "  theme_classic(base_size = 13) +",
    '  theme(plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey30"))',
    if (!is.null(f3)) sprintf('p <- p + facet_wrap(vars(.data[["%s"]]))', f3) else character(0),
    if (isTRUE(wrap_x)) "p <- p + scale_x_discrete(labels = scales::label_wrap(13))" else character(0),
    "",
    sprintf('ggsave("%s.pdf", p, width = %s, height = %s)', fig_stub,
            sprintf("%.4g", width), sprintf("%.4g", height)),
    'cat("Reproduced figure.\\n")'
  )
}
