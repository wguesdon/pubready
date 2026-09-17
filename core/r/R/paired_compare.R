# Recipe: paired before/after comparison. Each subject (--id) is measured under
# two conditions (--x); a line connects the subject's two points. The paired
# t-test or the Wilcoxon signed-rank test is chosen from the normality of the
# within-subject differences, unless forced with --test.

recipe_paired_compare <- function(df, spec) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id
  if (is.null(x) || is.null(y)) stop("paired_compare needs --x (condition) and --y (value)")
  if (is.null(idc) || !nzchar(idc)) stop("paired_compare needs --id (the subject column)")
  for (cn in c(idc, x, y)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  df <- df[, c(idc, x, y)]
  df[[y]] <- suppressWarnings(as.numeric(df[[y]]))
  n_rows0 <- nrow(df)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  # Preserve the order the conditions appear in (so pre stays left of post),
  # rather than sorting alphabetically.
  df[[x]] <- factor(df[[x]], levels = unique(as.character(df[[x]])))
  lv <- levels(df[[x]])
  if (length(lv) != 2) {
    stop(sprintf("paired_compare needs exactly 2 conditions in '%s'; found %d (%s)",
                 x, length(lv), paste(lv, collapse = ", ")))
  }

  a <- df[df[[x]] == lv[1], c(idc, y)]; names(a) <- c("id", "v1")
  b <- df[df[[x]] == lv[2], c(idc, y)]; names(b) <- c("id", "v2")
  a <- a[!duplicated(a$id), ]; b <- b[!duplicated(b$id), ]
  m <- merge(a, b, by = "id")
  n_ids <- length(unique(df[[idc]]))
  dropped <- n_ids - nrow(m)
  if (nrow(m) < 2) stop("paired_compare needs at least 2 complete pairs")
  clean_steps <- c(
    sprintf("Kept %d complete pair(s) across the two conditions.", nrow(m)),
    if (dropped > 0) sprintf("Dropped %d subject(s) missing one condition.", dropped) else NULL)

  diffs <- m$v2 - m$v1
  norm  <- shapiro_safe(diffs)
  method <- tolower(spec$test$method %||% "auto")
  family <- if (method %in% c("auto", "")) {
    if (isTRUE(norm$normal)) "t" else "wilcoxon"
  } else if (method %in% c("t", "t_test", "paired_t", "student_t")) "t"
    else if (method %in% c("wilcoxon", "wilcox", "signed_rank", "mwu")) "wilcoxon"
    else stop(sprintf("unknown paired test '%s' (use t or wilcoxon)", method))

  if (family == "t") {
    tt <- stats::t.test(m$v2, m$v1, paired = TRUE)
    stat <- unname(tt$statistic); pval <- tt$p.value; dfree <- unname(tt$parameter)
    eff <- mean(diffs) / stats::sd(diffs); eff_name <- "Cohen's dz"; label <- "paired t-test"
  } else {
    # exact = NULL lets R take the exact route when no tie and no zero remain,
    # which is the rule the Python engine repeats.
    wt <- suppressWarnings(stats::wilcox.test(m$v2, m$v1, paired = TRUE))
    stat <- unname(wt$statistic); pval <- wt$p.value; dfree <- NA_real_
    eff <- rank_biserial_paired(diffs)
    eff_name <- "rank-biserial r"; label <- "Wilcoxon signed-rank test"
  }

  plong <- rbind(
    data.frame(id = m$id, cond = lv[1], v = m$v1, stringsAsFactors = FALSE),
    data.frame(id = m$id, cond = lv[2], v = m$v2, stringsAsFactors = FALSE))
  plong$cond <- factor(plong$cond, levels = lv)

  pal <- spec$appearance$palette %||% pubready_palette[seq_len(2)]
  star <- if (is.na(pval)) "ns" else if (pval < 0.001) "***" else
          if (pval < 0.01) "**" else if (pval < 0.05) "*" else "ns"
  blabel <- spec$appearance$bracket$label %||% "p.signif"
  brk <- if (blabel == "p.format") paste0("p = ", formatC(pval, format = "g", digits = 2)) else star
  ypos <- max(plong$v, na.rm = TRUE) + 0.06 * diff(range(plong$v, na.rm = TRUE))

  p <- ggplot2::ggplot(plong, ggplot2::aes(x = cond, y = v)) +
    ggplot2::geom_line(ggplot2::aes(group = id), colour = "grey70",
                       linewidth = 0.4, alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(colour = cond), size = 2.4, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal) +
    ggpubr::geom_bracket(xmin = lv[1], xmax = lv[2], y.position = ypos,
                         label = brk, tip.length = 0.02) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% y,
                  title = spec$appearance$title) +
    pub_theme(spec$appearance$theme)

  stats_df <- data.frame(
    outcome = y, condition1 = lv[1], condition2 = lv[2], n_pairs = nrow(m),
    test = label, statistic = stat, df = dfree, p_value = pval,
    mean_diff = mean(diffs), effect_size_name = eff_name, effect_size = eff,
    shapiro_p_diff = norm$p, stringsAsFactors = FALSE)

  test_meta <- list(
    name = label, statistic = stat, p_value = pval,
    n_pairs = nrow(m), mean_difference = mean(diffs),
    effect_size = list(name = eff_name, value = eff),
    assumptions = list(normality_of_differences = list(test = "shapiro-wilk", p = norm$p)))

  resolved <- list(family = family, label = label)
  methods  <- .paired_methods(label, lv, nrow(m), pval, eff_name, eff, spec)
  build_script <- function(in_name, fig_stub) .paired_emit_script(spec, family, lv, in_name, fig_stub)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(y), df_used = m, clean_steps = clean_steps,
       qc = list(quantity = "within-subject differences",
                 panels = list(differences = diffs)),
       width = 3.8, height = 4.2)
}

.paired_methods <- function(label, lv, n, pval, eff_name, eff, spec) {
  v    <- pkg_versions(c("ggplot2", "ggpubr"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  y_lab <- spec$appearance$y_label %||% spec$data$y
  paste0(
    sprintf("%s was compared between %s and %s within each subject with a %s (n = %d pairs). ",
            cap_first(y_lab), lv[1], lv[2], label, n),
    "The test was chosen from the Shapiro-Wilk normality of the within-subject differences. ",
    sprintf("Effect size is reported as %s = %s. Significance was set at P < 0.05. ",
            eff_name, formatC(eff, format = "f", digits = 2)),
    sprintf("Drawn in R %s with ggplot2 %s and ggpubr %s.", rver, v$ggplot2, v$ggpubr))
}

.paired_emit_script <- function(spec, family, lv, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id
  test_line <- if (family == "t") {
    "res <- t.test(m$v2, m$v1, paired = TRUE)"
  } else "res <- suppressWarnings(wilcox.test(m$v2, m$v1, paired = TRUE))"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container.",
    "suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('raw[["%s"]] <- as.numeric(raw[["%s"]]); raw <- raw[complete.cases(raw[, c("%s","%s","%s")]), ]',
            y, y, idc, x, y),
    sprintf('a <- raw[raw[["%s"]] == "%s", c("%s","%s")]; names(a) <- c("id","v1")', x, lv[1], idc, y),
    sprintf('b <- raw[raw[["%s"]] == "%s", c("%s","%s")]; names(b) <- c("id","v2")', x, lv[2], idc, y),
    "m <- merge(a[!duplicated(a$id), ], b[!duplicated(b$id), ], by = \"id\")",
    test_line,
    "print(res)",
    sprintf('pl <- rbind(data.frame(id=m$id, cond="%s", v=m$v1), data.frame(id=m$id, cond="%s", v=m$v2))',
            lv[1], lv[2]),
    sprintf('pl$cond <- factor(pl$cond, levels = c("%s","%s"))', lv[1], lv[2]),
    "p <- ggplot(pl, aes(cond, v)) + geom_line(aes(group = id), colour = \"grey70\") +",
    "  geom_point(aes(colour = cond), size = 2.4) +",
    sprintf('  labs(x = "%s", y = "%s") + theme_classic() + theme(legend.position = "none")', x, y),
    sprintf('ggsave("%s.pdf", p, width = 3.8, height = 4.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
