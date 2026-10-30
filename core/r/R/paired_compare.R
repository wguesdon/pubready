# Recipe: paired before/after comparison. Each subject (--id) is measured under
# two conditions (--x); a line connects the subject's two points. The paired
# t-test or the Wilcoxon signed-rank test is chosen from the normality of the
# within-subject differences, unless forced with --test.

# The --test values that paired_compare knows (F63).
.paired_compare_tests <- c("auto", "t", "t_test", "paired_t", "student_t", "wilcoxon", "wilcox",
                           "signed_rank", "mwu")

recipe_paired_compare <- function(df, spec) {
  check_test_method("paired_compare", spec$test$method, .paired_compare_tests)
  x <- spec$data$x; y <- spec$data$y; idc <- spec$data$id
  if (is.null(x) || is.null(y)) stop("paired_compare needs --x (condition) and --y (value)")
  if (is.null(idc) || !nzchar(idc)) stop("paired_compare needs --id (the subject column)")
  for (cn in c(idc, x, y)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  df <- df[, c(idc, x, y)]
  df[[y]] <- suppressWarnings(as.numeric(df[[y]]))
  n_rows0 <- nrow(df)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  # One level rule in both engines: the --levels order (--levels pre,post keeps
  # pre left of post), else numeric order when every label is a number, else
  # code-point order. The row order of the file no longer sets the sign.
  lab <- level_labels(df[[x]])
  lv <- level_order(lab, unlist(spec$data$levels), "condition")
  df[[x]] <- factor(lab, levels = lv)
  if (length(lv) != 2) {
    stop(sprintf("paired_compare needs exactly 2 conditions in '%s'; found %d (%s)",
                 x, length(lv), paste(lv, collapse = ", ")))
  }

  # One value per subject and condition: a second row would make the result
  # depend on the row order, so the run stops and names the subjects.
  dup <- duplicated(data.frame(level_labels(df[[idc]]), df[[x]]))
  if (any(dup)) {
    stop(sprintf("paired_compare needs one value per subject and condition; '%s' repeats %s within a condition",
                 idc, .paired_subject_list(unique(level_labels(df[[idc]][dup])))))
  }
  a <- df[df[[x]] == lv[1], c(idc, y)]; names(a) <- c("id", "v1")
  b <- df[df[[x]] == lv[2], c(idc, y)]; names(b) <- c("id", "v2")
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
  # When the differences fail the normality check but the signed-rank test
  # cannot reach P < 0.05 with this count of pairs, auto keeps the paired t-test.
  rank_min_p <- rank_test_min_p(nrow(m), "signed_rank")
  kept <- method %in% c("auto", "") && !isTRUE(norm$normal) &&
    rank_test_cannot_reach(rank_min_p) && !is_zero_variance(diffs)
  family <- if (method %in% c("auto", "")) {
    if (isTRUE(norm$normal) || kept) "t" else "wilcoxon"
  } else if (method %in% c("t", "t_test", "paired_t", "student_t")) "t"
    else if (method %in% c("wilcoxon", "wilcox", "signed_rank", "mwu")) "wilcoxon"
    else stop(sprintf("unknown paired test '%s' (use t or wilcoxon)", method))
  if (family == "t" && is_zero_variance(diffs)) {
    stop(zero_variance_message("paired t-test", "the within-subject differences have",
                               "use --test wilcoxon"))
  }

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
  # The same star scale as every other recipe and the Python engine: **** for
  # P < 1e-4.
  star <- p_stars(pval)
  blabel <- spec$appearance$bracket$label %||% "p.signif"
  brk <- if (blabel == "p.format") paste0("p = ", formatC(pval, format = "g", digits = 2)) else star
  ypos <- max(plong$v, na.rm = TRUE) + 0.06 * diff(range(plong$v, na.rm = TRUE))

  p <- ggplot2::ggplot(plong, ggplot2::aes(x = cond, y = v)) +
    ggplot2::geom_line(ggplot2::aes(group = id), colour = "grey70",
                       linewidth = 0.4, alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(colour = cond), size = 2.4, alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_x_discrete(labels = stats::setNames(wrap_labels(lv), lv)) +
    ggplot2::scale_y_continuous(labels = plain_number_labels) +
    ggpubr::geom_bracket(xmin = lv[1], xmax = lv[2], y.position = ypos,
                         label = brk, tip.length = 0.02) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% y,
                  title = spec$appearance$title) +
    pub_theme(spec$appearance$theme)

  stats_df <- data.frame(
    outcome = y, condition1 = lv[1], condition2 = lv[2], contrast = paste(lv[2], "-", lv[1]),
    n_pairs = nrow(m),
    test = label, statistic = stat, df = dfree, p_value = pval,
    mean_diff = mean(diffs), effect_size_name = eff_name, effect_size = eff,
    shapiro_p_diff = norm$p, stringsAsFactors = FALSE)

  test_meta <- list(
    name = label, statistic = stat, p_value = pval,
    n_pairs = nrow(m), mean_difference = mean(diffs),
    effect_size = list(name = eff_name, value = eff),
    assumptions = list(normality_of_differences = list(test = "shapiro-wilk", p = norm$p)),
    rank_test_attainable_p = rank_min_p)

  resolved <- list(family = family, label = label)
  floor_txt <- if (kept || family == "wilcoxon") rank_floor_text(rank_min_p, kept) else ""
  methods  <- .paired_methods(label, lv, nrow(m), pval, eff_name, eff, spec, floor_txt)
  build_script <- list(fn = .paired_emit_script, args = list(spec = spec, family = family, lv = lv))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(y), df_used = m, clean_steps = clean_steps,
       qc = list(quantity = "within-subject differences",
                 panels = list(differences = diffs)),
       width = 3.8, height = 4.2)
}

# The subjects of a message: "the subject 'S02'", or the first ten subjects and
# the count of the others. The Python engine writes the same words.
.paired_subject_list <- function(ids) {
  shown <- utils::head(ids, 10)
  txt <- sprintf("the subject%s %s", if (length(ids) > 1) "s" else "", quote_join(shown))
  if (length(ids) > 10) txt <- sprintf("%s and %d more", txt, length(ids) - 10)
  txt
}

.paired_methods <- function(label, lv, n, pval, eff_name, eff, spec, floor_txt = "") {
  v    <- pkg_versions(c("ggplot2", "ggpubr"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  y_lab <- spec$appearance$y_label %||% spec$data$y
  paste0(
    sprintf("%s was compared between %s and %s within each subject with a %s (n = %d pairs). ",
            cap_first(y_lab), lv[1], lv[2], label, n),
    if (tolower(spec$test$method %||% "auto") == "auto") {
      "The test was chosen from the Shapiro-Wilk normality of the within-subject differences. "
    } else "The test was set by the user. ",
    if (nzchar(floor_txt)) paste0(floor_txt, " ") else "",
    sprintf("Effect size is reported as %s = %s (%s minus %s). Significance was set at P < 0.05. ",
            eff_name, formatC(eff, format = "f", digits = 2), lv[2], lv[1]),
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
