# Recipe: compare proportions of a categorical outcome across groups. Input is
# raw rows with a group column (--x) and an outcome column (--y). A chi-squared
# test of independence is used, or Fisher's exact test when any expected cell is
# below 5 (or when forced). For a 2x2 table an odds ratio is reported. The figure
# is a 100% stacked bar of the outcome proportions within each group.

# The --test values that proportions knows (F63).
.proportions_tests <- c("auto", "chisq", "chi", "chi_square", "chisquare", "fisher", "fisher_exact")

recipe_proportions <- function(df, spec) {
  check_test_method("proportions", spec$test$method, .proportions_tests)
  x <- spec$data$x; y <- spec$data$y
  if (is.null(x) || is.null(y)) stop("proportions needs --x (group) and --y (outcome)")
  for (cn in c(x, y)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  df <- df[, c(x, y)]
  n0 <- nrow(df)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  # The labels of both engines: a whole number reads "6", not "6.0".
  df[[x]] <- level_labels(df[[x]]); df[[y]] <- level_labels(df[[y]])
  df[[x]] <- factor(df[[x]], levels = .prop_group_order(df[[x]], unlist(spec$data$levels)))
  df[[y]] <- factor(df[[y]], levels = unique(as.character(df[[y]])))
  clean_steps <- if (nrow(df) < n0) {
    sprintf("Dropped %d row(s) with a missing group or outcome (%d -> %d).", n0 - nrow(df), n0, nrow(df))
  } else "No cleaning applied; input used as-is."

  tab <- table(df[[x]], df[[y]])
  if (nrow(tab) < 2 || ncol(tab) < 2) stop("proportions needs at least 2 groups and 2 outcomes")

  method  <- tolower(spec$test$method %||% "auto")
  small   <- any(suppressWarnings(stats::chisq.test(tab)$expected) < 5)
  use_fisher <- if (method %in% c("auto", "")) small else method %in% c("fisher", "fisher_exact")
  if (method %in% c("chisq", "chi", "chi_square", "chisquare")) use_fisher <- FALSE

  is2x2 <- all(dim(tab) == 2)
  or <- NA_real_; or_ci <- c(NA_real_, NA_real_)
  if (use_fisher) {
    ft <- stats::fisher.test(tab)
    pval <- ft$p.value; stat <- NA_real_; dfree <- NA_real_
    test_label <- "Fisher's exact test"
    if (is2x2) { or <- unname(ft$estimate); or_ci <- as.numeric(ft$conf.int) }
  } else {
    ct <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
    pval <- ct$p.value; stat <- unname(ct$statistic); dfree <- unname(ct$parameter)
    test_label <- "Pearson's chi-squared test"
    # A zero cell leaves the sample odds ratio at 0 or infinity, so it is not
    # given (NA), and the methods text says why. The Python engine does the same.
    if (is2x2 && all(tab > 0)) or <- (tab[1, 1] * tab[2, 2]) / (tab[1, 2] * tab[2, 1])
  }

  prop <- prop.table(tab, margin = 1)
  cnt  <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(cnt) <- c("group", "outcome", "n")
  cnt$proportion <- as.data.frame(prop, stringsAsFactors = FALSE)$Freq
  cnt$group   <- factor(cnt$group, levels = levels(df[[x]]))
  cnt$outcome <- factor(cnt$outcome, levels = levels(df[[y]]))

  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = nlevels(cnt$outcome))
  psub <- sprintf("%s, P = %s", test_label, formatC(pval, format = "g", digits = 2))

  p <- ggplot2::ggplot(cnt, ggplot2::aes(x = group, y = proportion, fill = outcome)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(proportion >= 0.04,
                                                   scales::percent(proportion, accuracy = 1), "")),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3, colour = "white") +
    ggplot2::scale_x_discrete(labels = stats::setNames(wrap_labels(levels(cnt$group)), levels(cnt$group))) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(x = spec$appearance$x_label %||% x,
                  y = spec$appearance$y_label %||% "proportion",
                  fill = y, title = spec$appearance$title, subtitle = psub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  stats_df <- cnt[, c("group", "outcome", "n", "proportion")]
  stats_df$omnibus_test <- test_label
  stats_df$omnibus_statistic <- stat
  stats_df$omnibus_df1 <- dfree
  stats_df$omnibus_df2 <- NA_real_
  stats_df$omnibus_p <- pval
  stats_df$odds_ratio <- or
  test_meta <- list(name = test_label, statistic = stat, df = dfree, p_value = pval,
                    odds_ratio = or, or_conf_low = or_ci[1], or_conf_high = or_ci[2],
                    n = sum(tab), any_expected_below_5 = small)
  resolved <- list(family = "proportions", label = test_label)
  methods  <- .prop_methods(test_label, nrow(tab), ncol(tab), sum(tab), pval, or, is2x2,
                            stat, dfree, rownames(tab), colnames(tab),
                            zero_cell = is2x2 && !use_fisher && any(tab == 0))
  build_script <- list(fn = .prop_emit_script, args = list(spec = spec, use_fisher = use_fisher))

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(y), df_used = df, clean_steps = clean_steps,
       width = max(3.6, category_width(nrow(tab), wrap_labels(levels(cnt$group)), 0) +
                          legend_width(levels(cnt$outcome), y)),
       height = 4.2)
}

.prop_methods <- function(test_label, n_grp, n_out, n, pval, or, is2x2, stat, dfree, groups, outcomes,
                          zero_cell = FALSE) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  # The odds ratio is the odds of the first outcome in the first group over the
  # odds in the second group, so the text names the outcome and both groups.
  or_txt <- if (is2x2 && is.finite(or)) {
    sprintf(" The odds ratio of %s for %s versus %s was %s.", outcomes[1], groups[1], groups[2],
            formatC(or, format = "f", digits = 2))
  } else if (zero_cell) {
    " The odds ratio is not given, because a cell of the 2 x 2 table is zero."
  } else ""
  p_txt <- sprintf("P = %s", formatC(pval, format = "g", digits = 2))
  if (is.finite(stat)) p_txt <- sprintf("chi-square(%d) = %.2f, %s", as.integer(dfree), stat, p_txt)
  paste0(
    sprintf("The distribution of the outcome (%d categories) across %d groups (n = %d) was tested with %s (%s).",
            n_out, n_grp, n, test_label, p_txt),
    or_txt,
    " Proportions are shown as a 100% stacked bar per group.",
    sprintf(" Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.prop_emit_script <- function(spec, use_fisher, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  test_line <- if (use_fisher) "print(fisher.test(tab))" else "print(chisq.test(tab, correct = FALSE))"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container.",
    "suppressPackageStartupMessages({ library(ggplot2); library(scales) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('raw <- raw[complete.cases(raw[, c("%s","%s")]), ]', x, y),
    sprintf('tab <- table(raw[["%s"]], raw[["%s"]])', x, y),
    test_line,
    "pd <- as.data.frame(prop.table(tab, 1)); names(pd) <- c(\"group\", \"outcome\", \"prop\")",
    "p <- ggplot(pd, aes(group, prop, fill = outcome)) + geom_col(width = 0.7, colour = \"white\") +",
    "  geom_text(aes(label = percent(prop, accuracy = 1)), position = position_stack(vjust = 0.5), colour = \"white\") +",
    sprintf('  scale_y_continuous(labels = percent) + labs(x = "%s", y = "proportion") + theme_classic()', x),
    sprintf('ggsave("%s.pdf", p, width = 4.2, height = 4.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}

# Group order: the order of --levels when it is given, else first appearance in
# the data. A group of the data that --levels does not list stops the run; a
# listed group that no row holds is left out.
.prop_group_order <- function(g, levels = NULL) {
  u <- unique(as.character(g))
  if (length(levels) == 0) return(u)
  levels <- unique(trimws(as.character(levels)))
  miss <- setdiff(u, levels)
  if (length(miss) > 0) {
    stop(sprintf("--levels does not list the group(s) %s of the data",
                 paste(sprintf("'%s'", miss), collapse = ", ")))
  }
  levels[levels %in% u]
}
