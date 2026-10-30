# Recipe: the composition of each sample, for example the cell-type proportions
# of each sample of a single-cell experiment, compared between conditions.
#
# Input is long: one row per sample (--id) and category (--x) with a count
# (--y), and an optional condition (--group). Without --y each row is one
# observation and the recipe counts the rows. The counts are summed per sample
# and category and turned into proportions within each sample.
#
# The figure is one 100% stacked bar per sample, the samples ordered by group
# and then by name, one strip per group. With two groups each category is
# compared with the Mann-Whitney U test on the per-sample proportions, with three
# or more with the Kruskal-Wallis test, and the P values are adjusted across the
# categories with the Benjamini-Hochberg method. Proportions are compositional
# (one category rising forces the others down), so the tests are exploratory.

.comp_slot_in <- function(n_samples) if (n_samples <= 20) 0.3 else 0.22

# Wrap a label at a width, the same lines as the Python engine's
# textwrap.fill(break_on_hyphens = False): strwrap keeps lines shorter than its
# width, so the width is one more than the longest line wanted.
.comp_wrap <- function(s, width) paste(strwrap(s, width = width + 1), collapse = "\n")

# Natural order of the sample names: the text before a trailing number, then
# the number, then the full name, so that S2 comes before S10.
.comp_natural_order <- function(s) {
  has_num <- grepl("[0-9]+$", s)
  prefix <- sub("[0-9]+$", "", s)
  num <- rep(NA_real_, length(s))
  num[has_num] <- as.numeric(sub("^.*?([0-9]+)$", "\\1", s[has_num], perl = TRUE))
  order(prefix, num, s, method = "radix")
}

.comp_signif <- function(p) {
  if (is.na(p)) "" else if (p < 1e-4) "****" else if (p < 1e-3) "***"
  else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

# The test of one category on the per-sample proportions. A category whose
# proportion is the same in every sample has no test.
.comp_test_category <- function(j, P, g, lv) {
  v <- P[, j]
  if (length(unique(v)) < 2) return(c(NA_real_, NA_real_))
  if (length(lv) == 2) {
    wt <- suppressWarnings(stats::wilcox.test(v[g == lv[1]], v[g == lv[2]]))
    c(unname(wt$statistic), wt$p.value)
  } else {
    kt <- stats::kruskal.test(v, factor(g, levels = lv))
    c(unname(kt$statistic), kt$p.value)
  }
}

# The rows of the stats table for one category: one row per group.
.comp_category_rows <- function(j, P, g, lv, cats, tests, test_label) {
  out <- NULL
  for (gl in lv) {
    v <- P[g == gl, j]
    out <- rbind(out, data.frame(
      category = cats[j], group = gl, n_samples = length(v),
      mean_proportion = mean(v),
      sd_proportion = if (length(v) >= 2) stats::sd(v) else NA_real_,
      median_proportion = stats::median(v),
      test = test_label, statistic = tests$statistic[j], p_value = tests$p_value[j],
      p_adj = tests$p_adj[j], significant = isTRUE(tests$p_adj[j] < 0.05),
      stringsAsFactors = FALSE))
  }
  out
}

recipe_composition <- function(df, spec) {
  id <- spec$data$id
  cat_col <- spec$data$x
  cnt <- spec$data$y
  grp <- spec$data$group
  has_count <- !is.null(cnt) && nzchar(cnt)
  has_group <- !is.null(grp) && nzchar(grp)
  if (is.null(id) || !nzchar(id) || is.null(cat_col) || !nzchar(cat_col)) {
    stop("composition needs --id (the sample) and --x (the category); --y (a count) and --group are optional")
  }
  for (cn in c(id, cat_col, if (has_count) cnt, if (has_group) grp)) {
    if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))
  }

  steps <- character(0)
  n0 <- nrow(df)
  d <- data.frame(id = as.character(df[[id]]), category = as.character(df[[cat_col]]),
                  group = if (has_group) as.character(df[[grp]]) else "all",
                  stringsAsFactors = FALSE)
  if (has_count) {
    yv <- suppressWarnings(as.numeric(df[[cnt]]))
    coerced <- sum(is.na(yv) & !is.na(df[[cnt]]))
    if (coerced > 0) steps <- c(steps, sprintf("Coerced column '%s' to numeric; %d value(s) became NA.", cnt, coerced))
    d$count <- yv
  } else {
    d$count <- 1
    steps <- c(steps, "No count column was given, so each row counts as one observation.")
  }
  keep <- !is.na(d$id) & !is.na(d$category) & !is.na(d$group) & is.finite(d$count)
  if (sum(!keep) > 0) {
    steps <- c(steps, sprintf("Dropped %d row(s) with a missing sample, category, group or count (%d -> %d rows).",
                              sum(!keep), n0, sum(keep)))
  }
  d <- d[keep, , drop = FALSE]
  if (any(d$count < 0)) stop("composition needs counts of zero or more")
  if (has_group) {
    per_sample <- unique(d[, c("id", "group")])
    if (anyDuplicated(per_sample$id)) {
      stop(sprintf("each sample must belong to one group; sample '%s' has more than one",
                   per_sample$id[anyDuplicated(per_sample$id)]))
    }
  }

  cats <- unique(d$category)
  ids <- unique(d$id)
  # The summed counts, samples in rows and categories in columns.
  M <- unclass(stats::xtabs(count ~ factor(id, levels = ids) + factor(category, levels = cats), data = d))
  M <- matrix(as.numeric(M), nrow = length(ids))
  totals <- rowSums(M)
  if (any(totals == 0)) {
    steps <- c(steps, sprintf("Dropped %d sample(s) with a total count of zero.", sum(totals == 0)))
    M <- M[totals > 0, , drop = FALSE]
    ids <- ids[totals > 0]
    totals <- totals[totals > 0]
  }
  # A category with no count in any sample has no proportion to test: its P
  # would be NA and it would still count in the BH adjustment, so it is dropped
  # and logged (F88). The Python engine has the same rule.
  empty <- colSums(M) == 0
  if (any(empty)) {
    steps <- c(steps, sprintf("Dropped %d category(ies) with a total count of zero: %s.",
                              sum(empty), .pca_name_list(cats[empty])))
    M <- M[, !empty, drop = FALSE]
    cats <- cats[!empty]
  }
  if (length(cats) < 2) stop("composition needs at least 2 categories")
  if (length(ids) < 1) stop("composition needs at least 1 sample with a count above zero")
  g_of <- d$group[match(ids, d$id)]
  lv <- sort(unique(g_of), method = "radix")
  ord <- integer(0)
  for (gl in lv) {
    in_g <- which(g_of == gl)
    ord <- c(ord, in_g[.comp_natural_order(ids[in_g])])
  }
  ids <- ids[ord]; g_of <- g_of[ord]; M <- M[ord, , drop = FALSE]; totals <- totals[ord]
  P <- M / totals
  if (length(steps) == 0) steps <- "No cleaning applied; input used as-is."

  n_by <- as.integer(table(factor(g_of, levels = lv)))
  run_test <- has_group && length(lv) >= 2 && all(n_by >= 2)
  test_label <- if (!run_test) "none (descriptive)"
                else if (length(lv) == 2) "Mann-Whitney U test" else "Kruskal-Wallis test"
  tests <- data.frame(statistic = rep(NA_real_, length(cats)), p_value = NA_real_, p_adj = NA_real_)
  if (run_test) {
    res <- vapply(seq_along(cats), .comp_test_category, numeric(2), P = P, g = g_of, lv = lv)
    tests$statistic <- res[1, ]
    tests$p_value <- res[2, ]
    tests$p_adj <- stats::p.adjust(tests$p_value, method = "BH")
  }
  stats_df <- do.call(rbind, lapply(seq_along(cats), .comp_category_rows, P = P, g = g_of, lv = lv,
                                    cats = cats, tests = tests, test_label = test_label))
  rownames(stats_df) <- NULL

  dims <- .comp_dims(ids, g_of, lv, cats, tests, spec)
  p <- .comp_plot(P, ids, g_of, lv, cats, tests, run_test, has_group, test_label, dims, spec, id, cat_col)

  long <- data.frame(sample = rep(ids, times = length(cats)), group = rep(g_of, times = length(cats)),
                     category = rep(cats, each = length(ids)), count = as.vector(M),
                     proportion = as.vector(P), stringsAsFactors = FALSE)
  test_meta <- list(
    name = if (run_test) sprintf("%s per category with Benjamini-Hochberg adjustment", test_label)
           else "descriptive (no test)",
    n_samples = length(ids), n_categories = length(cats),
    groups = if (has_group) stats::setNames(as.list(n_by), lv) else NULL,
    p_adjust = if (run_test) "BH" else "none",
    n_significant = sum(tests$p_adj < 0.05, na.rm = TRUE),
    categories = stats::setNames(lapply(seq_along(cats), .comp_meta_row, tests = tests), cats))
  resolved <- list(family = "composition", label = test_meta$name)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .comp_methods(ids, cats, lv, n_by, has_group, run_test, test_label, tests, spec, id, cat_col),
       build_script = list(fn = .comp_emit_script, args = list(spec = spec, has_count = has_count,
                                                                has_group = has_group)),
       label = slugify(cat_col), df_used = long, clean_steps = steps,
       width = dims$width, height = dims$height)
}

.comp_meta_row <- function(j, tests) {
  list(statistic = tests$statistic[j], p_value = tests$p_value[j], p_adj = tests$p_adj[j])
}

# The legend label of one category: wrapped at 14 characters, with the
# significance stars of its adjusted P value.
.comp_legend_label <- function(j, cats, tests) {
  stars <- .comp_signif(tests$p_adj[j])
  paste0(.comp_wrap(cats[j], 14), if (nzchar(stars)) paste0(" ", stars) else "")
}

# The figure size and the label layout, from the sample count, the label
# lengths and the legend height. The Python engine uses the same rules.
.comp_dims <- function(ids, g_of, lv, cats, tests, spec) {
  n <- length(ids)
  slot <- .comp_slot_in(n)
  rotate <- n > 6 || max(nchar(ids)) > 8
  lab_h <- if (rotate) 0.075 * min(max(nchar(ids)), 40) else 0.25
  leg_labels <- vapply(seq_along(cats), .comp_legend_label, character(1), cats = cats, tests = tests)
  leg_lines <- sum(nchar(gsub("[^\n]", "", leg_labels))) + length(cats)
  leg_h <- 0.19 * leg_lines + 0.08 * length(cats) + 0.4
  height <- max(4.4, 3.4 + lab_h + if (n > 20) 0.8 else 0)
  leg_cols <- if (leg_h > height - 0.6) 2 else 1
  if (leg_cols == 2) height <- max(height, leg_h / 2 + 1.0)
  leg_w <- leg_cols * 1.45 + 0.3
  width <- max(4.8, 1.0 + slot * n + 0.12 * (length(lv) - 1) + leg_w)
  strip_chars <- vapply(lv, .comp_strip_chars, numeric(1), g_of = g_of, slot = slot)
  list(width = width, height = height, rotate = rotate, leg_cols = leg_cols,
       small_x = n > 20, strip_chars = strip_chars)
}

.comp_strip_chars <- function(gl, g_of, slot) max(8, floor(sum(g_of == gl) * slot / 0.08))

.comp_plot <- function(P, ids, g_of, lv, cats, tests, run_test, has_group, test_label, dims, spec, id, cat_col) {
  # The strip names its column, so a group level cannot read as a category.
  strip_text <- if (has_group) paste0(spec$data$group, ": ", lv) else lv
  strips <- vapply(seq_along(lv), .comp_strip_label, character(1), lv = strip_text,
                   chars = dims$strip_chars)
  long <- data.frame(
    sample = factor(rep(ids, times = length(cats)), levels = ids),
    panel = factor(strips[match(rep(g_of, times = length(cats)), lv)], levels = strips),
    category = factor(rep(cats, each = length(ids)), levels = cats),
    proportion = as.vector(P))
  pal <- pub_palette(length(cats), spec$appearance$palette)
  names(pal) <- cats
  labs_leg <- vapply(seq_along(cats), .comp_legend_label, character(1), cats = cats, tests = tests)
  n_sig <- sum(tests$p_adj < 0.05, na.rm = TRUE)
  sub <- if (run_test) {
    sprintf("%s per category, BH-adjusted\n%d of %d categories with P adj < 0.05%s",
            test_label, n_sig, length(cats), if (n_sig > 0) " (* in the legend)" else "")
  } else if (has_group) {
    sprintf("Composition of %d samples\nno test, because a group has fewer than 2 samples", length(ids))
  } else sprintf("Composition of %d samples (descriptive, no test)", length(ids))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = sample, y = proportion, fill = category)) +
    ggplot2::geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = pal, labels = labs_leg,
                               guide = ggplot2::guide_legend(ncol = dims$leg_cols)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), breaks = c(0, 0.25, 0.5, 0.75, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(x = spec$appearance$x_label %||% id, y = spec$appearance$y_label %||% "Proportion",
                  fill = cat_col, title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   legend.text = ggplot2::element_text(size = 9, lineheight = 0.95),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.key.size = grid::unit(0.9, "lines"),
                   legend.key.spacing.y = grid::unit(0.15, "lines"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 10, lineheight = 1.0,
                                                      margin = ggplot2::margin(2, 0, 4, 0)),
                   panel.spacing.x = grid::unit(0.5, "lines"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
  if (isTRUE(dims$rotate)) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = if (isTRUE(dims$small_x)) 8.5 else 10.4))
  }
  if (has_group) {
    p <- p + ggplot2::facet_grid(cols = ggplot2::vars(panel), scales = "free_x", space = "free_x")
  }
  p
}

.comp_strip_label <- function(i, lv, chars) .comp_wrap(lv[i], chars[i])

.comp_methods <- function(ids, cats, lv, n_by, has_group, run_test, test_label, tests, spec, id, cat_col) {
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  v <- pkg_versions("ggplot2")
  intro <- sprintf(paste0("The counts of %d categories of %s were summed within each of %d samples (%s)",
                          " and converted to proportions of the sample total."),
                   length(cats), cat_col, length(ids), id)
  grp_txt <- if (has_group) {
    sprintf(" The samples fall into %d groups (%s).", length(lv),
            paste(sprintf("%s: n = %d", lv, n_by), collapse = ", "))
  } else ""
  test_txt <- if (run_test) {
    n_sig <- sum(tests$p_adj < 0.05, na.rm = TRUE)
    sprintf(paste0(" For each category, the per-sample proportions were compared between the groups with",
                   " the %s%s, and the P values were adjusted across the %d categories with the",
                   " Benjamini-Hochberg method; %d categories had an adjusted P below 0.05."),
            test_label, if (length(lv) == 2) " (two-sided, with the exact null distribution when there are fewer than 50 samples per group and no ties)" else "",
            length(cats), n_sig)
  } else if (has_group) {
    " No test was run, because a group holds fewer than 2 samples."
  } else " The composition is described without a test."
  caveat <- paste0(" Proportions are compositional: they sum to one within each sample, so a rise in one",
                   " category forces the others down and the per-category tests are not independent.",
                   " The tests are therefore exploratory; a compositional method, for example a log-ratio",
                   " transform, is needed for a confirmatory claim.")
  paste0(intro, grp_txt, test_txt, if (has_group) caveat else "",
         sprintf(" Computed in R %s and drawn with ggplot2 %s.", rver, v$ggplot2))
}

.comp_emit_script <- function(spec, has_count, has_group, in_name, fig_stub) {
  id <- spec$data$id
  cat_col <- spec$data$x
  cnt_line <- if (has_count) sprintf('d$count <- as.numeric(raw[["%s"]])', spec$data$y) else "d$count <- 1"
  grp_line <- if (has_group) sprintf('d$group <- as.character(raw[["%s"]])', spec$data$group) else 'd$group <- "all"'
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2); library(scales) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(id = as.character(raw[["%s"]]), category = as.character(raw[["%s"]]))', id, cat_col),
    cnt_line,
    grp_line,
    "d <- d[stats::complete.cases(d), ]",
    "cats <- unique(d$category)",
    "tab <- xtabs(count ~ id + factor(category, levels = cats), data = d)",
    "P <- prop.table(tab, 1)",
    "P <- P[rowSums(tab) > 0, , drop = FALSE]",
    "g <- d$group[match(rownames(P), d$id)]",
    "lv <- sort(unique(g), method = \"radix\")",
    "if (length(lv) >= 2) for (j in seq_along(cats)) {",
    "  v <- as.numeric(P[, j])",
    "  if (length(unique(v)) < 2) next",
    "  if (length(lv) == 2) print(suppressWarnings(wilcox.test(v[g == lv[1]], v[g == lv[2]])))",
    "  else print(kruskal.test(v, factor(g, levels = lv)))",
    "}",
    "",
    "long <- data.frame(sample = rep(rownames(P), times = length(cats)), group = rep(g, times = length(cats)),",
    "                   category = factor(rep(cats, each = nrow(P)), levels = cats), proportion = as.vector(P))",
    "p <- ggplot(long, aes(x = sample, y = proportion, fill = category)) +",
    "  geom_col(width = 0.85, colour = \"white\", linewidth = 0.2) +",
    "  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.01))) +",
    "  facet_grid(cols = vars(group), scales = \"free_x\", space = \"free_x\") +",
    sprintf('  labs(x = "%s", y = "Proportion", fill = "%s") +', spec$appearance$x_label %||% id, cat_col),
    "  theme_classic(base_size = 13) +",
    "  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = element_blank())",
    "",
    sprintf('ggsave("%s.pdf", p, width = max(4.8, 2.5 + 0.3 * nrow(P)), height = 4.8)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
