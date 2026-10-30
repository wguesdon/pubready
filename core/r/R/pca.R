# Recipe: principal component analysis scatter. Input has rows = samples, numeric
# feature columns, and a grouping column (--group); an optional --label names the
# points. Features are centered and scaled, PC1 vs PC2 is drawn colored by group
# with 95% confidence ellipses, and a PERMANOVA p-value for group separation is
# put on the plot.

# One-way PERMANOVA (Anderson 2001) on Euclidean distances of the scaled data.
# The total sum of squares is fixed under permutation, so only the within-group
# sum of squares changes. A design with 5000 distinct labelings or fewer, such as
# two groups of four, is tested exactly by enumerating every labeling, so both
# engines report the same P. A larger design takes 999 seeded random
# permutations, and its P differs between the engines in the third digit.
.pm_max_exact <- 5000

.pm_ss_group <- function(lv, mat, g) {
  xi <- mat[g == lv, , drop = FALSE]
  if (nrow(xi) < 1) return(0)
  sum(sweep(xi, 2, colMeans(xi))^2)
}

.pm_fstat <- function(g, mat, lv, ss_t) {
  ssw <- sum(vapply(lv, .pm_ss_group, numeric(1), mat = mat, g = g))
  ((ss_t - ssw) / (length(lv) - 1)) / (ssw / (nrow(mat) - length(lv)))
}

# Every distinct way to split the row indices into groups of the given sizes.
.pm_assignments <- function(idx, sizes) {
  if (length(sizes) == 1) return(list(list(idx)))
  out <- list()
  for (first in utils::combn(idx, sizes[1], simplify = FALSE)) {
    for (rest in .pm_assignments(setdiff(idx, first), sizes[-1])) {
      out[[length(out) + 1]] <- c(list(first), rest)
    }
  }
  out
}

.pm_labels <- function(assignment, lv, n) {
  g <- character(n)
  for (k in seq_along(lv)) g[assignment[[k]]] <- lv[k]
  g
}

.pm_shuffle <- function(i, g) sample(g)

.permanova <- function(X, grp, nperm = 999, seed = 1) {
  X <- as.matrix(X); g <- as.character(grp); lv <- levels(grp)
  sizes <- as.integer(table(factor(g, levels = lv)))
  ss_t <- sum(sweep(X, 2, colMeans(X))^2)
  obs <- .pm_fstat(g, X, lv, ss_t)
  n_labelings <- exp(lfactorial(length(g)) - sum(lfactorial(sizes)))
  if (n_labelings <= .pm_max_exact) {
    labs <- lapply(.pm_assignments(seq_along(g), sizes), .pm_labels, lv = lv, n = length(g))
    fs <- vapply(labs, .pm_fstat, numeric(1), mat = X, lv = lv, ss_t = ss_t)
    return(list(F = obs, p = mean(fs >= obs - 1e-9 * abs(obs)), method = "exact",
                permutations = length(fs)))
  }
  set.seed(seed)
  fs <- vapply(lapply(seq_len(nperm), .pm_shuffle, g = g), .pm_fstat, numeric(1),
               mat = X, lv = lv, ss_t = ss_t)
  list(F = obs, p = (sum(fs >= obs) + 1) / (nperm + 1), method = "monte carlo",
       permutations = nperm)
}

# Column names in quotes, the first 20 and a count of the rest. The Python
# engine writes the same list.
.pca_name_list <- function(cols, max_n = 20) {
  q <- sprintf("'%s'", utils::head(cols, max_n))
  if (length(cols) > max_n) q <- c(q, sprintf("and %d more", length(cols) - max_n))
  paste(q, collapse = ", ")
}

# Data row numbers, counted from 1 below the header: the first 10 and a count
# of the rest. The Python engine writes the same list.
.pca_row_list <- function(rows, max_n = 10) {
  txt <- paste(utils::head(rows, max_n), collapse = ", ")
  if (length(rows) > max_n) txt <- sprintf("%s and %d more", txt, length(rows) - max_n)
  paste(if (length(rows) == 1) "row" else "rows", txt)
}

# The feature columns of the PCA. With --features, the listed columns, each of
# which must exist and be numeric. Else every numeric column other than the
# group and the --id and --label columns. coerce_numeric_columns() decides what
# is numeric: a column with one text cell such as "<0.1" stays a feature, and
# the cell becomes NA (F42). Returns the table, the features and the steps that
# name the coerced cells and the columns left out.
.pca_feature_columns <- function(df, grp_col, skip, features) {
  if (length(features) > 0) {
    absent <- setdiff(features, names(df))
    if (length(absent) > 0) {
      stop(sprintf("pca: --features names column(s) not found in data: %s",
                   paste(sprintf("'%s'", absent), collapse = ", ")))
    }
    co <- coerce_numeric_columns(df, features)
    if (length(co$text) > 0) {
      stop(sprintf("pca: --features names column(s) that are not numeric: %s",
                   paste(sprintf("'%s'", co$text), collapse = ", ")))
    }
    return(list(df = co$df, features = features, steps = co$steps))
  }
  skip <- intersect(skip[nzchar(skip)], names(df))
  cand <- setdiff(names(df), c(grp_col, skip))
  co <- coerce_numeric_columns(df, cand)
  steps <- co$steps
  if (length(skip) > 0) {
    steps <- c(steps, sprintf("Left out the --id and --label column(s) %s.",
                              paste(sprintf("'%s'", skip), collapse = ", ")))
  }
  if (length(co$text) > 0) {
    steps <- c(steps, sprintf("Left out the text column(s) %s.", .pca_name_list(co$text)))
  }
  list(df = co$df, features = co$numeric, steps = steps)
}

recipe_pca <- function(df, spec) {
  grp_col <- spec$data$group %||% spec$data$fill %||% spec$data$x
  if (is.null(grp_col) || !nzchar(grp_col)) stop("pca needs --group (the grouping column)")
  if (!grp_col %in% names(df)) stop(sprintf("column '%s' not found in data", grp_col))
  lab_col <- spec$data$label
  id_col <- spec$data$id

  # The features: the --features list, else every numeric column other than the
  # group, --id and --label. The data log names each feature, so a numeric
  # sample id that slipped in is visible (F47).
  sel <- .pca_feature_columns(df, grp_col, c(id_col, lab_col), spec$data$features)
  df <- sel$df
  feat_cols <- sel$features
  if (length(feat_cols) < 2) stop("pca needs at least 2 numeric feature columns")
  mat <- as.matrix(df[, feat_cols, drop = FALSE])
  # One rule in both engines: drop each sample with a missing group, then each
  # sample with a missing feature value, then each feature with no variance in
  # the samples that remain. The data log names each drop (F48).
  grp_all <- level_labels(df[[grp_col]])
  no_grp <- is.na(grp_all) | !nzchar(trimws(grp_all))
  no_feat <- !no_grp & !stats::complete.cases(mat)
  keepr <- !no_grp & !no_feat
  mat <- mat[keepr, , drop = FALSE]
  keepc <- (apply(mat, 2, stats::sd) > 0) %in% TRUE
  mat <- mat[, keepc, drop = FALSE]
  grp_chr <- grp_all[keepr]
  grp <- factor(grp_chr, levels = level_order(grp_chr))
  if (nlevels(grp) < 2) stop("pca needs at least 2 groups")
  labs <- if (!is.null(lab_col) && nzchar(lab_col)) as.character(df[[lab_col]])[keepr] else NULL
  clean_steps <- c(
    sprintf("Used %d samples and %d numeric features: %s.", nrow(mat), ncol(mat),
            .pca_name_list(colnames(mat))),
    sel$steps,
    if (any(no_grp)) sprintf("Dropped %d sample(s) with a missing group in '%s': %s.",
                             sum(no_grp), grp_col, .pca_row_list(which(no_grp))) else NULL,
    if (any(no_feat)) sprintf("Dropped %d sample(s) with a missing feature value: %s.",
                              sum(no_feat), .pca_row_list(which(no_feat))) else NULL,
    if (sum(!keepc) > 0) sprintf("Dropped %d zero-variance feature(s): %s.", sum(!keepc),
                                 .pca_name_list(feat_cols[!keepc])) else NULL)
  feat_cols <- colnames(mat)
  df_used <- df[keepr, , drop = FALSE]

  pc <- stats::prcomp(mat, center = TRUE, scale. = TRUE)
  ve <- (pc$sdev^2) / sum(pc$sdev^2) * 100
  scores <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = grp)
  if (!is.null(labs)) scores$lab <- labs

  perm <- .permanova(scale(mat), grp, nperm = 999, seed = 1)

  pal <- spec$appearance$palette %||% pubready_palette
  pal <- rep(pal, length.out = nlevels(grp))
  psub <- sprintf("PERMANOVA: F = %.2f, p = %s", perm$F, formatC(perm$p, format = "g", digits = 2))

  p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = group)) +
    ggplot2::geom_point(size = 2.4, alpha = 0.9) +
    ggplot2::stat_ellipse(ggplot2::aes(group = group), type = "norm", level = 0.95, linewidth = 0.5) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
                  colour = grp_col, title = spec$appearance$title, subtitle = psub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(legend.position = "right",
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))
  if (!is.null(labs)) {
    p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = lab), size = 2.6,
                                      colour = "grey30", max.overlaps = 15, seed = 1)
  }

  k <- min(5, length(ve))
  stats_df <- data.frame(PC = paste0("PC", seq_len(k)),
                         variance_explained_pct = round(ve[seq_len(k)], 3),
                         omnibus_test = "PERMANOVA", omnibus_statistic = perm$F,
                         omnibus_df1 = nlevels(grp) - 1, omnibus_df2 = nrow(mat) - nlevels(grp),
                         omnibus_p = perm$p, omnibus_method = perm$method,
                         omnibus_permutations = perm$permutations, stringsAsFactors = FALSE)
  test_meta <- list(name = "PCA with PERMANOVA",
                    pc1_variance_pct = ve[1], pc2_variance_pct = ve[2],
                    n_samples = nrow(mat), n_features = ncol(mat),
                    permanova = list(statistic_F = perm$F, p_value = perm$p, method = perm$method,
                                     permutations = perm$permutations))
  resolved <- list(family = "pca", label = "PCA")
  methods  <- .pca_methods(nrow(mat), ncol(mat), ve, perm, nlevels(grp))
  build_script <- function(in_name, fig_stub) .pca_emit_script(grp_col, feat_cols, in_name, fig_stub)

  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps, width = 5.2, height = 4.6)
}

.pm_method_text <- function(perm) {
  if (identical(perm$method, "exact")) {
    sprintf("exact test over all %d distinct labelings", as.integer(perm$permutations))
  } else sprintf("%d seeded random permutations", as.integer(perm$permutations))
}

.pca_methods <- function(n_samp, n_feat, ve, perm, n_grp) {
  v    <- pkg_versions("ggplot2")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  paste0(
    sprintf("Principal component analysis was run on %d samples by %d centered and scaled features. ",
            n_samp, n_feat),
    sprintf("PC1 and PC2 explained %.1f%% and %.1f%% of the variance and are shown with 95%% normal-theory confidence ellipses per group (drawn for groups of four samples or more). ",
            ve[1], ve[2]),
    sprintf("Group separation was tested with a one-way PERMANOVA on Euclidean distances (%s): F = %.2f, p = %s. ",
            .pm_method_text(perm), perm$F, formatC(perm$p, format = "g", digits = 2)),
    sprintf("Drawn in R %s with ggplot2 %s.", rver, v$ggplot2))
}

.pca_emit_script <- function(grp_col, feat_cols, in_name, fig_stub) {
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container.",
    "suppressPackageStartupMessages(library(ggplot2))",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    "# The features of the figure; a cell that does not read as a number is NA,",
    "# and a sample with a missing feature or group is dropped, as in the figure.",
    sprintf("feats <- %s", paste(deparse(feat_cols, width.cutoff = 500L), collapse = "")),
    "num <- raw[feats]",
    "for (f in feats) num[[f]] <- suppressWarnings(as.numeric(as.character(num[[f]])))",
    sprintf('g <- as.character(raw[["%s"]])', grp_col),
    'keep <- complete.cases(num) & !is.na(g) & trimws(g) != ""',
    "num <- num[keep, , drop = FALSE]; grp <- factor(g[keep])",
    "pc <- prcomp(num, center = TRUE, scale. = TRUE)",
    "ve <- pc$sdev^2 / sum(pc$sdev^2) * 100",
    "d <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], group = grp)",
    "p <- ggplot(d, aes(PC1, PC2, colour = group)) + geom_point(size = 2.4) +",
    "  stat_ellipse(type = \"norm\", level = 0.95) +",
    '  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic()',
    sprintf('ggsave("%s.pdf", p, width = 5.2, height = 4.6)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
