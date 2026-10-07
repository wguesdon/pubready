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

recipe_pca <- function(df, spec) {
  grp_col <- spec$data$group %||% spec$data$fill %||% spec$data$x
  if (is.null(grp_col) || !nzchar(grp_col)) stop("pca needs --group (the grouping column)")
  if (!grp_col %in% names(df)) stop(sprintf("column '%s' not found in data", grp_col))
  lab_col <- spec$data$label

  feat_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], c(grp_col, lab_col))
  if (length(feat_cols) < 2) stop("pca needs at least 2 numeric feature columns")
  mat <- as.matrix(df[, feat_cols, drop = FALSE])
  keepc <- apply(mat, 2, function(v) stats::sd(v, na.rm = TRUE) > 0)
  mat <- mat[, keepc, drop = FALSE]
  keepr <- stats::complete.cases(mat)
  mat <- mat[keepr, , drop = FALSE]
  grp_chr <- as.character(df[[grp_col]])[keepr]
  grp <- factor(grp_chr, levels = sort(unique(grp_chr), method = "radix"))
  if (nlevels(grp) < 2) stop("pca needs at least 2 groups")
  labs <- if (!is.null(lab_col) && nzchar(lab_col)) as.character(df[[lab_col]])[keepr] else NULL
  clean_steps <- c(
    sprintf("Used %d samples and %d numeric features.", nrow(mat), ncol(mat)),
    if (sum(!keepc) > 0) sprintf("Dropped %d zero-variance feature(s).", sum(!keepc)) else NULL,
    if (sum(!keepr) > 0) sprintf("Dropped %d sample(s) with missing values.", sum(!keepr)) else NULL)
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
                         variance_explained_pct = round(ve[seq_len(k)], 3))
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
    'num <- raw[, sapply(raw, is.numeric), drop = FALSE]',
    sprintf('grp <- factor(raw[["%s"]])', grp_col),
    "num <- num[, apply(num, 2, function(v) sd(v) > 0), drop = FALSE]",
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
