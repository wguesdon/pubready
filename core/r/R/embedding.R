# Recipe: UMAP or t-SNE scatter coloured by cluster, from coordinates that the
# user computed upstream (Seurat, scanpy). pubready does not compute the
# embedding, because the two engines cannot give the same UMAP. Input is one row
# per cell: --x and --y are the two embedding coordinates and --group the
# cluster. Each cluster is labelled at its median position, and the stats give
# the size, the median position and the mean silhouette width of each cluster on
# the 2D coordinates. A numeric --group with more than 12 values (and not a set
# of at most 50 whole-number cluster ids) is a continuous score: it is coloured
# on the viridis scale and has no cluster labels.

.em_max_legend   <- 12
.em_max_int_ids  <- 50
.em_max_sil      <- 5000
.em_max_points   <- 50000

.em_is_continuous <- function(v) {
  if (!is.numeric(v)) return(FALSE)
  u <- unique(v[is.finite(v)])
  if (length(u) <= .em_max_legend) return(FALSE)
  !(all(u == round(u)) && length(u) <= .em_max_int_ids)
}

# Cluster levels: numbers in numeric order (also when written as text), other
# names in C-locale order.
.em_levels <- function(v) {
  if (is.numeric(v)) return(as.character(sort(unique(v))))
  u <- unique(as.character(v))
  num <- suppressWarnings(as.numeric(u))
  if (!anyNA(num)) return(u[order(num)])
  sort(u, method = "radix")
}

.em_fmt_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

# Mean silhouette width of each cluster and overall, from cluster::silhouette on
# the Euclidean distances of the 2D coordinates. Above .em_max_sil cells, every
# k-th cell in input order is used.
.em_silhouette <- function(xy, codes, k_levels) {
  n <- nrow(xy)
  step <- if (n > .em_max_sil) as.integer(ceiling(n / .em_max_sil)) else 1L
  idx <- seq(1L, n, by = step)
  out <- list(per_cluster = rep(NA_real_, k_levels), overall = NA_real_, n = length(idx), step = step)
  sil <- cluster::silhouette(codes[idx], stats::dist(xy[idx, , drop = FALSE]))
  if (!inherits(sil, "silhouette")) return(out)
  w <- sil[, "sil_width"]
  cl <- sil[, "cluster"]
  m <- tapply(w, factor(cl, levels = seq_len(k_levels)), mean)
  out$per_cluster <- as.numeric(m)
  out$overall <- mean(w)
  out
}

# The figure size, from the aspect ratio of the data, so that the equal-aspect
# panel fills the figure; the legend width grows with the longest cluster name.
.em_dims <- function(x, y, legend, continuous, lv) {
  r <- min(2, max(0.5, diff(range(y)) / max(diff(range(x)), 1e-12)))
  panel_w <- min(4.4, 5.6 / r)
  leg_w <- if (!legend) 0 else if (continuous) 1.2 else 0.8 + 0.085 * max(nchar(lv))
  c(panel_w + 0.7 + leg_w, panel_w * r + 1.1)
}

.em_size <- function(n) min(1.6, max(0.3, 1.6 * sqrt(500 / n)))

.em_alpha <- function(n) if (n <= 1000) 0.9 else max(0.5, 0.9 * (1000 / n)^0.2)

recipe_embedding <- function(df, spec) {
  xc <- spec$data$x; yc <- spec$data$y; gc <- spec$data$group
  if (is.null(xc) || !nzchar(xc) || is.null(yc) || !nzchar(yc)) {
    stop("embedding needs --x and --y (the two embedding coordinates, for example UMAP_1 and UMAP_2)")
  }
  if (is.null(gc) || !nzchar(gc)) stop("embedding needs --group (the cluster column)")
  for (cn in c(xc, yc, gc)) if (!cn %in% names(df)) stop(sprintf("column '%s' not found in data", cn))

  n0 <- nrow(df)
  x <- suppressWarnings(as.numeric(df[[xc]])); y <- suppressWarnings(as.numeric(df[[yc]]))
  g <- df[[gc]]
  keep <- is.finite(x) & is.finite(y) & !is.na(g)
  if (is.numeric(g)) keep <- keep & is.finite(g)
  x <- x[keep]; y <- y[keep]; g <- g[keep]
  n <- length(x)
  if (n < 2) stop("embedding needs at least 2 cells with both coordinates and a group")
  clean_steps <- c(
    if (sum(!keep) > 0) sprintf("Dropped %d row(s) with a missing coordinate or group (%d -> %d).",
                                sum(!keep), n0, n))
  if (length(clean_steps) == 0) clean_steps <- "No cleaning applied; input used as-is."

  continuous <- .em_is_continuous(g)
  step_draw <- if (n > .em_max_points) as.integer(ceiling(n / .em_max_points)) else 1L
  draw <- seq(1L, n, by = step_draw)
  pt <- .em_size(length(draw)); al <- .em_alpha(length(draw))
  x_lab <- spec$appearance$x_label %||% xc
  y_lab <- spec$appearance$y_label %||% yc

  if (continuous) {
    stats_df <- data.frame(cluster = "all", n = n, percent = 100, median_x = stats::median(x),
                           median_y = stats::median(y), mean_silhouette = NA_real_,
                           stringsAsFactors = FALSE)
    sil <- list(overall = NA_real_, n = 0L, step = 1L)
    lv <- character(0)
    pdat <- data.frame(x = x[draw], y = y[draw], v = as.numeric(g[draw]))
    sub <- sprintf("n = %s cells; colour = %s", .em_fmt_count(n), gc)
    p <- ggplot2::ggplot(pdat, ggplot2::aes(x = x, y = y, colour = v)) +
      ggplot2::geom_point(size = pt, alpha = al, stroke = 0, shape = 16) +
      ggplot2::scale_colour_viridis_c(name = gc) +
      ggplot2::guides(colour = ggplot2::guide_colourbar(barheight = grid::unit(7, "lines")))
    legend <- TRUE
  } else {
    lv <- .em_levels(g)
    gl <- as.character(g)
    codes <- match(gl, lv)
    counts <- tabulate(codes, nbins = length(lv))
    mx <- as.numeric(tapply(x, factor(codes, levels = seq_along(lv)), stats::median))
    my <- as.numeric(tapply(y, factor(codes, levels = seq_along(lv)), stats::median))
    sil <- if (length(lv) >= 2) .em_silhouette(cbind(x, y), codes, length(lv)) else
      list(per_cluster = NA_real_, overall = NA_real_, n = 0L, step = 1L)
    stats_df <- data.frame(cluster = lv, n = counts, percent = 100 * counts / n, median_x = mx,
                           median_y = my, mean_silhouette = sil$per_cluster, stringsAsFactors = FALSE)
    pal <- pub_palette(length(lv), spec$appearance$palette)
    pdat <- data.frame(x = x[draw], y = y[draw], cl = factor(gl[draw], levels = lv))
    cent <- data.frame(x = mx, y = my, lab = lv, stringsAsFactors = FALSE)
    legend <- length(lv) <= .em_max_legend
    sub <- sprintf("n = %s cells; %d clusters; mean silhouette = %s", .em_fmt_count(n), length(lv),
                   if (is.na(sil$overall)) "NA" else formatC(sil$overall, format = "f", digits = 2))
    p <- ggplot2::ggplot(pdat, ggplot2::aes(x = x, y = y, colour = cl)) +
      ggplot2::geom_point(size = pt, alpha = al, stroke = 0, shape = 16) +
      ggplot2::scale_colour_manual(values = pal, name = gc, drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1))) +
      ggrepel::geom_label_repel(data = cent, ggplot2::aes(x = x, y = y, label = lab),
                                inherit.aes = FALSE, size = if (length(lv) > 20) 2.6 else 3.1,
                                fill = grDevices::adjustcolor("white", alpha.f = 0.9),
                                label.size = 0, label.r = grid::unit(0.2, "lines"),
                                label.padding = grid::unit(0.15, "lines"), box.padding = 0.1,
                                point.size = NA, min.segment.length = Inf, force = 0.6,
                                max.overlaps = Inf, seed = 1)
  }

  p <- p +
    ggplot2::coord_equal() +
    ggplot2::labs(x = x_lab, y = y_lab, title = spec$appearance$title, subtitle = sub) +
    pub_theme(spec$appearance$theme) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   legend.position = if (legend) "right" else "none",
                   legend.text = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(size = 11),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, colour = "grey30"))

  dims <- .em_dims(x, y, legend, continuous, lv)
  test_meta <- list(
    name = "Embedding scatter", n_cells = n, n_clusters = length(lv), continuous = continuous,
    mean_silhouette = sil$overall, silhouette_n = sil$n, silhouette_step = sil$step,
    n_drawn = length(draw), drawing_step = step_draw)
  resolved <- list(family = "embedding", label = test_meta$name)
  df_used <- stats::setNames(data.frame(x, y, g, stringsAsFactors = FALSE), c(xc, yc, gc))
  list(plot = p, stats = stats_df, test_meta = test_meta, resolved = resolved,
       methods = .em_methods(test_meta, xc, yc, gc),
       build_script = list(fn = .em_emit_script,
                           args = list(xc = xc, yc = yc, gc = gc, continuous = continuous)),
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       width = dims[1], height = dims[2])
}

.em_methods <- function(tm, xc, yc, gc) {
  v <- pkg_versions(c("ggplot2", "ggrepel", "cluster"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  draw_txt <- if (tm$drawing_step > 1) {
    sprintf("For the size of the vector outputs, one in every %d cells in input order is drawn (%s of %s); the counts and the medians use all cells. ",
            tm$drawing_step, .em_fmt_count(tm$n_drawn), .em_fmt_count(tm$n_cells))
  } else ""
  if (isTRUE(tm$continuous)) {
    return(paste0(
      sprintf("The embedding coordinates %s and %s of %s cells were computed upstream and are drawn as given; pubready does not compute the embedding. ",
              xc, yc, .em_fmt_count(tm$n_cells)),
      sprintf("Points are coloured by the continuous value %s on the viridis scale; no cluster statistic applies. ", gc),
      draw_txt,
      sprintf("Drawn in R %s with ggplot2 %s.", rver, v$ggplot2)))
  }
  sub_txt <- if (tm$silhouette_step > 1) {
    sprintf("Because there were %s cells, the silhouette was computed on one in every %d cells in input order (%s cells). ",
            .em_fmt_count(tm$n_cells), tm$silhouette_step, .em_fmt_count(tm$silhouette_n))
  } else ""
  paste0(
    sprintf("The embedding coordinates %s and %s of %s cells were computed upstream and are drawn as given; pubready does not compute the embedding. ",
            xc, yc, .em_fmt_count(tm$n_cells)),
    sprintf("Points are coloured by %s (%d clusters), and each cluster is labelled at its median position. ", gc, tm$n_clusters),
    "The mean silhouette width of each cluster was computed from the Euclidean distances between cells on the two embedding coordinates (cluster::silhouette). ",
    "It describes how well the clusters separate in the embedding, not in the original feature space, because an embedding such as UMAP or t-SNE distorts distances. ",
    sub_txt,
    sprintf("The overall mean silhouette width was %s. ",
            if (is.na(tm$mean_silhouette)) "not defined" else formatC(tm$mean_silhouette, format = "f", digits = 3)),
    draw_txt,
    sprintf("Drawn in R %s with ggplot2 %s and ggrepel %s; silhouettes from cluster %s.",
            rver, v$ggplot2, v$ggrepel, v$cluster))
}

.em_emit_script <- function(xc, yc, gc, continuous, in_name, fig_stub) {
  colour_line <- if (continuous) "  scale_colour_viridis_c() +" else "  guides(colour = guide_legend(override.aes = list(size = 3))) +"
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({ library(ggplot2); library(cluster) })",
    "",
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name),
    sprintf('d <- data.frame(x = as.numeric(raw[["%s"]]), y = as.numeric(raw[["%s"]]), g = raw[["%s"]])', xc, yc, gc),
    "d <- d[is.finite(d$x) & is.finite(d$y) & !is.na(d$g), ]",
    if (!continuous) c(
      "d$g <- factor(d$g)",
      "step <- max(1, ceiling(nrow(d) / 5000)); i <- seq(1, nrow(d), by = step)",
      "s <- silhouette(as.integer(d$g[i]), dist(d[i, c(\"x\", \"y\")]))",
      "print(summary(s)$clus.avg.widths)",
      "cent <- aggregate(cbind(x, y) ~ g, data = d, FUN = median)") else character(0),
    "p <- ggplot(d, aes(x, y, colour = g)) + geom_point(size = 0.6, shape = 16) +",
    colour_line,
    if (!continuous) '  geom_label(data = cent, aes(label = g), colour = "black", fill = "white", label.size = 0) +' else character(0),
    sprintf('  coord_equal() + labs(x = "%s", y = "%s") + theme_classic(base_size = 13) +', xc, yc),
    "  theme(axis.text = element_blank(), axis.ticks = element_blank())",
    sprintf('ggsave("%s.pdf", p, width = 6.6, height = 5.2)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
}
