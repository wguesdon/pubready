# The pubplot house style. One clean publication theme for v1.

pubplot_palette <- c("#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A")

theme_pubplot <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line   = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      axis.ticks  = ggplot2::element_line(colour = "black"),
      axis.text   = ggplot2::element_text(colour = "black"),
      axis.title  = ggplot2::element_text(colour = "black"),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(10, 14, 10, 10),
      legend.position = "none"
    )
}

# Add the distribution geom for the comparison recipes, plus optional points.
# geom is one of box | violin | bar | raincloud. `x` is the grouping column name;
# fill is mapped per layer so the significance-bracket layer is not affected.
add_dist_geom <- function(p, x, geom = "box", show_points = TRUE) {
  fill <- ggplot2::aes(fill = .data[[x]])
  p <- if (geom == "violin") {
    p + ggplot2::geom_violin(fill, trim = FALSE, width = 0.7, alpha = 0.9)
  } else if (geom == "bar") {
    p +
      ggplot2::stat_summary(fill, fun = mean, geom = "bar", width = 0.6,
                            alpha = 0.9, colour = "black", linewidth = 0.3) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.2)
  } else if (geom == "raincloud") {
    p +
      ggdist::stat_halfeye(fill, adjust = 0.6, width = 0.6, .width = 0,
                           justification = -0.2, point_colour = NA, alpha = 0.85) +
      ggplot2::geom_boxplot(fill, width = 0.12, outlier.shape = NA, alpha = 0.9)
  } else {
    p + ggplot2::geom_boxplot(fill, width = 0.6, outlier.shape = NA, alpha = 0.9)
  }
  if (isTRUE(show_points)) {
    jit <- if (geom == "raincloud") 0.06 else 0.12
    p <- p + ggplot2::geom_jitter(fill, width = jit, size = 1.5, alpha = 0.7,
                                  shape = 21, stroke = 0.3)
  }
  p
}

# Resolve the ggplot theme from the spec. "prism" gives the GraphPad Prism look
# via ggprism::theme_prism; anything else is the house style. A custom --palette
# still applies on top, because the recipes set the fill scale themselves.
pub_theme <- function(theme = "pubplot_house", base_size = 13) {
  if (identical(theme, "prism")) {
    ggprism::theme_prism(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  } else {
    theme_pubplot(base_size)
  }
}
