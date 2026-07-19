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
