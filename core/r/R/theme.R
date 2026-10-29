# The pubready house style. One clean publication theme for v1.

# The first six colours are the house palette. The fourteen after them keep a
# figure with up to twenty groups distinct; pub_palette() cycles past twenty.
pubready_palette <- c("#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A",
                      "#56B4E9", "#E78AC3", "#8C6D31", "#66C2A5", "#F28E2B", "#A6D854",
                      "#1F3A5F", "#B07AA1", "#9C755F", "#BAB0AC", "#17BECF", "#FF9DA7",
                      "#4E7A27", "#D4A6C8")

# n colours from a custom palette or the house palette, cycled when n is larger
# than the palette. Mirrors palette() of the Python engine.
pub_palette <- function(n, custom = NULL) {
  cols <- if (length(custom)) custom else pubready_palette
  rep_len(cols, n)
}

theme_pubready <- function(base_size = 13) {
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

# Wrap one label at the spaces so that no line is longer than `width`
# characters. A word longer than `width` stays whole on its own line. The Python
# engine has the same greedy rule (theme.wrap_label), so both engines break a
# label at the same places.
wrap_label <- function(s, width = 12) {
  s <- as.character(s)
  if (is.na(s) || nchar(s) <= width) return(s)
  words <- strsplit(s, " ", fixed = TRUE)[[1]]
  words <- words[nzchar(words)]
  if (length(words) == 0) return(s)
  lines <- character(0)
  cur <- words[1]
  for (w in words[-1]) {
    if (nchar(cur) + 1 + nchar(w) <= width) {
      cur <- paste(cur, w)
    } else {
      lines <- c(lines, cur)
      cur <- w
    }
  }
  paste(c(lines, cur), collapse = "\n")
}

# Vectorised wrap_label, for the category labels of an axis or a legend.
wrap_labels <- function(x, width = 12) {
  vapply(as.character(x), wrap_label, character(1), width = width, USE.NAMES = FALSE)
}

# The length in characters of the longest line of a wrapped label.
longest_line <- function(s) {
  max(nchar(strsplit(as.character(s), "\n", fixed = TRUE)[[1]]), 0L)
}

# The figure width of a category axis: 1.1 in per group plus 1 in for the axis,
# and never less than `min_width`. A wrapped label whose longest line (a long
# word) is wider than the slot widens every slot to 0.085 in per character.
# Python uses the same formula.
category_width <- function(n, labels = NULL, min_width = 4.2) {
  chars <- if (length(labels)) max(vapply(labels, longest_line, integer(1))) else 0L
  max(min_width, max(1.1, 0.085 * chars) * n + 1.0)
}

# The width of a legend right of the panel: the key plus 0.075 in per character
# of its longest line or of its title. Python uses the same formula.
legend_width <- function(labels, title = "") {
  chars <- max(vapply(as.character(labels), longest_line, integer(1)), nchar(title), 0L)
  0.55 + 0.075 * chars
}

# Format axis breaks as plain numbers with a comma between the thousands, never
# in scientific notation. The count of decimals is the smallest one that shows
# every break exactly (at most 6). Python has the same rule in
# theme.PlainNumberFormatter.
plain_number_labels <- function(x) {
  fin <- x[is.finite(x)]
  digits <- 0L
  while (digits < 6L && length(fin) > 0 &&
         any(abs(round(fin, digits) - fin) > 1e-9 * pmax(1, abs(fin)))) {
    digits <- digits + 1L
  }
  out <- formatC(x + 0, format = "f", digits = digits, big.mark = ",")
  out[!is.finite(x)] <- NA_character_
  out
}

# Log-axis ticks of the 1-2-5 sequence for a plot with a reference line `ref`
# (the null value of a ratio, 1). The ticks inside the data range are kept. When
# the data reach past `ref` on one side but no tick lies on that side, the next
# value of the sequence on that side is added, so the axis reads on both sides
# of the line. More than 7 ticks are thinned to the powers of ten. Python has
# the same rule in theme.log_ticks_ref().
log_ticks_ref <- function(lo, hi, ref = 1) {
  seq125 <- sort(as.vector(outer(c(1, 2, 5), 10^(-4:4))))
  tol <- 1e-9
  lo2 <- min(lo, ref); hi2 <- max(hi, ref)
  ticks <- seq125[seq125 >= lo2 * (1 - tol) & seq125 <= hi2 * (1 + tol)]
  if (hi > ref && !any(ticks > ref * (1 + tol))) {
    ticks <- c(ticks, min(seq125[seq125 > ref * (1 + tol)]))
  }
  if (lo < ref && !any(ticks < ref * (1 - tol))) {
    ticks <- c(max(seq125[seq125 < ref * (1 - tol)]), ticks)
  }
  if (length(ticks) > 7) {
    ticks <- ticks[abs(log10(ticks) - round(log10(ticks))) < 1e-9]
  }
  sort(unique(ticks))
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
pub_theme <- function(theme = "pubready_house", base_size = 13) {
  if (identical(theme, "prism")) {
    ggprism::theme_prism(base_size = base_size) +
      ggplot2::theme(legend.position = "none")
  } else {
    theme_pubready(base_size)
  }
}
