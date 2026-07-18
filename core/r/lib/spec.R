# The figure spec: one structure that both `plot` (from CLI args) and
# `render` (from plot_config.yaml) build, then hand to a recipe.

default_appearance <- function() {
  list(
    theme       = "pubplot_house",
    palette     = NULL,            # NULL means use the house palette
    x_label     = NULL,
    y_label     = NULL,
    title       = NULL,
    y_limits    = list(NULL, NULL),
    show_points = TRUE,
    geom        = "box",           # box | violin | bar
    bracket     = list(show = TRUE, label = "p.signif"),
    scale       = "row",           # heatmap: row | column | none
    cluster     = "both"           # heatmap: both | rows | columns | none
  )
}

# Build a spec from parsed command-line options.
spec_from_opt <- function(opt) {
  ap <- default_appearance()
  if (!is.null(opt$geom))    ap$geom    <- opt$geom
  if (!is.null(opt$xlab))    ap$x_label <- opt$xlab
  if (!is.null(opt$ylab))    ap$y_label <- opt$ylab
  if (!is.null(opt$title))   ap$title   <- opt$title
  if (!is.null(opt$palette)) ap$palette <- trimws(strsplit(opt$palette, ",")[[1]])
  if (!is.null(opt$scale))   ap$scale   <- opt$scale
  if (!is.null(opt$cluster)) ap$cluster <- opt$cluster
  list(
    engine = "r",
    recipe = opt$recipe,
    data   = list(file = basename(opt$data), x = opt$x, y = opt$y,
                  fill = opt$fill, facet = opt$facet,
                  time = opt$time, event = opt$event,
                  annotation = opt$annotation,
                  covariates = if (!is.null(opt$covariates)) {
                    trimws(strsplit(opt$covariates, ",")[[1]])
                  } else NULL),
    test   = list(
      method   = opt$test %||% "auto",
      paired   = isTRUE(opt$paired),
      p_adjust = opt$p_adjust %||% "none"
    ),
    appearance = ap
  )
}

# YAML 1.1 reads a bare `y` key as the boolean TRUE, so a hand-edited
# `y: volume` arrives as a list element named "TRUE". Map the axis keys back so
# scientists can edit plot_config.yaml naturally without quoting `y`.
fix_axis_keys <- function(d) {
  if (is.null(d) || is.null(names(d))) return(d)
  nm <- names(d)
  if ("TRUE" %in% nm && !"y" %in% nm)  nm[nm == "TRUE"]  <- "y"
  if ("FALSE" %in% nm && !"n" %in% nm) nm[nm == "FALSE"] <- "n"
  names(d) <- nm
  d
}

# Build a spec from a plot_config.yaml file, filling any missing appearance keys.
spec_from_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  cfg$data <- fix_axis_keys(cfg$data)
  cfg$appearance <- modifyList(default_appearance(), cfg$appearance %||% list())
  cfg$test <- modifyList(list(method = "auto", paired = FALSE, p_adjust = "none"),
                         cfg$test %||% list())
  cfg
}

# Write a spec back out as plot_config.yaml.
write_config <- function(spec, path) {
  yaml::write_yaml(spec, path)
}
