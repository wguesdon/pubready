# The figure spec: one structure that both `plot` (from CLI args) and
# `render` (from plot_config.yaml) build, then hand to a recipe.

default_appearance <- function() {
  list(
    theme        = "pubready_house",
    palette      = NULL,           # NULL means use the house palette
    x_label      = NULL,
    y_label      = NULL,
    title        = NULL,
    y_limits     = list(NULL, NULL),
    show_points  = TRUE,
    geom         = "box",          # box | violin | bar | raincloud
    bracket      = list(show = TRUE, label = "p.signif"),
    scale        = "row",          # heatmap: row | column | none
    cluster      = "both",         # heatmap: both | rows | columns | none
    fc_cutoff    = 1.0,            # volcano: |log2FC| threshold
    p_cutoff     = 0.05,           # volcano: p threshold
    top_n        = 15,             # volcano / enrichment: labels / terms to show
    pd_threshold = 20,             # spider: RECIST progressive disease line
    pr_threshold = -30             # spider: RECIST partial response line
  )
}

# Split a list option such as --levels into its items. The items are split on
# "|" when the value holds one, so an item can hold a comma, for example
# --levels "edema, no diuretics|edema despite diuretics". Else they are split on
# ",". NULL stays NULL.
split_list_option <- function(x) {
  if (is.null(x)) return(NULL)
  sep <- if (grepl("|", x, fixed = TRUE)) "|" else ","
  trimws(strsplit(x, sep, fixed = TRUE)[[1]])
}

# Build a spec from parsed command-line options.
spec_from_opt <- function(opt) {
  ap <- default_appearance()
  if (!is.null(opt$geom))    ap$geom    <- opt$geom
  if (!is.null(opt$theme))   ap$theme   <- opt$theme
  if (!is.null(opt$xlab))    ap$x_label <- opt$xlab
  if (!is.null(opt$ylab))    ap$y_label <- opt$ylab
  if (!is.null(opt$title))   ap$title   <- opt$title
  if (!is.null(opt$palette)) ap$palette <- trimws(strsplit(opt$palette, ",")[[1]])
  if (!is.null(opt$scale))     ap$scale     <- opt$scale
  if (!is.null(opt$cluster))   ap$cluster   <- opt$cluster
  if (!is.null(opt$fc_cutoff)) ap$fc_cutoff <- opt$fc_cutoff
  if (!is.null(opt$p_cutoff))  ap$p_cutoff  <- opt$p_cutoff
  if (!is.null(opt$top_n))     ap$top_n     <- opt$top_n
  if (!is.null(opt$pd_threshold)) ap$pd_threshold <- opt$pd_threshold
  if (!is.null(opt$pr_threshold)) ap$pr_threshold <- opt$pr_threshold
  list(
    engine = "r",
    recipe = opt$recipe,
    data   = list(file = basename(opt$data), x = opt$x, y = opt$y,
                  fill = opt$fill, facet = opt$facet,
                  time = opt$time, event = opt$event,
                  annotation = opt$annotation,
                  id = opt$id, group = opt$group, label = opt$label,
                  positive = opt$positive, se = opt$se,
                  ci_low = opt$ci_low, ci_high = opt$ci_high,
                  measure = opt$measure,
                  ref = opt$ref, cause = opt$cause,
                  levels = split_list_option(opt$levels),
                  visit_levels = split_list_option(opt$visit_levels),
                  features = split_list_option(opt$features),
                  markers = split_list_option(opt$markers),
                  covariates = split_list_option(opt$covariates)),
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
