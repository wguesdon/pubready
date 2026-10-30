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

# The options of figkit plot and figkit diagnose, one list for both commands, so
# diagnose accepts every option that plot accepts. An option that the user does
# not give stays NULL, and spec_from_opt() fills its default. The Python engine
# (spec.add_plot_arguments) defines the same options.
plot_option_list <- function() {
  chr <- c("recipe", "data", "x", "y", "fill", "facet", "time", "event", "covariates",
           "annotation", "id", "group", "label", "scale", "cluster")
  opts <- lapply(paste0("--", chr), optparse::make_option, type = "character", default = NULL)
  c(opts, list(
    optparse::make_option("--fc_cutoff", type = "double",   default = NULL),
    optparse::make_option("--p_cutoff",  type = "double",   default = NULL),
    optparse::make_option("--top_n",     type = "integer",  default = NULL),
    optparse::make_option("--pd_threshold", type = "double", default = NULL),
    optparse::make_option("--pr_threshold", type = "double", default = NULL),
    optparse::make_option("--positive", type = "character", default = NULL),
    optparse::make_option("--se",       type = "character", default = NULL),
    optparse::make_option("--ci_low",   type = "character", default = NULL),
    optparse::make_option("--ci_high",  type = "character", default = NULL),
    optparse::make_option("--measure",  type = "character", default = NULL),
    optparse::make_option("--markers",  type = "character", default = NULL),
    optparse::make_option("--ref",      type = "character", default = NULL),
    optparse::make_option("--cause",    type = "character", default = NULL),
    optparse::make_option("--levels",   type = "character", default = NULL),
    optparse::make_option("--visit_levels", type = "character", default = NULL),
    optparse::make_option("--features", type = "character", default = NULL),
    optparse::make_option("--test",     type = "character", default = NULL),
    optparse::make_option("--paired",   action = "store_true", default = FALSE),
    optparse::make_option("--geom",     type = "character", default = NULL),
    optparse::make_option("--theme",    type = "character", default = NULL),
    optparse::make_option("--xlab",     type = "character", default = NULL),
    optparse::make_option("--ylab",     type = "character", default = NULL),
    optparse::make_option("--title",    type = "character", default = NULL),
    optparse::make_option("--palette",  type = "character", default = NULL),
    optparse::make_option("--sheet",    type = "character", default = NULL),
    optparse::make_option("--stamp",    type = "character", default = NULL),
    # Read by the figkit dispatcher for routing; ignored here.
    optparse::make_option("--engine",   type = "character", default = NULL),
    optparse::make_option("--out",      type = "character", default = "pubready_output")
  ))
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
      paired   = isTRUE(opt$paired)
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
  cfg$test <- modifyList(list(method = "auto", paired = FALSE), cfg$test %||% list())
  # --p_adjust was removed: no recipe applied it. An old config drops the key.
  cfg$test$p_adjust <- NULL
  cfg
}

# The path of a file that a plot_config.yaml names. A relative path is found in
# the folder of the config first, so a bundle renders from any folder, and else
# in the current folder, as before. NULL stays NULL.
resolve_config_path <- function(path, config) {
  if (is.null(path) || !nzchar(path) || grepl("^(/|~|[A-Za-z]:)", path)) return(path)
  near <- file.path(dirname(config), path)
  if (file.exists(near)) near else path
}

# Write a spec back out as plot_config.yaml.
write_config <- function(spec, path) {
  yaml::write_yaml(spec, path)
}
