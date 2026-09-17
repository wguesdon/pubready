# Recipe: UpSet plot of set intersections from a binary membership matrix.
# Input is a CSV where each row is an element and each column is a set; an
# optional leading id column names the elements. Cells are 1/0 or TRUE/FALSE.
# The plot is descriptive (intersection sizes), so there is no hypothesis test;
# the stats table lists set sizes and every intersection size. R draws it with
# ComplexHeatmap's UpSet (make_comb_mat), which needs no extra package.

recipe_upset <- function(df, spec) {
  first_is_id <- !is.numeric(df[[1]]) && !is.logical(df[[1]])
  if (first_is_id) {
    ids   <- as.character(df[[1]])
    setdf <- df[, -1, drop = FALSE]
  } else {
    ids   <- as.character(seq_len(nrow(df)))
    setdf <- df
  }
  if (ncol(setdf) < 2) stop("upset needs at least 2 set columns")

  as_member <- function(v) {
    if (is.logical(v)) return(v)
    vn <- suppressWarnings(as.numeric(v))
    if (all(is.na(vn) | vn %in% c(0, 1))) return(!is.na(vn) & vn == 1)
    tolower(as.character(v)) %in% c("1", "true", "yes", "y", "t")
  }
  bin <- as.data.frame(lapply(setdf, as_member))
  names(bin)    <- names(setdf)
  rownames(bin) <- ids

  n_elements <- nrow(bin)
  in_any  <- rowSums(as.matrix(bin)) > 0
  dropped <- sum(!in_any)
  bin <- bin[in_any, , drop = FALSE]
  clean_steps <- if (dropped > 0) {
    sprintf("Dropped %d element(s) that belonged to no set (%d -> %d).",
            dropped, n_elements, nrow(bin))
  } else "No cleaning applied; every element belongs to at least one set."
  df_used <- data.frame(element = rownames(bin), bin, check.names = FALSE)

  m <- ComplexHeatmap::make_comb_mat(bin)
  draw_fn <- function() ComplexHeatmap::draw(ComplexHeatmap::UpSet(m))

  ss <- ComplexHeatmap::set_size(m)
  cs <- ComplexHeatmap::comb_size(m)
  cd <- ComplexHeatmap::comb_degree(m)
  cn <- ComplexHeatmap::comb_name(m)
  set_names <- names(ss)
  combo_sets <- vapply(cn, function(code) {
    bits <- strsplit(code, "")[[1]] == "1"
    paste(set_names[bits], collapse = " & ")
  }, character(1))
  ord <- order(cs, decreasing = TRUE)

  stats_df <- rbind(
    data.frame(type = "set", sets = set_names, degree = NA_integer_,
               size = as.integer(ss), stringsAsFactors = FALSE),
    data.frame(type = "intersection", sets = combo_sets[ord],
               degree = as.integer(cd[ord]), size = as.integer(cs[ord]),
               stringsAsFactors = FALSE))

  test_meta <- list(
    name = "UpSet plot (ComplexHeatmap)",
    n_sets = length(ss), n_elements = nrow(bin),
    n_intersections = length(cs), largest_intersection = max(cs))

  resolved <- list(family = "upset", label = "UpSet plot")
  methods  <- .up_methods(length(ss), nrow(bin), length(cs), set_names)
  build_script <- function(in_name, fig_stub) .up_emit_script(in_name, fig_stub, first_is_id)

  list(plot = NULL, draw = draw_fn, stats = stats_df, test_meta = test_meta,
       resolved = resolved, methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       width  = max(6.0, 0.45 * length(cs) + 2.5),
       height = max(3.8, 0.45 * length(ss) + 2.8))
}

.up_methods <- function(n_sets, n_elem, n_inter, set_names) {
  v    <- pkg_versions("ComplexHeatmap")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  paste0(
    sprintf("Set intersections among %d sets (%s) across %d elements were drawn as an UpSet plot with ComplexHeatmap. ",
            n_sets, paste(set_names, collapse = ", "), n_elem),
    sprintf("The plot shows the %d observed intersections and each set's total size; ", n_inter),
    "it is descriptive, so no statistical test is applied.",
    sprintf(" Rendered in R %s with ComplexHeatmap %s.", rver, v$ComplexHeatmap))
}

.up_emit_script <- function(in_name, fig_stub, first_is_id) {
  read_line <- if (first_is_id) {
    sprintf('raw <- read.csv("%s", check.names = FALSE, row.names = 1)', in_name)
  } else {
    sprintf('raw <- read.csv("%s", check.names = FALSE)', in_name)
  }
  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction. Run inside the pinned pubready container",
    "# (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages(library(ComplexHeatmap))",
    "",
    read_line,
    "bin <- as.data.frame(lapply(raw, function(v) as.numeric(v) == 1))",
    "bin <- bin[rowSums(as.matrix(bin)) > 0, , drop = FALSE]",
    "m <- make_comb_mat(bin)",
    "",
    sprintf('cairo_pdf("%s.pdf", width = 7, height = 4.5)', fig_stub),
    "draw(UpSet(m)); dev.off()",
    'cat("Reproduced figure.\\n")'
  )
}
