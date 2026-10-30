# Recipe: UpSet plot of set intersections from a binary membership matrix.
# Input is a CSV where each row is an element and each column is a set; an
# optional leading id column names the elements. Cells are 1/0 or TRUE/FALSE.
# The plot is descriptive (intersection sizes), so there is no hypothesis test;
# the stats table lists set sizes and every intersection size. R draws it with
# ComplexHeatmap's UpSet (make_comb_mat), which needs no extra package.
#
# The figure shows the top_n largest intersections (15 by default, --top_n sets
# it) and says in its title how many it hides. The intersections are ordered by
# size, then by degree, then by their set names, and the sets by size; the
# Python engine uses the same order, so both engines show the same bars.

# One membership rule in both engines (F49). A cell is a member for 1, TRUE, T,
# yes or y, and not a member for 0, FALSE, F, no or n, in any case. An empty cell
# is not a member, and the data log counts it. Any other value, such as a count
# of 2 or an "x", stops the run and names the column. The Python engine has the
# same rule and message.
.up_yes <- c("1", "true", "t", "yes", "y")
.up_no <- c("0", "false", "f", "no", "n")

# The cells of one column as text: a whole number with no decimal point, as
# level_labels() writes it, and NA for an empty cell.
.up_cells <- function(v) {
  s <- if (is.logical(v)) ifelse(v, "TRUE", "FALSE")
       else if (is.numeric(v)) level_labels(v)
       else trimws(as.character(v))
  s[!is.na(s) & !nzchar(s)] <- NA_character_
  s
}

.up_as_member <- function(v, col) {
  s <- .up_cells(v)
  code <- tolower(s)
  bad <- unique(s[!is.na(code) & !code %in% c(.up_yes, .up_no)])
  if (length(bad) > 0) {
    stop(sprintf("upset: column '%s' holds the value(s) %s, which are not a membership code; use 1/0, TRUE/FALSE or yes/no, and leave a cell empty for not a member",
                 col, paste(sprintf("'%s'", utils::head(bad, 5)), collapse = ", ")))
  }
  list(member = !is.na(code) & code %in% .up_yes, n_empty = sum(is.na(code)))
}

.up_combo_sets <- function(code, set_names) {
  bits <- strsplit(code, "")[[1]] == "1"
  paste(set_names[bits], collapse = " & ")
}

# The subtitle that says how many intersections the figure hides.
.up_subtitle <- function(n_shown, n_total) {
  if (n_shown >= n_total) return(NULL)
  sprintf("The %d largest of %d intersections; %d smaller ones are not shown",
          n_shown, n_total, n_total - n_shown)
}

# The width of the set-name column in inches, from the longest name.
.up_label_inches <- function(set_names) 0.3 + 0.08 * max(nchar(set_names), 4)

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

  bin <- data.frame(row.names = seq_len(nrow(setdf)))
  empty_steps <- character(0)
  for (cn in names(setdf)) {
    mem <- .up_as_member(setdf[[cn]], cn)
    bin[[cn]] <- mem$member
    if (mem$n_empty > 0) {
      empty_steps <- c(empty_steps, sprintf("Read %d empty cell(s) of column '%s' as not a member.",
                                            mem$n_empty, cn))
    }
  }
  rownames(bin) <- ids

  n_elements <- nrow(bin)
  in_any  <- rowSums(as.matrix(bin)) > 0
  dropped <- sum(!in_any)
  bin <- bin[in_any, , drop = FALSE]
  clean_steps <- c(
    empty_steps,
    if (dropped > 0) sprintf("Dropped %d element(s) that belonged to no set (%d -> %d).",
                             dropped, n_elements, nrow(bin)))
  if (length(clean_steps) == 0) {
    clean_steps <- "No cleaning applied; every element belongs to at least one set."
  }
  df_used <- data.frame(element = rownames(bin), bin, check.names = FALSE)

  m <- ComplexHeatmap::make_comb_mat(bin)
  ss <- ComplexHeatmap::set_size(m)
  cs <- ComplexHeatmap::comb_size(m)
  cd <- ComplexHeatmap::comb_degree(m)
  cn <- ComplexHeatmap::comb_name(m)
  set_names <- names(ss)
  combo_sets <- vapply(cn, .up_combo_sets, character(1), set_names = set_names)
  # Size first, then degree, then the set names in the C locale, so that ties
  # sort the same way in both engines.
  ord <- order(-cs, cd, combo_sets, method = "radix")
  top_n <- max(1L, as.integer(spec$appearance$top_n %||% 15))
  shown <- ord[seq_len(min(length(ord), top_n))]

  # Subsetting a combination matrix keeps the full set sizes.
  m_shown <- m[shown]
  comb_order <- match(cn[shown], ComplexHeatmap::comb_name(m_shown))
  set_order <- match(set_names[order(-ss, seq_along(ss), method = "radix")],
                     ComplexHeatmap::set_name(m_shown))
  subtitle <- .up_subtitle(length(shown), length(cs))
  title <- spec$appearance$title
  if (!is.null(title) && !is.null(subtitle)) title <- paste(title, subtitle, sep = "\n")
  if (is.null(title)) title <- subtitle

  up <- ComplexHeatmap::UpSet(
    m_shown, comb_order = comb_order, set_order = set_order,
    row_names_gp = grid::gpar(fontsize = 10),
    top_annotation = ComplexHeatmap::upset_top_annotation(
      m_shown, add_numbers = TRUE, numbers_rot = 0,
      numbers_gp = grid::gpar(fontsize = 8),
      axis_param = list(gp = grid::gpar(fontsize = 8)),
      annotation_name_gp = grid::gpar(fontsize = 9),
      height = grid::unit(3.2, "cm")),
    right_annotation = ComplexHeatmap::upset_right_annotation(
      m_shown, add_numbers = TRUE,
      numbers_gp = grid::gpar(fontsize = 8),
      axis_param = list(labels_rot = 0, gp = grid::gpar(fontsize = 8)),
      annotation_name_gp = grid::gpar(fontsize = 9),
      width = grid::unit(2.6, "cm")))
  draw_fn <- if (is.null(title)) {
    purrr::partial(ComplexHeatmap::draw, up)
  } else {
    purrr::partial(ComplexHeatmap::draw, up, column_title = title,
                   column_title_gp = grid::gpar(fontsize = 10, col = "grey30"))
  }

  stats_df <- rbind(
    data.frame(type = "set", sets = set_names, degree = NA_integer_,
               size = as.integer(ss), stringsAsFactors = FALSE),
    data.frame(type = "intersection", sets = unname(combo_sets[ord]),
               degree = as.integer(cd[ord]), size = as.integer(cs[ord]),
               stringsAsFactors = FALSE))

  test_meta <- list(
    name = "UpSet plot (ComplexHeatmap)",
    n_sets = length(ss), n_elements = nrow(bin),
    n_intersections = length(cs), n_intersections_shown = length(shown),
    largest_intersection = max(cs))

  resolved <- list(family = "upset", label = "UpSet plot")
  methods  <- .up_methods(length(ss), nrow(bin), length(cs), length(shown), set_names)
  label_in <- .up_label_inches(set_names)
  width  <- max(6.0, label_in + 0.32 * length(shown) + 1.9)
  height <- max(3.8, 0.3 * length(ss) + 2.9)
  build_script <- list(fn = .up_emit_script,
                       args = list(first_is_id = first_is_id, top_n = top_n,
                                   width = width, height = height))

  list(plot = NULL, draw = draw_fn, stats = stats_df, test_meta = test_meta,
       resolved = resolved, methods = methods, build_script = build_script,
       label = slugify(tools::file_path_sans_ext(spec$data$file)),
       df_used = df_used, clean_steps = clean_steps,
       width = width, height = height)
}

.up_methods <- function(n_sets, n_elem, n_inter, n_shown, set_names) {
  v    <- pkg_versions("ComplexHeatmap")
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  shown_txt <- if (n_shown < n_inter) {
    sprintf("The plot shows the %d largest of the %d observed intersections and each set's total size; the stats table lists them all; ",
            n_shown, n_inter)
  } else {
    sprintf("The plot shows the %d observed intersections and each set's total size; ", n_inter)
  }
  paste0(
    sprintf("Set intersections among %d sets (%s) across %d elements were drawn as an UpSet plot with ComplexHeatmap. ",
            n_sets, paste(set_names, collapse = ", "), n_elem),
    shown_txt,
    "it is descriptive, so no statistical test is applied.",
    sprintf(" Rendered in R %s with ComplexHeatmap %s.", rver, v$ComplexHeatmap))
}

.up_emit_script <- function(first_is_id, top_n, width, height, in_name, fig_stub) {
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
    "# 1, TRUE, T, yes or y is a member; any other code or an empty cell is not.",
    "bin <- raw",
    "for (cn in names(raw)) {",
    "  s <- tolower(trimws(as.character(raw[[cn]])))",
    "  bin[[cn]] <- !is.na(s) & s %in% c(\"1\", \"true\", \"t\", \"yes\", \"y\")",
    "}",
    "bin <- as.data.frame(bin[rowSums(bin) > 0, , drop = FALSE])",
    "m <- make_comb_mat(bin)",
    "cs <- comb_size(m); cd <- comb_degree(m); cn <- comb_name(m)",
    "# The largest intersections: size, then degree, then the set names.",
    "nm <- character(length(cn))",
    "for (i in seq_along(cn)) nm[i] <- paste(set_name(m)[strsplit(cn[i], \"\")[[1]] == \"1\"], collapse = \" & \")",
    "ord <- order(-cs, cd, nm, method = \"radix\")",
    sprintf("shown <- ord[seq_len(min(length(ord), %d))]", top_n),
    "m2 <- m[shown]",
    "up <- UpSet(m2, comb_order = match(cn[shown], comb_name(m2)),",
    "            set_order = match(set_name(m)[order(-set_size(m), seq_along(set_size(m)), method = \"radix\")], set_name(m2)),",
    "            top_annotation = upset_top_annotation(m2, add_numbers = TRUE, numbers_rot = 0),",
    "            right_annotation = upset_right_annotation(m2, add_numbers = TRUE,",
    "                                                      axis_param = list(labels_rot = 0)))",
    "",
    sprintf('cairo_pdf("%s.pdf", width = %.2f, height = %.2f)', fig_stub, width, height),
    "draw(up); invisible(dev.off())",
    'cat("Reproduced figure.\\n")'
  )
}
