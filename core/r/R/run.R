# Drive one recipe end to end: read, clean, run the recipe, write the bundle.

container_info <- function() {
  list(
    image          = Sys.getenv("PUBREADY_IMAGE", "localhost/pubready:3.0.0"),
    image_id       = Sys.getenv("PUBREADY_IMAGE_ID", ""),
    digest         = Sys.getenv("PUBREADY_IMAGE_DIGEST", ""),
    podman_version = Sys.getenv("PUBREADY_PODMAN_VERSION", "")
  )
}

run_recipe <- function(spec, raw_input, out_root, sheet = NULL, stamp = NULL) {
  rd <- read_tidy_logged(raw_input, sheet)
  df <- rd$df

  fn_name <- paste0("recipe_", spec$recipe)
  if (!exists(fn_name, mode = "function")) {
    stop(sprintf("unknown recipe '%s'", spec$recipe))
  }
  # Each recipe cleans and validates its own inputs, and returns the data it
  # actually used plus the cleaning steps, so the bundle records them.
  res <- match.fun(fn_name)(df, spec)

  bdir <- write_bundle(
    spec           = spec,
    resolved       = res$resolved,
    df_used        = res$df_used,
    raw_input_path = raw_input,
    plot           = res$plot,
    stats_df       = res$stats,
    test_meta      = res$test_meta,
    data_log_steps = res$clean_steps,
    out_root       = out_root,
    container      = container_info(),
    git_commit     = Sys.getenv("PUBREADY_GIT_COMMIT", "unknown"),
    methods_text   = res$methods,
    build_script   = res$build_script,
    base_label     = res$label %||% spec$data$y %||% spec$data$time %||%
                     spec$data$x %||% spec$recipe,
    draw           = res$draw,
    extra_inputs   = res$extra_inputs,
    fig_width      = res$width %||% 3.8,
    fig_height     = res$height %||% 4.0,
    qc             = res$qc,
    stamp          = stamp,
    read_steps     = rd$steps
  )
  cat(sprintf("Wrote bundle: %s\n", bdir))
  invisible(bdir)
}

# Data inspection for the `inspect` command. The Python engine (run.py) applies
# the same rules and prints the same words.
#
# - Each column gets one of the type names numeric, text, logical, date, other.
# - A text or logical column with 1 to 20 distinct values, and a numeric column
#   with 2 to 6, lists its levels in the order of level_order(), each with its
#   count of rows.
# - A column with 2 to 6 levels is a candidate group. A numeric one is marked as
#   a numeric code.
# - A column is a possible subject id when each of its values occurs at most once
#   in each level of another column, half of its values or more occur in two of
#   those levels or more, and it holds more distinct values than that column.
#   The other column holds 2 to 20 levels and repeats each level at least twice
#   on average; the first 10 such columns are checked.
# - A numeric column with more than 6 distinct values that is not a possible
#   subject id is a candidate outcome.
INSPECT_GROUP_MAX <- 6
INSPECT_LEVELS_MAX <- 20
INSPECT_BY_MAX <- 10

# The type name of a column in the inspect report.
inspect_type <- function(v) {
  if (is.logical(v)) return("logical")
  if (is.numeric(v)) return("numeric")
  if (is.character(v) || is.factor(v)) return("text")
  if (inherits(v, c("Date", "POSIXt"))) return("date")
  "other"
}

# The text labels of a column for the level and id rules: level_labels() for a
# number, "TRUE" and "FALSE" for a logical value, NA for a missing value.
inspect_labels <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (inherits(v, c("Date", "POSIXt"))) v <- format(v)
  level_labels(v)
}

# True when a column of the given type and distinct count lists its levels.
inspect_has_levels <- function(type, n_unique) {
  (type %in% c("text", "logical") && n_unique >= 1 && n_unique <= INSPECT_LEVELS_MAX) ||
    (type == "numeric" && n_unique >= 2 && n_unique <= INSPECT_GROUP_MAX)
}

# One column of the inspect report.
inspect_column <- function(nm, df) {
  v <- df[[nm]]
  type <- inspect_type(v)
  lab <- inspect_labels(v)
  present <- lab[!is.na(lab)]
  n_unique <- length(unique(present))
  out <- list(name = nm, type = type, n_missing = sum(is.na(v)), n_unique = n_unique)
  if (inspect_has_levels(type, n_unique)) {
    lv <- level_order(present)
    out$levels <- as.list(lv)
    out$level_counts <- as.list(as.integer(table(factor(present, levels = lv))))
  }
  if (type == "numeric" && any(is.finite(v))) {
    out$min  <- min(v, na.rm = TRUE)
    out$max  <- max(v, na.rm = TRUE)
    out$mean <- mean(v, na.rm = TRUE)
  }
  out
}

# The id record of column `a` across the levels of column `b`, or NULL when `a`
# is not a possible subject id there.
inspect_id_candidate <- function(a_lab, b_lab, a_name, b_name) {
  keep <- !is.na(a_lab) & !is.na(b_lab)
  a <- a_lab[keep]
  b <- b_lab[keep]
  n_a <- length(unique(a))
  if (n_a < 2 || n_a <= length(unique(b))) return(NULL)
  if (anyDuplicated(data.frame(a, b))) return(NULL)
  tab <- table(a, b)
  n_levels <- rowSums(tab > 0)
  if (mean(n_levels >= 2) < 0.5) return(NULL)
  list(column = a_name, across = b_name, n_subjects = nrow(tab),
       complete = all(n_levels == ncol(tab)))
}

inspect_data <- function(df) {
  cols <- lapply(names(df), inspect_column, df = df)
  labs <- lapply(df, inspect_labels)
  n_unique <- vapply(cols, `[[`, integer(1), "n_unique")
  by <- which(n_unique >= 2 & n_unique <= INSPECT_LEVELS_MAX & n_unique * 2 <= nrow(df))
  by <- utils::head(by, INSPECT_BY_MAX)
  ids <- list()
  for (i in seq_along(cols)) {
    for (j in by) {
      if (i == j) next
      cand <- inspect_id_candidate(labs[[i]], labs[[j]], cols[[i]]$name, cols[[j]]$name)
      if (!is.null(cand)) {
        ids[[length(ids) + 1]] <- cand
        break
      }
    }
  }
  id_names <- vapply(ids, `[[`, character(1), "column")
  groups <- character(0)
  outcomes <- character(0)
  for (cl in cols) {
    if (!is.null(cl$levels) && cl$n_unique >= 2 && cl$n_unique <= INSPECT_GROUP_MAX) {
      groups <- c(groups, cl$name)
    }
    if (cl$type == "numeric" && cl$n_unique > INSPECT_GROUP_MAX && !cl$name %in% id_names) {
      outcomes <- c(outcomes, cl$name)
    }
  }
  list(n_rows = nrow(df), n_cols = ncol(df), columns = cols,
       group_candidates = as.list(groups), outcome_candidates = as.list(outcomes),
       id_candidates = ids)
}

# A wide table (a matrix of genes, for example) would print one line per column
# and a candidate list hundreds of names long. The text report shows the first
# columns in full and summarises the rest; --format json keeps every column.
.inspect_max_cols <- 25
.inspect_max_names <- 20

.inspect_names <- function(nms) {
  if (length(nms) <= .inspect_max_names) return(paste(nms, collapse = ", "))
  sprintf("%s, and %d more", paste(nms[seq_len(.inspect_max_names)], collapse = ", "),
          length(nms) - .inspect_max_names)
}

# The levels of a column as the text report prints them: each level in quotes
# with its count of rows.
inspect_levels_text <- function(cl) {
  paste(sprintf("\"%s\" (%d)", unlist(cl$levels), unlist(cl$level_counts)), collapse = ", ")
}

# One line of the candidate list of the text report.
inspect_candidate_line <- function(label, names) {
  cat(sprintf("%-33s%s\n", label, .inspect_names(names)))
}

print_inspect <- function(info) {
  cat(sprintf("Rows: %d   Columns: %d\n\n", info$n_rows, info$n_cols))
  for (i in seq_len(min(length(info$columns), .inspect_max_cols))) {
    cl <- info$columns[[i]]
    line <- sprintf("  %-20s %-8s  missing=%d  unique=%d",
                    cl$name, cl$type, cl$n_missing, cl$n_unique)
    if (!is.null(cl$levels)) {
      line <- paste0(line, "  levels: ", inspect_levels_text(cl))
    }
    if (!is.null(cl$mean)) {
      line <- paste0(line, sprintf("  [%.3g .. %.3g, mean %.3g]", cl$min, cl$max, cl$mean))
    }
    cat(line, "\n", sep = "")
  }
  hidden <- length(info$columns) - .inspect_max_cols
  if (hidden > 0) {
    cat(sprintf("  ... and %d more column(s); use --format json for every column.\n", hidden))
  }
  cat("\n")
  types <- vapply(info$columns, `[[`, character(1), "type")
  names(types) <- vapply(info$columns, `[[`, character(1), "name")
  groups <- unlist(info$group_candidates)
  if (length(groups)) {
    is_code <- types[groups] == "numeric"
    groups[is_code] <- paste0(groups[is_code], " (numeric code)")
    inspect_candidate_line("Candidate grouping columns (x):", groups)
  }
  if (length(info$outcome_candidates)) {
    inspect_candidate_line("Candidate outcome columns (y):", unlist(info$outcome_candidates))
  }
  if (length(info$id_candidates)) {
    ids <- vapply(info$id_candidates, inspect_id_text, character(1))
    inspect_candidate_line("Possible subject id columns:", ids)
  }
}

# One possible subject id as the text report prints it.
inspect_id_text <- function(cand) {
  sprintf("%s (repeats across the levels of %s)", cand$column, cand$across)
}
