# Drive one recipe end to end: read, clean, run the recipe, write the bundle.

container_info <- function() {
  list(
    image          = Sys.getenv("PUBREADY_IMAGE", "localhost/pubready:1.0.0"),
    image_id       = Sys.getenv("PUBREADY_IMAGE_ID", ""),
    digest         = Sys.getenv("PUBREADY_IMAGE_DIGEST", ""),
    podman_version = Sys.getenv("PUBREADY_PODMAN_VERSION", "")
  )
}

run_recipe <- function(spec, raw_input, out_root, sheet = NULL, stamp = NULL) {
  df <- read_tidy(raw_input, sheet)

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
    stamp          = stamp
  )
  cat(sprintf("Wrote bundle: %s\n", bdir))
  invisible(bdir)
}

# Generic data inspection used by the `inspect` command.
inspect_data <- function(df) {
  cols <- lapply(names(df), function(nm) {
    v <- df[[nm]]
    is_num <- is.numeric(v)
    out <- list(
      name      = nm,
      type      = class(v)[1],
      n_missing = sum(is.na(v)),
      n_unique  = length(unique(v[!is.na(v)]))
    )
    if (!is_num && out$n_unique <= 20 && out$n_unique > 0) {
      out$levels <- as.list(sort(unique(as.character(v[!is.na(v)]))))
    }
    if (is_num && any(is.finite(v))) {
      out$min  <- min(v, na.rm = TRUE)
      out$max  <- max(v, na.rm = TRUE)
      out$mean <- mean(v, na.rm = TRUE)
    }
    out
  })
  list(n_rows = nrow(df), n_cols = ncol(df), columns = cols)
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

print_inspect <- function(info) {
  cat(sprintf("Rows: %d   Columns: %d\n\n", info$n_rows, info$n_cols))
  group_like <- character(0)
  numeric_cols <- character(0)
  for (i in seq_along(info$columns)) {
    c <- info$columns[[i]]
    if (!is.null(c$levels) && c$n_unique >= 2 && c$n_unique <= 6) group_like <- c(group_like, c$name)
    if (!is.null(c$mean)) numeric_cols <- c(numeric_cols, c$name)
    if (i > .inspect_max_cols) next
    line <- sprintf("  %-20s %-10s  missing=%d  unique=%d",
                    c$name, c$type, c$n_missing, c$n_unique)
    if (!is.null(c$levels)) {
      line <- paste0(line, "  levels: ", paste(unlist(c$levels), collapse = ", "))
    }
    if (!is.null(c$mean)) {
      line <- paste0(line, sprintf("  [%.3g .. %.3g, mean %.3g]", c$min, c$max, c$mean))
    }
    cat(line, "\n")
  }
  hidden <- length(info$columns) - .inspect_max_cols
  if (hidden > 0) {
    cat(sprintf("  ... and %d more column(s); use --format json for every column.\n", hidden))
  }
  cat("\n")
  if (length(group_like)) cat("Candidate grouping columns (x): ", .inspect_names(group_like), "\n")
  if (length(numeric_cols)) cat("Candidate outcome columns  (y): ", .inspect_names(numeric_cols), "\n")
}
