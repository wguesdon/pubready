# The artifact bundle writer. One self-contained folder per figure:
# figure, standalone script, stats, input copy, data log, plot config,
# session info, manifest, methods paragraph, and REPRODUCE.md.

timestamp_utc <- function() format(Sys.time(), "%Y_%m_%d_%H%M%S", tz = "UTC")
iso_utc       <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

slugify <- function(s) {
  s <- tolower(as.character(s))
  s <- gsub("[^a-z0-9]+", "_", s)
  gsub("^_|_$", "", s)
}

# The OS of the container, from /etc/os-release (e.g. "Ubuntu 22.04.5 LTS").
os_pretty <- function() {
  f <- "/etc/os-release"
  if (!file.exists(f)) return(NA_character_)
  ln <- grep("^PRETTY_NAME=", readLines(f, warn = FALSE), value = TRUE)
  if (length(ln) == 0) return(NA_character_)
  gsub('"', "", sub("^PRETTY_NAME=", "", ln[1]))
}

# Packages the generated script loads directly, via library()/require()/pkg::.
script_direct_pkgs <- function(script_lines) {
  txt <- paste(script_lines, collapse = "\n")
  grab <- function(re) {
    m <- regmatches(txt, gregexpr(re, txt, perl = TRUE))[[1]]
    if (length(m) == 0) return(character(0))
    sub(re, "\\1", m, perl = TRUE)
  }
  libs <- grab("(?:library|require)\\(([A-Za-z][A-Za-z0-9.]*)")
  ns   <- grab("([A-Za-z][A-Za-z0-9.]*)::")
  unique(c(libs, ns))
}

# Write the package versions that produced this figure, next to the figure and
# script. The packages the script loads directly come first (the answer to "what
# did this figure use"), then the complete installed set, since R sessionInfo()
# lists only the loaded packages.
write_packages <- function(path, script_lines = NULL) {
  ip   <- installed.packages()[, c("Package", "Version"), drop = FALSE]
  ver  <- stats::setNames(ip[, "Version"], ip[, "Package"])
  rver <- paste(R.version$major, R.version$minor, sep = ".")
  lines <- c("# Package versions for this figure",
             sprintf("# R %s | %s", rver, os_pretty() %||% "unknown OS"), "")
  direct <- if (!is.null(script_lines)) script_direct_pkgs(script_lines) else character(0)
  direct <- direct[direct %in% names(ver)]
  if (length(direct)) {
    lines <- c(lines, "## Used by this script",
               sprintf("%s==%s", direct, ver[direct]), "")
  }
  all_sorted <- ip[order(tolower(ip[, "Package"])), , drop = FALSE]
  lines <- c(lines, "## Full environment (every installed package)", "",
             sprintf("%s==%s", all_sorted[, "Package"], all_sorted[, "Version"]))
  writeLines(lines, path)
}

save_figure <- function(plot, dir, stub, width = 3.8, height = 4.0, draw = NULL) {
  # draw is a function that renders to the active device (base-graphics style,
  # e.g. ComplexHeatmap); otherwise plot is a ggplot/grob saved with ggsave.
  if (!is.null(draw)) {
    grDevices::cairo_pdf(file.path(dir, paste0(stub, ".pdf")), width = width, height = height)
    draw(); grDevices::dev.off()
    grDevices::png(file.path(dir, paste0(stub, ".png")), width = width, height = height,
                   units = "in", res = 300, bg = "white")
    draw(); grDevices::dev.off()
    grDevices::svg(file.path(dir, paste0(stub, ".svg")), width = width, height = height)
    draw(); grDevices::dev.off()
    return(invisible())
  }
  ggplot2::ggsave(file.path(dir, paste0(stub, ".pdf")), plot,
                  width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(dir, paste0(stub, ".png")), plot,
                  width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(stub, ".svg")), plot,
                  width = width, height = height)
}

# The read rule of io.R, as the data log states it. The Python engine writes the
# same text.
READ_RULE_TEXT <- paste(
  "Each cell was read as text, and the spaces around each value were removed.",
  "Only an empty cell or \"NA\" was read as missing. A column was read as numeric",
  "when each present cell is a number with \".\" as the decimal mark, and as",
  "logical when each present cell is TRUE or FALSE.")

# NaN becomes NA, so a missing value has one form in the stats table.
nan_to_na <- function(v) {
  if (is.double(v)) v[is.nan(v)] <- NA_real_
  v
}

# Write a stats table in the form that both engines share: NA for a missing
# value, TRUE and FALSE, Inf and -Inf, and a whole number with no decimal part.
write_stats_csv <- function(stats_df, path) {
  stats_df[] <- lapply(stats_df, nan_to_na)
  readr::write_csv(stats_df, path, na = "NA")
}

write_data_log <- function(path, raw_path, in_name, sums, df_used, steps,
                           read_steps = character(0)) {
  types <- vapply(df_used, function(v) class(v)[1], character(1))
  if (length(read_steps) == 0) read_steps <- "No value changed at read."
  lines <- c(
    "# Data provenance log",
    "",
    sprintf("- Original file: `%s`", basename(raw_path)),
    sprintf("- Stored in bundle as: `%s`", in_name),
    sprintf("- md5: `%s`", sums$md5),
    sprintf("- sha256: `%s`", sums$sha256),
    sprintf("- Analyzed rows: %d", nrow(df_used)),
    sprintf("- Columns: %d", ncol(df_used)),
    "",
    "## Read rule",
    "",
    READ_RULE_TEXT,
    "",
    paste0("- ", read_steps),
    "",
    "## Column types (analyzed data)",
    "",
    paste0("- `", names(types), "`: ", types),
    "",
    "## Cleaning steps",
    "",
    paste0("- ", steps)
  )
  writeLines(lines, path)
}

# The plot options that the manifest records, in this order, with the spec
# field that holds each one. Each option is recorded, with null when the run did
# not set it. The Python engine records the same options in the same order.
MANIFEST_DATA_ARGS <- c("x", "y", "fill", "facet", "time", "event", "covariates",
                        "annotation", "id", "group", "label", "positive", "se",
                        "ci_low", "ci_high", "measure", "markers", "ref", "cause",
                        "levels", "visit_levels", "features")
MANIFEST_APPEARANCE_ARGS <- c(scale = "scale", cluster = "cluster",
                              fc_cutoff = "fc_cutoff", p_cutoff = "p_cutoff",
                              top_n = "top_n", pd_threshold = "pd_threshold",
                              pr_threshold = "pr_threshold", geom = "geom",
                              theme = "theme", xlab = "x_label", ylab = "y_label",
                              title = "title", palette = "palette")
MANIFEST_LIST_ARGS <- c("covariates", "markers", "levels", "visit_levels",
                        "features", "palette")

# The arguments block of the manifest: every plot option, named as the command
# line names it.
manifest_arguments <- function(spec) {
  out <- list()
  for (k in MANIFEST_DATA_ARGS) out[k] <- list(spec$data[[k]])
  out["test"] <- list(spec$test$method %||% "auto")
  out["paired"] <- list(isTRUE(spec$test$paired))
  for (k in names(MANIFEST_APPEARANCE_ARGS)) {
    out[k] <- list(spec$appearance[[MANIFEST_APPEARANCE_ARGS[[k]]]])
  }
  # A list option stays a JSON array when it holds one item.
  for (k in MANIFEST_LIST_ARGS) {
    if (!is.null(out[[k]])) out[[k]] <- I(as.character(unlist(out[[k]])))
  }
  out
}

build_manifest <- function(spec, resolved, in_name, orig_name, sums, n_rows,
                           test_meta, container, git_commit, created, fig_stub) {
  list(
    pubready_version    = PUBREADY_VERSION,
    pubready_git_commit = git_commit,
    created_utc        = created,
    engine             = spec$engine,
    recipe             = spec$recipe,
    recipe_version     = "1",
    arguments          = manifest_arguments(spec),
    input = list(file = in_name, original_name = orig_name,
                 md5 = sums$md5, sha256 = sums$sha256, n_rows = n_rows),
    statistical_test = test_meta,
    container = container,
    environment = list(
      language         = "R",
      language_version = paste(R.version$major, R.version$minor, sep = "."),
      packages         = pkg_versions(c("ggplot2", "ggpubr", "rstatix", "ggprism"))
    ),
    outputs = list(paste0(fig_stub, ".pdf"), paste0(fig_stub, ".png"),
                   paste0(fig_stub, ".svg"))
  )
}

# The JSON text of a manifest. Each number keeps its full precision, because
# some values, such as the standard-curve intercept, are only in the manifest.
manifest_json <- function(manifest) {
  jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null",
                   digits = NA)
}

# The distribution geom of an emitted comparison script, as add_dist_geom() draws
# it: box | violin | bar | raincloud, plus the points unless show_points is
# FALSE. `prep` holds the lines that set the raincloud offsets on df before the
# plot; `layers` holds the ggplot2 layer lines, each ending in " +".
script_geom_lines <- function(x, geom = "box", show_points = TRUE) {
  fill <- sprintf('aes(fill = .data[["%s"]])', x)
  layers <- switch(geom,
    violin = sprintf("  geom_violin(%s, trim = FALSE, width = 0.7, alpha = 0.9) +", fill),
    bar = c(
      sprintf('  stat_summary(%s, fun = mean, geom = "bar", width = 0.6, alpha = 0.9, colour = "black", linewidth = 0.3) +', fill),
      '  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +'),
    raincloud = c(
      sprintf("  ggdist::stat_halfeye(%s, adjust = 0.6, width = 0.6, .width = 0, justification = -0.2, point_colour = NA, alpha = 0.85) +", fill),
      sprintf("  geom_boxplot(%s, width = 0.12, outlier.shape = NA, alpha = 0.9) +", fill)),
    sprintf("  geom_boxplot(%s, width = 0.6, outlier.shape = NA, alpha = 0.9) +", fill))
  prep <- character(0)
  if (isTRUE(show_points)) {
    # The golden-ratio offsets of rain_offsets() by row in each group: the rain
    # of a raincloud falls left of the box, other points spread +/- 0.12.
    lim <- if (geom == "raincloud") c(-0.18, 0.12) else c(-0.12, 0.24)
    prep <- c(
      sprintf('point_k <- ave(seq_len(nrow(df)), df[["%s"]], FUN = seq_along)', x),
      sprintf('df$.point_x <- as.integer(df[["%s"]]) + %s + %s * ((point_k * 0.6180339887498949) %%%% 1)',
              x, format(lim[1]), format(lim[2])))
    layers <- c(layers, sprintf(
      '  geom_point(aes(x = .point_x, fill = .data[["%s"]]), size = 1.5, alpha = 0.7, shape = 21, stroke = 0.3) +', x))
  }
  list(prep = prep, layers = layers)
}

# A recipe returns its script builder either as a function of (in_name,
# fig_stub) or as list(fn = <function>, args = list(...)). The list form lets a
# recipe name a module-level emitter and its arguments instead of a closure.
call_build_script <- function(build_script, in_name, fig_stub) {
  if (is.function(build_script)) return(build_script(in_name, fig_stub))
  do.call(build_script$fn, c(build_script$args,
                             list(in_name = in_name, fig_stub = fig_stub)))
}

reproduce_md <- function(container, git_commit, script_name) {
  c(
    "# Reproduce this figure",
    "",
    "This bundle is self-contained. To regenerate the figure in the exact",
    "environment it was made in:",
    "",
    "## 1. Get the pinned environment",
    "",
    sprintf("- Container image: `%s`", container$image %||% "localhost/pubready:4.0.0"),
    sprintf("- Image id: `%s`", container$image_id %||% "unknown"),
    sprintf("- Image digest: `%s`", if (nzchar(container$digest %||% "")) container$digest else "n/a (local build)"),
    sprintf("- pubready commit: `%s`", git_commit),
    "",
    "## 2. Run the standalone script",
    "",
    "From inside this folder:",
    "",
    "```bash",
    "podman run --rm --security-opt label=disable -v \"$PWD\":/work -w /work \\",
    sprintf("  %s Rscript %s", container$image %||% "localhost/pubready:4.0.0", script_name),
    "```",
    "",
    "The script reads the bundled input copy, reruns the same test, and redraws",
    "the figure. The stats are recomputed, so the figure and the numbers cannot",
    "drift apart.",
    "",
    "See `manifest_*.json` for the full package versions and the test result,",
    "and `session_info.txt` for the complete environment."
  )
}

# One column of a report table as escaped text. The column becomes text before
# the missing values become "", because a factor has no "" level. A "|" becomes
# "\|", as in the Python engine.
report_escape_column <- function(column) {
  column <- as.character(nan_to_na(column))
  column[is.na(column)] <- ""
  gsub("|", "\\|", column, fixed = TRUE)
}

report_markdown_row <- function(row) {
  paste0("| ", paste(row, collapse = " | "), " |")
}

report_markdown_table <- function(frame, limit = 100L) {
  shown <- utils::head(frame, limit)
  shown[] <- lapply(shown, report_escape_column)
  header <- paste0("| ", paste(names(shown), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(shown)), collapse = " | "), " |")
  rows <- apply(shown, 1L, report_markdown_row)
  paste(c(header, rule, rows), collapse = "\n")
}

report_git_dirty <- function() {
  value <- Sys.getenv("PUBREADY_GIT_DIRTY", unset = "not_recorded")
  if (value == "true") return(TRUE)
  if (value == "false") return(FALSE)
  value
}

report_provenance_text <- function(value) {
  if (length(value) == 1L) return(as.character(value))
  jsonlite::toJSON(value, auto_unbox = TRUE, null = "null")
}

build_report_provenance <- function(spec, container, git_commit, created,
                                    in_name, sums) {
  list(
    repository = "pubready",
    repository_version = PUBREADY_VERSION,
    git_commit = git_commit,
    git_dirty = report_git_dirty(),
    llm_provider = Sys.getenv("PUBREADY_LLM_PROVIDER", unset = "none"),
    llm_model = Sys.getenv("PUBREADY_LLM_MODEL", unset = "none"),
    agent_harness = Sys.getenv("PUBREADY_AGENT_HARNESS", unset = "manual"),
    agent_harness_version = Sys.getenv(
      "PUBREADY_AGENT_HARNESS_VERSION",
      unset = "not_recorded"
    ),
    container_image = container$image,
    container_digest = container$digest %||% "not_recorded",
    input = list(file = in_name, md5 = sums$md5, sha256 = sums$sha256),
    command = paste("figkit plot --recipe", spec$recipe),
    analysis_options = spec,
    random_seeds = spec$seed %||% "not_applicable",
    session_record = "session_info.txt",
    created_utc = created,
    review_state = "not_applicable",
    reviewer = "none",
    quarto_version = "1.11.4"
  )
}

report_callout <- function(provenance) {
  if (provenance$llm_model == "none") {
    return("A manual command generated this report without an LLM model.")
  }
  paste0(
    "An LLM generated this report through the `",
    provenance$agent_harness,
    "` agent harness."
  )
}

write_quarto_report <- function(bdir, spec, stats_df, methods_text, fig_stub,
                                provenance) {
  jsonlite::write_json(
    provenance,
    file.path(bdir, "report_provenance.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    digits = NA
  )
  provenance_frame <- data.frame(
    field = names(provenance),
    value = vapply(provenance, report_provenance_text, character(1)),
    stringsAsFactors = FALSE
  )
  model <- jsonlite::toJSON(provenance$llm_model, auto_unbox = TRUE)
  harness <- jsonlite::toJSON(provenance$agent_harness, auto_unbox = TRUE)
  lines <- c(
    "---",
    'title: "Publication figure report"',
    "author:",
    paste0("  - name: ", model),
    "    affiliations:",
    paste0("      - name: ", harness),
    "format:",
    "  html:",
    "    embed-resources: true",
    "    toc: true",
    "---",
    "",
    "::: {.callout-important}",
    report_callout(provenance),
    "No separate human scientific review is part of this fixed recipe.",
    "The report records the selected recipe and does not claim a new method.",
    ":::",
    "",
    "## Figure",
    "",
    paste0("![Publication figure](", fig_stub, ".png)"),
    "",
    "## Statistical results",
    "",
    report_markdown_table(stats_df),
    "",
    "## Methods",
    "",
    methods_text,
    "",
    "## Plot configuration",
    "",
    "The complete plot configuration is in `plot_config.yaml`.",
    "",
    "## Provenance",
    "",
    report_markdown_table(provenance_frame),
    "",
    "The report source reads only files in this bundle."
  )
  writeLines(lines, file.path(bdir, "report.qmd"), useBytes = TRUE)
  old_directory <- getwd()
  on.exit(setwd(old_directory), add = TRUE)
  setwd(bdir)
  status <- system2(
    "quarto",
    c("render", "report.qmd", "--to", "html", "--output", "report.html")
  )
  if (status != 0L || !file.exists("report.html")) {
    stop("Quarto did not create report.html.")
  }
  c("report.qmd", "report.html", "report_provenance.json")
}

report_output_record <- function(path) {
  list(
    file = basename(path),
    bytes = unname(file.info(path)$size),
    sha256 = digest::digest(path, algo = "sha256", file = TRUE)
  )
}

write_bundle <- function(spec, resolved, df_used, raw_input_path, plot, stats_df,
                         test_meta, data_log_steps, out_root, container,
                         git_commit, methods_text, build_script, base_label,
                         draw = NULL, extra_inputs = NULL,
                         fig_width = 3.8, fig_height = 4.0, qc = NULL, stamp = NULL,
                         read_steps = character(0)) {
  # A fixed --stamp gives deterministic bundle and file names (used for the
  # committed reference outputs); otherwise names carry a real UTC timestamp.
  ts      <- if (!is.null(stamp) && nzchar(stamp)) stamp else timestamp_utc()
  env_created <- Sys.getenv("PUBREADY_CREATED", "")
  created <- if (nzchar(env_created)) env_created else iso_utc()
  base     <- slugify(base_label)
  bdir     <- file.path(out_root, sprintf("%s_%s_%s", spec$recipe, base, ts))
  dir.create(bdir, recursive = TRUE, showWarnings = FALSE)

  fig_stub <- sprintf("figure_%s_%s", base, ts)
  save_figure(plot, bdir, fig_stub, width = fig_width, height = fig_height, draw = draw)

  # Normality QC panel: a QQ per checked quantity, so the reader can judge the
  # assumption behind the test rather than trust a single Shapiro-Wilk p-value.
  if (!is.null(qc)) {
    qcres <- qc_normality(qc)
    d <- list(width = qcres$width, height = qcres$height)
    ggplot2::ggsave(file.path(bdir, sprintf("qc_normality_%s.png", ts)), qcres$plot,
                    width = d$width, height = d$height, dpi = 200, bg = "white")
    ggplot2::ggsave(file.path(bdir, sprintf("qc_normality_%s.pdf", ts)), qcres$plot,
                    width = d$width, height = d$height, device = grDevices::cairo_pdf)
  }

  write_stats_csv(stats_df, file.path(bdir, sprintf("stats_%s.csv", ts)))

  in_ext  <- tools::file_ext(raw_input_path)
  in_name <- sprintf("input_%s.%s", base, if (nzchar(in_ext)) in_ext else "csv")
  file.copy(raw_input_path, file.path(bdir, in_name), overwrite = TRUE)
  sums <- file_checksums(raw_input_path)
  script_lines <- call_build_script(build_script, in_name, fig_stub)

  # Copy any secondary input files (e.g. a heatmap annotation table) with their
  # original names so the standalone script and REPRODUCE step can find them.
  for (ex in extra_inputs) {
    if (!is.null(ex) && nzchar(ex) && file.exists(ex)) {
      file.copy(ex, file.path(bdir, basename(ex)), overwrite = TRUE)
    }
  }

  write_data_log(file.path(bdir, "data_log.md"),
                 raw_input_path, in_name, sums, df_used, data_log_steps, read_steps)

  # The config names the input copy and the annotation copy by the bare file
  # name; render finds them in the folder of the config.
  cfg <- spec
  cfg$data$file <- in_name
  if (!is.null(cfg$data$annotation) && file.exists(file.path(bdir, basename(cfg$data$annotation)))) {
    cfg$data$annotation <- basename(cfg$data$annotation)
  }
  write_config(cfg, file.path(bdir, "plot_config.yaml"))

  writeLines(capture.output(utils::sessionInfo()), file.path(bdir, "session_info.txt"))

  write_packages(file.path(bdir, sprintf("packages_%s.txt", ts)), script_lines)

  manifest <- build_manifest(spec, resolved, in_name, basename(raw_input_path),
                             sums, nrow(df_used), test_meta, container,
                             git_commit, created, fig_stub)
  manifest$outputs <- c(manifest$outputs, list(sprintf("packages_%s.txt", ts)))
  if (!is.null(qc)) {
    manifest$outputs <- c(manifest$outputs,
                          list(sprintf("qc_normality_%s.png", ts),
                               sprintf("qc_normality_%s.pdf", ts)))
  }
  writeLines(methods_text, file.path(bdir, sprintf("methods_%s.md", ts)))

  script_name <- sprintf("script_%s.R", ts)
  writeLines(script_lines, file.path(bdir, script_name))

  writeLines(reproduce_md(container, git_commit, script_name),
             file.path(bdir, "REPRODUCE.md"))

  provenance <- build_report_provenance(
    spec,
    container,
    git_commit,
    created,
    in_name,
    sums
  )
  report_outputs <- write_quarto_report(
    bdir,
    spec,
    stats_df,
    methods_text,
    fig_stub,
    provenance
  )
  manifest$outputs <- c(manifest$outputs, as.list(report_outputs))
  output_paths <- file.path(bdir, unlist(manifest$outputs))
  output_paths <- output_paths[file.exists(output_paths)]
  manifest$output_checksums <- lapply(output_paths, report_output_record)
  writeLines(manifest_json(manifest), file.path(bdir, sprintf("manifest_%s.json", ts)))

  bdir
}
