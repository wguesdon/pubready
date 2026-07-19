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

write_data_log <- function(path, raw_path, in_name, sums, df_used, steps) {
  types <- vapply(df_used, function(v) class(v)[1], character(1))
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

build_manifest <- function(spec, resolved, in_name, orig_name, sums, n_rows,
                           test_meta, container, git_commit, created, fig_stub) {
  list(
    pubplot_version    = PUBPLOT_VERSION,
    pubplot_git_commit = git_commit,
    created_utc        = created,
    engine             = spec$engine,
    recipe             = spec$recipe,
    recipe_version     = "1",
    arguments          = list(x = spec$data$x, y = spec$data$y,
                              test = spec$test$method, paired = spec$test$paired,
                              p_adjust = spec$test$p_adjust),
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

methods_paragraph <- function(resolved, stats_df, test_meta, spec) {
  y_label <- spec$appearance$y_label %||% spec$data$y
  v <- pkg_versions(c("rstatix", "ggpubr", "ggprism"))
  rver <- paste(R.version$major, R.version$minor, sep = ".")

  norm_txt <- "Normality was assessed with the Shapiro-Wilk test."
  var_txt <- if (resolved$family == "t") {
    if (isTRUE(resolved$var_equal)) " Equal variances were assumed."
    else " Equal variances were not assumed, and the Welch correction was applied."
  } else " A non-parametric test was used because the normality assumption was not met."
  eff_txt <- if (!is.null(test_meta$effect_size)) {
    sprintf(" Effect size is reported as %s.", test_meta$effect_size$name)
  } else ""
  adj_txt <- if (!identical(spec$test$p_adjust, "none")) {
    sprintf(" P-values were adjusted for multiple comparisons using the %s method.", spec$test$p_adjust)
  } else ""
  n_txt <- if (all(c("n1", "n2") %in% names(stats_df))) {
    sprintf(" Group sizes were n = %d and n = %d.", stats_df$n1[1], stats_df$n2[1])
  } else ""

  paste0(
    sprintf("%s was compared between the two groups with a %s.",
            cap_first(y_label), resolved$label),
    " ", norm_txt, var_txt, eff_txt, adj_txt,
    " Significance was set at P < 0.05.",
    sprintf(" Analyses were performed in R %s with rstatix %s; figures were produced with ggpubr %s and ggprism %s.",
            rver, v$rstatix, v$ggpubr, v$ggprism),
    n_txt
  )
}

# A standalone, readable script that reproduces the figure from the bundled
# input using the concrete resolved parameters.
emit_script <- function(spec, resolved, in_name, fig_stub) {
  x <- spec$data$x; y <- spec$data$y
  pal <- spec$appearance$palette %||% pubplot_palette[1:2]
  pal_str <- paste0("c(", paste(sprintf('"%s"', pal), collapse = ", "), ")")
  paired <- if (isTRUE(resolved$paired)) "TRUE" else "FALSE"
  geom <- spec$appearance$geom %||% "box"
  bracket_label <- spec$appearance$bracket$label %||% "p.signif"

  if (identical(bracket_label, "p.format")) {
    plabel_line <- 'stat$plabel <- paste0("p = ", formatC(stat$p, format = "g", digits = 2))'
    label_arg <- "plabel"
  } else {
    plabel_line <- character(0)
    label_arg <- "p.signif"
  }

  test_line <- if (resolved$family == "t") {
    sprintf('stat <- rstatix::t_test(df, %s ~ %s, paired = %s, var.equal = %s)',
            y, x, paired, if (isTRUE(resolved$var_equal)) "TRUE" else "FALSE")
  } else {
    sprintf('stat <- rstatix::wilcox_test(df, %s ~ %s, paired = %s)', y, x, paired)
  }
  geom_line <- if (geom == "violin") {
    sprintf('geom_violin(aes(fill = .data[["%s"]]), trim = FALSE, width = 0.7, alpha = 0.9)', x)
  } else {
    sprintf('geom_boxplot(aes(fill = .data[["%s"]]), width = 0.6, outlier.shape = NA, alpha = 0.9)', x)
  }

  c(
    "#!/usr/bin/env Rscript",
    "# Standalone reproduction of this figure. Run inside the pinned pubplot",
    "# container (see REPRODUCE.md) from this bundle folder.",
    "suppressPackageStartupMessages({",
    "  library(ggplot2); library(ggpubr); library(rstatix); library(readr)",
    "})",
    "",
    sprintf('df <- as.data.frame(read_csv("%s", show_col_types = FALSE))', in_name),
    sprintf('df[["%s"]] <- factor(df[["%s"]])', x, x),
    sprintf('df[["%s"]] <- as.numeric(df[["%s"]])', y, y),
    "",
    test_line,
    "stat <- rstatix::add_significance(stat)",
    sprintf('stat <- rstatix::add_xy_position(stat, x = "%s")', x),
    plabel_line,
    "",
    sprintf('p <- ggplot(df, aes(x = .data[["%s"]], y = .data[["%s"]])) +', x, y),
    sprintf("  %s +", geom_line),
    sprintf('  geom_jitter(aes(fill = .data[["%s"]]), width = 0.12, size = 1.6, alpha = 0.75, shape = 21, stroke = 0.3) +', x),
    sprintf("  scale_fill_manual(values = %s) +", pal_str),
    sprintf('  stat_pvalue_manual(stat, label = "%s", tip.length = 0.01, bracket.size = 0.4) +',
            label_arg),
    sprintf('  labs(x = "%s", y = "%s") +',
            spec$appearance$x_label %||% x, spec$appearance$y_label %||% y),
    "  theme_classic(base_size = 13) +",
    '  theme(legend.position = "none")',
    "",
    sprintf('ggsave("%s.pdf", p, width = 3.8, height = 4.0)', fig_stub),
    'cat("Reproduced figure.\\n")'
  )
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
    sprintf("- Container image: `%s`", container$image %||% "localhost/pubplot:0.4.0"),
    sprintf("- Image id: `%s`", container$image_id %||% "unknown"),
    sprintf("- Image digest: `%s`", if (nzchar(container$digest %||% "")) container$digest else "n/a (local build)"),
    sprintf("- pubplot commit: `%s`", git_commit),
    "",
    "## 2. Run the standalone script",
    "",
    "From inside this folder:",
    "",
    "```bash",
    sprintf("podman run --rm -v \"$PWD\":/work -w /work \\"),
    sprintf("  %s Rscript %s", container$image %||% "localhost/pubplot:0.4.0", script_name),
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

write_bundle <- function(spec, resolved, df_used, raw_input_path, plot, stats_df,
                         test_meta, data_log_steps, out_root, container,
                         git_commit, methods_text, build_script, base_label,
                         draw = NULL, extra_inputs = NULL,
                         fig_width = 3.8, fig_height = 4.0, stamp = NULL) {
  # A fixed --stamp gives deterministic bundle and file names (used for the
  # committed reference outputs); otherwise names carry a real UTC timestamp.
  ts      <- if (!is.null(stamp) && nzchar(stamp)) stamp else timestamp_utc()
  env_created <- Sys.getenv("PUBPLOT_CREATED", "")
  created <- if (nzchar(env_created)) env_created else iso_utc()
  base     <- slugify(base_label)
  bdir     <- file.path(out_root, sprintf("%s_%s_%s", spec$recipe, base, ts))
  dir.create(bdir, recursive = TRUE, showWarnings = FALSE)

  fig_stub <- sprintf("figure_%s_%s", base, ts)
  save_figure(plot, bdir, fig_stub, width = fig_width, height = fig_height, draw = draw)

  readr::write_csv(stats_df, file.path(bdir, sprintf("stats_%s.csv", ts)))

  in_ext  <- tools::file_ext(raw_input_path)
  in_name <- sprintf("input_%s.%s", base, if (nzchar(in_ext)) in_ext else "csv")
  file.copy(raw_input_path, file.path(bdir, in_name), overwrite = TRUE)
  sums <- file_checksums(raw_input_path)

  # Copy any secondary input files (e.g. a heatmap annotation table) with their
  # original names so the standalone script and REPRODUCE step can find them.
  for (ex in extra_inputs) {
    if (!is.null(ex) && nzchar(ex) && file.exists(ex)) {
      file.copy(ex, file.path(bdir, basename(ex)), overwrite = TRUE)
    }
  }

  write_data_log(file.path(bdir, "data_log.md"),
                 raw_input_path, in_name, sums, df_used, data_log_steps)

  cfg <- spec
  cfg$data$file <- in_name
  write_config(cfg, file.path(bdir, "plot_config.yaml"))

  writeLines(capture.output(utils::sessionInfo()), file.path(bdir, "session_info.txt"))

  manifest <- build_manifest(spec, resolved, in_name, basename(raw_input_path),
                             sums, nrow(df_used), test_meta, container,
                             git_commit, created, fig_stub)
  writeLines(
    jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 6),
    file.path(bdir, sprintf("manifest_%s.json", ts))
  )

  writeLines(methods_text, file.path(bdir, sprintf("methods_%s.md", ts)))

  script_name <- sprintf("script_%s.R", ts)
  writeLines(build_script(in_name, fig_stub), file.path(bdir, script_name))

  writeLines(reproduce_md(container, git_commit, script_name),
             file.path(bdir, "REPRODUCE.md"))

  bdir
}
