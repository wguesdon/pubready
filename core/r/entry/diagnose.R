#!/usr/bin/env Rscript
# Entry: figkit diagnose. Show the assumption evidence behind the test choice so
# a scientist can confirm or override it. Writes a normality QQ panel per checked
# quantity plus a machine-readable report, and prints the recommended test. It
# does not draw the final figure; run `figkit plot` once the test is agreed.
suppressPackageStartupMessages(library(optparse))

source(file.path(Sys.getenv("PUBREADY_CORE", "/opt/pubready/core"), "r", "bootstrap.R"))

# diagnose accepts every option of figkit plot, so the recipe sees the same spec.
opt <- parse_args(OptionParser(option_list = plot_option_list()))
for (req in c("recipe", "data")) if (is.null(opt[[req]])) stop(sprintf("missing required --%s", req))

spec <- spec_from_opt(opt)
df   <- read_tidy(opt$data, opt$sheet)

fn <- paste0("recipe_", spec$recipe)
if (!exists(fn, mode = "function")) stop(sprintf("unknown recipe '%s'", spec$recipe))
res <- match.fun(fn)(df, spec)

rec  <- res$resolved$label %||% res$test_meta$name %||% "(chosen at plot time)"
ts   <- opt$stamp %||% timestamp_utc()
base <- slugify(res$label %||% spec$data$y %||% spec$data$x %||% spec$recipe)
ddir <- file.path(opt$out, sprintf("diagnose_%s_%s_%s", spec$recipe, base, ts))
dir.create(ddir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Recipe:            %s\n", spec$recipe))
cat(sprintf("Recommended test:  %s\n", rec))

if (!is.null(res$qc)) {
  qcres <- qc_normality(res$qc)
  d <- list(width = qcres$width, height = qcres$height)
  png_path <- file.path(ddir, "qc_normality.png")
  ggplot2::ggsave(png_path, qcres$plot, width = d$width, height = d$height, dpi = 200, bg = "white")
  ggplot2::ggsave(file.path(ddir, "qc_normality.pdf"), qcres$plot,
                  width = d$width, height = d$height, device = grDevices::cairo_pdf)

  cat(sprintf("Normality check on: %s\n\n", qcres$quantity))
  pf <- function(x, g = FALSE) if (is.na(x)) "NA" else formatC(x, format = if (g) "g" else "f", digits = 3)
  for (i in seq_len(nrow(qcres$panels))) {
    r <- qcres$panels[i, ]
    cat(sprintf("  %-16s n=%-4d Shapiro p=%-8s skew=%-7s kurtosis=%s\n",
                r$panel, r$n, pf(r$shapiro_p, TRUE), pf(r$skewness), pf(r$kurtosis)))
  }
  cat("\nShapiro-Wilk is a default, not a verdict: it is underpowered at small n\n")
  cat("and over-rejects at large n. Read the QQ plot; mild deviation is usually fine.\n")
  cat(sprintf("\nQQ plot: %s\n", png_path))

  report <- list(
    recipe = spec$recipe, engine = "r", recommended_test = rec,
    normality_quantity = qcres$quantity,
    panels = lapply(seq_len(nrow(qcres$panels)), function(i) {
      r <- qcres$panels[i, ]
      list(panel = r$panel, n = r$n, shapiro_p = r$shapiro_p,
           skewness = r$skewness, kurtosis = r$kurtosis)
    }),
    qq_plot = "qc_normality.png"
  )
  writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = NA),
             file.path(ddir, "diagnostics.json"))
} else {
  cat("\nThis recipe does not choose its test from a normality check, so there is\n")
  cat("no QQ diagnostic. See reference/decision_tree.md for how its test is chosen.\n")
  writeLines(jsonlite::toJSON(list(recipe = spec$recipe, engine = "r",
                                   recommended_test = rec, normality_quantity = NA),
                              auto_unbox = TRUE, pretty = TRUE, na = "null"),
             file.path(ddir, "diagnostics.json"))
}

cat(sprintf("\nWrote diagnostics: %s\n", ddir))
cat("Confirm or override the test, then run `figkit plot` to draw the figure.\n")
