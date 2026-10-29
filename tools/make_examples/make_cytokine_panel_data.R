#!/usr/bin/env Rscript
# Generate the committed feature_panel example: nine serum cytokines measured
# in 16 subjects of each of three groups, in long form (one row per subject and
# cytokine). Concentrations are log-normal, as cytokine assays are, so the
# Shapiro-Wilk check rejects normality and the recipe picks the Kruskal-Wallis
# test. Six cytokines rise with disease severity and three do not. The group
# labels are long on purpose, to exercise the label wrap. Deterministic (seeded)
# and base R only. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_cytokine_panel_data.R
set.seed(20261029)

groups <- c("Healthy control", "Mild disease", "Severe disease")
n_per <- 16
cytokines <- data.frame(
  name  = c("IL-6", "TNF-alpha", "IL-1beta", "IL-10", "IFN-gamma", "IL-8", "IL-17A", "CXCL10", "MCP-1"),
  base  = c(2.0, 8.0, 1.5, 4.0, 6.0, 12.0, 3.0, 90.0, 150.0),
  shift = c(1.0, 0.6, 0.5, 0.45, 0.0, 0.7, 0.0, 0.9, 0.0))

# The rows of one cytokine: log-normal values whose mean rises with the group.
cytokine_rows <- function(i, cytokines, groups, n_per) {
  g <- rep(seq_along(groups), each = n_per)
  mu <- log(cytokines$base[i]) + cytokines$shift[i] * (g - 1)
  data.frame(subject = sprintf("S%02d", seq_along(g)), group = groups[g],
             cytokine = cytokines$name[i],
             concentration_pg_ml = round(exp(stats::rnorm(length(g), mu, 0.55)), 2),
             stringsAsFactors = FALSE)
}

out <- do.call(rbind, lapply(seq_len(nrow(cytokines)), cytokine_rows,
                             cytokines = cytokines, groups = groups, n_per = n_per))
utils::write.csv(out, "example/cytokine_panel.csv", row.names = FALSE)
cat(sprintf("wrote example/cytokine_panel.csv (%d rows, %d cytokines, %d groups)\n",
            nrow(out), nrow(cytokines), length(groups)))
