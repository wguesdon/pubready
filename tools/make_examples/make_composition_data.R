#!/usr/bin/env Rscript
# Generate the example dataset of the composition recipe: the cell-type counts
# of a single-cell experiment, 12 blood samples in two conditions. The treated
# samples hold more CD8 T cells and fewer B cells. Each sample draws its
# proportions from a Dirichlet distribution around the condition mean and its
# counts from a multinomial distribution. Seeded, so a rerun writes the same
# file. Run inside the pubready container from the repository root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_composition_data.R
set.seed(20261030)

cell_types <- c("CD4 T", "CD8 T", "B", "NK", "classical monocyte",
                "non-classical monocyte", "cDC", "pDC")
mean_control <- c(0.32, 0.14, 0.13, 0.10, 0.20, 0.05, 0.04, 0.02)
mean_treated <- c(0.30, 0.24, 0.06, 0.10, 0.19, 0.05, 0.04, 0.02)

# One sample: Dirichlet proportions (normalised gamma draws), multinomial counts.
draw_sample <- function(id, condition, mu) {
  a <- stats::rgamma(length(mu), shape = 60 * mu)
  n_cells <- sample(2500:6000, 1)
  counts <- as.vector(stats::rmultinom(1, n_cells, a / sum(a)))
  data.frame(sample = id, condition = condition, cell_type = cell_types, n_cells = counts,
             stringsAsFactors = FALSE)
}

df <- do.call(rbind, c(
  lapply(sprintf("ctrl_%02d", 1:6), draw_sample, condition = "control", mu = mean_control),
  lapply(sprintf("trt_%02d", 1:6), draw_sample, condition = "treated", mu = mean_treated)))

utils::write.csv(df, "example/composition.csv", row.names = FALSE)
cat(sprintf("wrote example/composition.csv (%d rows, %d samples)\n", nrow(df), length(unique(df$sample))))
