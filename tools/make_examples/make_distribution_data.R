#!/usr/bin/env Rscript
# Generate the example dataset of the distribution recipe: the area of single
# cells in three culture conditions. The control cells are close to normal, the
# serum-starved cells are right-skewed (log-normal), and the mixed culture is
# bimodal (two populations of different size). Seeded, so a rerun writes the
# same file. Run inside the pubready container from the repository root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_distribution_data.R
set.seed(20261029)

control <- stats::rnorm(90, mean = 260, sd = 38)
starved <- stats::rlnorm(80, meanlog = log(215), sdlog = 0.3)
mixed_large <- stats::runif(70) < 0.45
mixed <- ifelse(mixed_large, stats::rnorm(70, 380, 30), stats::rnorm(70, 220, 28))

df <- data.frame(
  condition = rep(c("control", "serum-starved", "mixed culture"), times = c(90, 80, 70)),
  cell_area_um2 = round(c(control, starved, mixed), 1),
  stringsAsFactors = FALSE)
df <- df[sample(nrow(df)), ]

utils::write.csv(df, "example/distribution.csv", row.names = FALSE)
cat(sprintf("wrote example/distribution.csv (%d rows)\n", nrow(df)))
