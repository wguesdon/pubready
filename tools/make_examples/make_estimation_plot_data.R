#!/usr/bin/env Rscript
# Generate an example for the estimation_plot recipe: the migration speed of
# cells in a wound-healing assay, one value per well, under a vehicle control and
# three compounds. Compound A slows the cells a little, compound B a lot, and
# compound C has no effect. The levels sort with "vehicle" last, so the command
# names the reference with --ref vehicle. Run once inside the pubready container
# from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_estimation_plot_data.R
set.seed(20261029)

design <- data.frame(treatment = c("vehicle", "compound_A", "compound_B", "compound_C"),
                     n = c(14, 12, 12, 10), mu = c(21.0, 15.0, 12.5, 21.6),
                     sd = c(3.8, 4.0, 3.2, 4.6))
rows <- list()
for (i in seq_len(nrow(design))) {
  r <- design[i, ]
  rows[[i]] <- data.frame(treatment = r$treatment,
                          speed_um_h = round(stats::rnorm(r$n, r$mu, r$sd), 1))
}
d <- do.call(rbind, rows)
d <- cbind(well = sprintf("W%02d", seq_len(nrow(d))), d)
write.csv(d, "example/estimation_plot.csv", row.names = FALSE)
cat(sprintf("wrote example/estimation_plot.csv (%d wells)\n", nrow(d)))
