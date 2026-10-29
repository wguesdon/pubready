#!/usr/bin/env Rscript
# Generate an example for the enzyme_kinetics recipe: the initial rate of an
# enzyme (nmol/min/mg) at eight substrate concentrations, in triplicate, with no
# inhibitor and with a competitive inhibitor at 2 and 10 uM. A competitive
# inhibitor raises the apparent Km (Km * (1 + I / Ki), Ki = 4 uM) and leaves Vmax
# unchanged, so the shared-Vmax test is not significant and the shared-Km test
# is. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_enzyme_kinetics_data.R
set.seed(20261101)

vmax <- 120; km <- 18; ki <- 4
substrate <- c(2.5, 5, 10, 20, 40, 80, 160, 320)
inhibitor <- c(0, 2, 10)
rows <- list()
for (i in inhibitor) {
  km_app <- km * (1 + i / ki)
  s <- rep(substrate, each = 3)
  v <- vmax * s / (km_app + s)
  rows[[length(rows) + 1]] <- data.frame(
    inhibitor = sprintf("%g uM", i), substrate_um = s,
    rate = round(v * (1 + stats::rnorm(length(s), 0, 0.04)) + stats::rnorm(length(s), 0, 1.2), 2))
}
d <- do.call(rbind, rows)
write.csv(d, "example/enzyme_kinetics.csv", row.names = FALSE)
cat(sprintf("wrote example/enzyme_kinetics.csv (%d rates)\n", nrow(d)))
