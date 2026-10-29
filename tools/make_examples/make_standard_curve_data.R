#!/usr/bin/env Rscript
# Generate an example for the standard_curve recipe: a sandwich ELISA plate read
# at 450 nm. Eight standards in duplicate (a blank and a two-fold series from
# 0.156 to 10 ng/mL) follow a four-parameter logistic curve. Twenty serum samples
# have no concentration: the recipe reads them off the curve. Two samples sit
# below the lowest standard and one above the highest, so they are flagged and
# not extrapolated. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_standard_curve_data.R
set.seed(20261031)

bottom <- 0.055; top <- 3.2; ec50 <- 2.4; hill <- 1.25
od <- function(conc) bottom + (top - bottom) / (1 + (ec50 / pmax(conc, 1e-12))^hill)
std_conc <- c(0, 10 / 2^(6:0))
std <- data.frame(sample = sprintf("STD%d_%s", rep(seq_along(std_conc), each = 2), c("a", "b")),
                  conc_ng_ml = rep(std_conc, each = 2))
std$absorbance <- round(od(std$conc_ng_ml) * (1 + stats::rnorm(nrow(std), 0, 0.03)) +
                          stats::rnorm(nrow(std), 0, 0.004), 3)
truth <- c(0.05, 0.09, exp(stats::runif(16, log(0.25), log(8))), 4.1, 18)
unk <- data.frame(sample = sprintf("Serum_%02d", seq_along(truth)), conc_ng_ml = NA)
unk$absorbance <- round(od(truth) * (1 + stats::rnorm(length(truth), 0, 0.03)), 3)
d <- rbind(std, unk)
write.csv(d, "example/standard_curve.csv", row.names = FALSE, na = "")
cat(sprintf("wrote example/standard_curve.csv (%d standards, %d unknowns)\n", nrow(std), nrow(unk)))
