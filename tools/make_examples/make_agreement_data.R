#!/usr/bin/env Rscript
# Generate a method comparison dataset for the bland_altman recipe: 40 samples
# measured with a reference assay and a new assay. The new assay reads about 3%
# high plus a constant 2 units, so the plot shows a small proportional bias. Run
# once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_agreement_data.R
set.seed(55)

n <- 40
truth <- runif(n, 40, 220)
reference_assay <- truth + rnorm(n, 0, 3)
new_assay <- 2 + 1.03 * truth + rnorm(n, 0, 3)
d <- data.frame(sample = sprintf("X%02d", seq_len(n)),
                reference_assay = round(reference_assay, 1),
                new_assay = round(new_assay, 1))
write.csv(d, "example/method_agreement.csv", row.names = FALSE)
cat(sprintf("wrote example/method_agreement.csv (%d samples)\n", n))
