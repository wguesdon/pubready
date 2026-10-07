#!/usr/bin/env Rscript
# Generate a diagnostic biomarker dataset for the roc_curve recipe: 60 patients
# with the disease and 80 healthy controls, measured with two correlated serum
# markers. marker_a separates the groups well and marker_b less so. The outcome
# levels sort with "healthy" second, so the example passes --positive disease.
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_roc_data.R
set.seed(34)

n_case <- 60; n_ctrl <- 80
status <- c(rep("disease", n_case), rep("healthy", n_ctrl))
shared <- rnorm(n_case + n_ctrl)
shift_a <- ifelse(status == "disease", 1.4, 0)
shift_b <- ifelse(status == "disease", 0.7, 0)
marker_a <- exp(1.0 + 0.45 * (shift_a + 0.6 * shared + 0.8 * rnorm(n_case + n_ctrl)))
marker_b <- exp(2.0 + 0.35 * (shift_b + 0.6 * shared + 0.8 * rnorm(n_case + n_ctrl)))
d <- data.frame(patient = sprintf("S%03d", seq_len(n_case + n_ctrl)), status = status,
                marker_a = round(marker_a, 2), marker_b = round(marker_b, 2))
write.csv(d, "example/biomarker_roc.csv", row.names = FALSE)
cat(sprintf("wrote example/biomarker_roc.csv (%d cases, %d controls)\n", n_case, n_ctrl))
