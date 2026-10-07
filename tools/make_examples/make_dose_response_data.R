#!/usr/bin/env Rscript
# Generate a cell viability dataset for the dose_response recipe: three compounds,
# eight 3-fold concentrations and three replicate wells each. Viability falls with
# the dose, so each curve reports an IC50. Run once inside the pubready container
# from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_dose_response_data.R
set.seed(21)

conc <- c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10)
compounds <- data.frame(compound = c("cmpd_A", "cmpd_B", "cmpd_C"),
                        ic50 = c(0.05, 0.3, 2.0), hill = c(-1.2, -1.0, -0.8),
                        stringsAsFactors = FALSE)
rows <- list()
for (i in seq_len(nrow(compounds))) {
  cp <- compounds[i, ]
  mu <- 5 + (100 - 5) / (1 + 10^((log10(cp$ic50) - log10(conc)) * cp$hill))
  for (rep in 1:3) {
    rows[[length(rows) + 1]] <- data.frame(
      compound = cp$compound, concentration_um = conc, replicate = rep,
      viability_pct = round(mu + rnorm(length(conc), 0, 4), 1))
  }
}
d <- do.call(rbind, rows)
d <- d[order(d$compound, d$concentration_um, d$replicate), ]
write.csv(d, "example/dose_response.csv", row.names = FALSE)
cat(sprintf("wrote example/dose_response.csv (%d rows, 3 compounds)\n", nrow(d)))
