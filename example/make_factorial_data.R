#!/usr/bin/env Rscript
# Generate the factorial example datasets deterministically. Run once inside the
# pubplot container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.3.0 \
#     Rscript example/make_factorial_data.R
set.seed(42)

# Two-way: genotype (wt/ko) x treatment (vehicle/low/high).
# ko responds to treatment, wt does not -> a genotype x treatment interaction.
tw <- expand.grid(rep = 1:8,
                  treatment = c("vehicle", "low", "high"),
                  genotype  = c("wt", "ko"),
                  stringsAsFactors = FALSE)
mu_tw <- 20 +
  ifelse(tw$genotype == "ko" & tw$treatment == "low",  6,  0) +
  ifelse(tw$genotype == "ko" & tw$treatment == "high", 12, 0) +
  ifelse(tw$treatment == "low",  1.0, 0) +
  ifelse(tw$treatment == "high", 0.5, 0)
tw$response <- round(mu_tw + rnorm(nrow(tw), 0, 2), 1)
write.csv(tw[, c("genotype", "treatment", "response")],
          "example/twoway_response.csv", row.names = FALSE)

# Three-way: genotype (wt/ko) x treatment (vehicle/drug) x sex (male/female).
th <- expand.grid(rep = 1:6,
                  sex       = c("male", "female"),
                  treatment = c("vehicle", "drug"),
                  genotype  = c("wt", "ko"),
                  stringsAsFactors = FALSE)
mu_th <- 10 +
  ifelse(th$genotype == "ko", 3, 0) +
  ifelse(th$treatment == "drug", 4, 0) +
  ifelse(th$genotype == "ko" & th$treatment == "drug", 3, 0) +
  ifelse(th$sex == "male", 1, 0)
th$response <- round(mu_th + rnorm(nrow(th), 0, 1.5), 1)
write.csv(th[, c("genotype", "treatment", "sex", "response")],
          "example/threeway_response.csv", row.names = FALSE)

cat("wrote example/twoway_response.csv and example/threeway_response.csv\n")
