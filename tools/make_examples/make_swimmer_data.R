#!/usr/bin/env Rscript
# Generate a swimmer dataset for the swimmer recipe: 24 patients in two arms, one
# row per patient, with the months on treatment, whether the treatment is still
# ongoing at the data cut-off, the month of the first response and the month of
# progression. Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_swimmer_data.R
set.seed(77)

n <- 24
arm <- rep(c("control", "treated"), each = n / 2)
months <- round(ifelse(arm == "treated", rexp(n, 1 / 14), rexp(n, 1 / 7)) + 1, 1)
ongoing <- rbinom(n, 1, ifelse(arm == "treated", 0.4, 0.15))
responder <- rbinom(n, 1, ifelse(arm == "treated", 0.6, 0.25)) == 1
response_month <- ifelse(responder, round(pmin(months, runif(n, 1.5, 4)), 1), NA)
progression_month <- ifelse(ongoing == 0 & rbinom(n, 1, 0.7) == 1, months, NA)
d <- data.frame(patient = sprintf("PT-%02d", seq_len(n)), arm = arm, months = months,
                ended = 1L - ongoing, response_month = response_month,
                progression_month = progression_month)
write.csv(d, "example/swimmer_patients.csv", row.names = FALSE)
cat(sprintf("wrote example/swimmer_patients.csv (%d patients, %d ongoing)\n", n, sum(ongoing)))
