#!/usr/bin/env Rscript
# Generate an external validation cohort for the calibration recipe: 500
# patients with a one-year event (yes or no) and the predicted risks of two
# models. risk_clinical is well calibrated. risk_extended was overfitted in its
# development data, so its predictions are too extreme (calibration slope below
# 1) and too high on average (observed/expected below 1).
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_calibration_data.R
set.seed(4242)

n <- 500
age <- rnorm(n, 0, 1)
marker <- rnorm(n, 0, 1)
lp_true <- -1.1 + 0.8 * age + 0.6 * marker
event <- rbinom(n, 1, plogis(lp_true))
lp_clinical <- lp_true + rnorm(n, 0, 0.35)
lp_extended <- 0.3 + 1.6 * lp_true + rnorm(n, 0, 0.25)
d <- data.frame(patient = sprintf("V%03d", seq_len(n)),
                event_1y = ifelse(event == 1, "yes", "no"),
                risk_clinical = round(plogis(lp_clinical), 4),
                risk_extended = round(plogis(lp_extended), 4))
write.csv(d, "example/calibration.csv", row.names = FALSE)
cat(sprintf("wrote example/calibration.csv (%d patients, %d events)\n", n, sum(event)))
