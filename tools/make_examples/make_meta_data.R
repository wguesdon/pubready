#!/usr/bin/env Rscript
# Generate a meta-analysis dataset for the forest_meta recipe: nine trials that
# report a hazard ratio with its 95% CI, as a paper reports them. The true log
# hazard ratio is -0.25 with a between-trial SD of 0.12. Run once inside the
# pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.5.0 \
#     Rscript tools/make_examples/make_meta_data.R
set.seed(89)

k <- 9
events <- round(runif(k, 40, 400))
se <- 2 / sqrt(events)
theta <- rnorm(k, -0.25, 0.12)
yi <- rnorm(k, theta, se)
z <- qnorm(0.975)
d <- data.frame(trial = sprintf("Trial %s", LETTERS[seq_len(k)]),
                year = sort(sample(2004:2024, k)),
                events = events,
                hazard_ratio = round(exp(yi), 3),
                ci_low = round(exp(yi - z * se), 3),
                ci_high = round(exp(yi + z * se), 3))
write.csv(d, "example/meta_trials.csv", row.names = FALSE)
cat(sprintf("wrote example/meta_trials.csv (%d trials)\n", k))
