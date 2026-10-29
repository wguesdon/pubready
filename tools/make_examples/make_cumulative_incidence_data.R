#!/usr/bin/env Rscript
# Generate a competing-risks dataset for the cumulative_incidence recipe: 200
# patients after a stem cell transplant, randomised to a standard or an
# intensified conditioning regimen. Each patient is followed until relapse
# (status 1), death without relapse (status 2) or censoring (status 0). The
# intensified regimen lowers the relapse hazard and raises the non-relapse
# mortality, which is the pattern where Kaplan-Meier overestimates each risk.
# Times are rounded to 0.1 month, so some event times are tied.
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_cumulative_incidence_data.R
set.seed(2718)

n_arm <- 100
arm <- rep(c("standard", "intensified"), each = n_arm)
h_relapse <- ifelse(arm == "standard", 0.024, 0.012)
h_death <- ifelse(arm == "standard", 0.006, 0.013)
t_relapse <- rexp(2 * n_arm, h_relapse)
t_death <- rexp(2 * n_arm, h_death)
# Staggered entry gives an administrative censoring time between 18 and 72
# months, and a small hazard of loss to follow-up.
t_cens <- pmin(runif(2 * n_arm, 18, 72), rexp(2 * n_arm, 0.004))
months <- pmin(t_relapse, t_death, t_cens)
status <- ifelse(months == t_cens, 0L, ifelse(months == t_relapse, 1L, 2L))
months <- pmax(0.1, round(months, 1))

d <- data.frame(patient = sprintf("P%03d", seq_len(2 * n_arm)), arm = arm,
                months = months, status = status)
write.csv(d, "example/cumulative_incidence.csv", row.names = FALSE)
cat(sprintf("wrote example/cumulative_incidence.csv (%d patients: %d relapse, %d death, %d censored)\n",
            nrow(d), sum(status == 1), sum(status == 2), sum(status == 0)))
