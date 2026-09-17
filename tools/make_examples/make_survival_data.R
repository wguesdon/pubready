#!/usr/bin/env Rscript
# Generate a survival example dataset deterministically. Run once inside the
# pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.4.0 \
#     Rscript tools/make_examples/make_survival_data.R
set.seed(7)

make_arm <- function(rate, arm, n) {
  t_event <- rexp(n, rate)
  t_cens  <- runif(n, 0, 24)          # administrative censoring at 24 months
  data.frame(
    arm   = arm,
    time  = round(pmin(t_event, t_cens), 1),
    event = as.integer(t_event <= t_cens),   # 1 = event, 0 = censored
    stringsAsFactors = FALSE
  )
}

# Treated arm has a lower event rate -> better survival (log-rank significant).
df <- rbind(make_arm(0.12, "control", 60),
            make_arm(0.05, "treated", 60))

# Covariates for the Cox model.
df$age   <- round(rnorm(nrow(df), 60, 10))
df$sex   <- sample(c("male", "female"), nrow(df), replace = TRUE)
df$stage <- sample(c("I", "II", "III"), nrow(df), replace = TRUE, prob = c(0.4, 0.35, 0.25))

write.csv(df[, c("time", "event", "arm", "age", "sex", "stage")],
          "example/survival_trial.csv", row.names = FALSE)

cat(sprintf("wrote example/survival_trial.csv (%d rows, %d events)\n",
            nrow(df), sum(df$event)))
