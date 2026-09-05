#!/usr/bin/env Rscript
# Generate the spider plot example data deterministically: the target lesion sum
# of 20 patients across five visits, in two arms. One row per measurement, which
# is the shape a trial database exports. Run once inside the pubplot container
# from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.1 \
#     Rscript tools/make_examples/make_spider_data.R
#
# The design puts the response in the treated arm and the growth in the control
# arm, and it stops three patients early, so the file exercises the dropout that
# a trial always has. The values are a sum of target lesion diameters in mm,
# which is what RECIST 1.1 measures.
set.seed(11)

n_patients <- 20
weeks      <- c(0, 8, 16, 24, 32)
arms       <- rep(c("control", "treated"), each = n_patients / 2)
patients   <- sprintf("PT-%02d", 1:n_patients)

# The baseline sum of diameters, between 40 and 90 mm.
baseline <- round(runif(n_patients, 40, 90), 1)

# The weekly trend of each patient, as a fraction of the baseline for each visit.
# The control arm grows and five of its ten patients cross the +20% line. The
# treated arm shrinks and seven of its ten cross the -30% line.
trend <- ifelse(arms == "control",
                runif(n_patients,  0.030,  0.105),
                runif(n_patients, -0.160, -0.015))

# Three patients leave the study after week 16.
last_visit <- rep(5L, n_patients)
last_visit[c(3, 8, 15)] <- 3L

rows <- do.call(rbind, lapply(seq_len(n_patients), function(i) {
  visits <- seq_len(last_visit[i])
  drift  <- rnorm(length(visits), 0, 1.6)
  value  <- baseline[i] * (1 + trend[i] * weeks[visits] / 8) + drift
  value[1] <- baseline[i]
  data.frame(patient = patients[i], arm = arms[i], week = weeks[visits],
             target_lesion_mm = round(pmax(value, 1), 1),
             stringsAsFactors = FALSE)
}))

write.csv(rows, "example/spider_response.csv", row.names = FALSE)

change <- ave(rows$target_lesion_mm, rows$patient,
              FUN = function(v) 100 * (v - v[1]) / v[1])
best <- tapply(change, rows$patient, min)
worst <- tapply(change, rows$patient, max)
cat(sprintf("wrote example/spider_response.csv (%d patients, %d rows)\n",
            n_patients, nrow(rows)))
cat(sprintf("partial response %d, progressive disease %d, stable %d\n",
            sum(best <= -30), sum(best > -30 & worst >= 20),
            sum(best > -30 & worst < 20)))
