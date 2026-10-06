#!/usr/bin/env Rscript
# Generate a tumour growth dataset for the timecourse recipe: three arms of eight
# mice, with the tumour volume measured on days 0 to 21. The vehicle tumours grow
# fastest. One vehicle mouse is removed after day 17, so the file carries the
# missing visit that a study always has. Run once inside the pubready container
# from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.5.0 \
#     Rscript tools/make_examples/make_timecourse_data.R
set.seed(13)

days <- c(0, 3, 7, 10, 14, 17, 21)
arms <- data.frame(arm = c("vehicle", "drug_low", "drug_high"),
                   rate = c(0.11, 0.08, 0.05), stringsAsFactors = FALSE)
rows <- list()
for (i in seq_len(nrow(arms))) {
  for (m in 1:8) {
    start <- rnorm(1, 100, 12)
    rate <- arms$rate[i] * exp(rnorm(1, 0, 0.12))
    vol <- start * exp(rate * days) * exp(rnorm(length(days), 0, 0.06))
    rows[[length(rows) + 1]] <- data.frame(mouse = paste0(arms$arm[i], "_", m), arm = arms$arm[i],
                                           day = days, volume_mm3 = round(vol, 1))
  }
}
d <- do.call(rbind, rows)
d <- d[!(d$mouse == "vehicle_3" & d$day == 21), ]
write.csv(d, "example/tumor_growth.csv", row.names = FALSE)
cat(sprintf("wrote example/tumor_growth.csv (%d rows, 24 mice)\n", nrow(d)))
