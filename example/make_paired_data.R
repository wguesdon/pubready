#!/usr/bin/env Rscript
# Generate a paired before/after dataset for the paired_compare recipe: each
# subject is measured under two conditions. Run once inside the pubplot container
# from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.0 \
#     Rscript example/make_paired_data.R
set.seed(8)

n <- 22
subject <- sprintf("P%02d", seq_len(n))
baseline <- rnorm(n, mean = 100, sd = 15)          # per-subject level
pre  <- baseline + rnorm(n, 0, 4)
post <- baseline - rnorm(n, 12, 5)                 # treatment lowers the value

long <- data.frame(
  subject   = rep(subject, 2),
  condition = rep(c("pre", "post"), each = n),
  value     = round(c(pre, post), 1),
  stringsAsFactors = FALSE)
write.csv(long, "example/paired_response.csv", row.names = FALSE)
cat(sprintf("wrote example/paired_response.csv (%d subjects, 2 conditions)\n", n))
