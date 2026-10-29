#!/usr/bin/env Rscript
# Generate a response dataset for the sankey recipe: 80 patients assessed at
# weeks 8, 16 and 24, each with a response category at each assessment. The
# categories follow a Markov chain: responses deepen over time and progressive
# disease is absorbing. Most patients with progressive disease leave the study,
# and a few others miss an assessment, so the later visits hold fewer patients.
# The rows are sorted by week and by response, so the categories first appear in
# their clinical order, which the recipe keeps.
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_sankey_data.R
set.seed(1618)

states <- c("complete response", "partial response", "stable disease", "progressive disease")
start <- c(0.05, 0.30, 0.45, 0.20)
move <- rbind(c(0.85, 0.10, 0.00, 0.05),
              c(0.20, 0.55, 0.10, 0.15),
              c(0.00, 0.15, 0.50, 0.35),
              c(0.00, 0.00, 0.00, 1.00))
n <- 80
weeks <- c(8, 16, 24)
rows <- NULL
for (i in seq_len(n)) {
  s <- sample(4, 1, prob = start)
  present <- TRUE
  for (w in weeks) {
    if (w > weeks[1]) {
      s <- sample(4, 1, prob = move[s, ])
      leave <- if (s == 4) 0.6 else 0.05
      if (present && stats::runif(1) < leave) present <- FALSE
    }
    if (present || w == weeks[1]) {
      rows <- rbind(rows, data.frame(patient = sprintf("R%02d", i), week = w,
                                     response = states[s], rank = s))
    }
  }
}
rows <- rows[order(rows$week, rows$rank, rows$patient), c("patient", "week", "response")]
write.csv(rows, "example/sankey.csv", row.names = FALSE)
cat(sprintf("wrote example/sankey.csv (%d rows, %d patients)\n", nrow(rows), n))
