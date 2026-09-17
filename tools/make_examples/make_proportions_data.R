#!/usr/bin/env Rscript
# Generate a raw contingency dataset (one row per patient) for the proportions
# recipe: treatment arm versus clinical response. Run once inside the pubready
# container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.4.0 \
#     Rscript tools/make_examples/make_proportions_data.R
set.seed(5)

make_arm <- function(arm, n, p_resp) {
  data.frame(arm = arm,
             response = ifelse(runif(n) < p_resp, "responder", "non-responder"),
             stringsAsFactors = FALSE)
}
df <- rbind(
  make_arm("control", 60, 0.30),
  make_arm("treated", 64, 0.62))
df <- df[sample(nrow(df)), ]

write.csv(df, "example/response_by_arm.csv", row.names = FALSE)
cat(sprintf("wrote example/response_by_arm.csv (%d patients)\n", nrow(df)))
