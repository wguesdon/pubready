#!/usr/bin/env Rscript
# Generate the UpSet example data deterministically: a binary set-membership
# matrix. First column = element id; each remaining column is a set, with 1 when
# the element belongs to the set and 0 otherwise. Run once inside the pubplot
# container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.0 \
#     Rscript tools/make_examples/make_upset_data.R
set.seed(7)

n <- 80
genes <- sprintf("gene_%03d", 1:n)
membership <- data.frame(
  gene         = genes,
  apoptosis    = rbinom(n, 1, 0.45),
  inflammation = rbinom(n, 1, 0.40),
  metabolism   = rbinom(n, 1, 0.35),
  cell_cycle   = rbinom(n, 1, 0.30)
)

write.csv(membership, "example/set_membership.csv", row.names = FALSE)
cat(sprintf("wrote example/set_membership.csv (%d elements x %d sets)\n", n, 4))
