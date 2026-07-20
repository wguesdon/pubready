#!/usr/bin/env Rscript
# Generate a samples-by-features table with a group structure for the pca recipe:
# rows are samples, numeric columns are features, plus a group column. Run once
# inside the pubplot container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.0 \
#     Rscript tools/make_examples/make_pca_data.R
set.seed(31)

groups   <- rep(c("healthy", "mild", "severe"), each = 12)
n        <- length(groups)
n_feat   <- 40
# Two latent axes shift with group; the rest is noise, so groups separate on PCA.
axis1 <- c(rep(-2, 12), rep(0, 12), rep(2.5, 12)) + rnorm(n, 0, 0.5)
axis2 <- c(rep(1, 12), rep(-1.5, 12), rep(0.5, 12)) + rnorm(n, 0, 0.5)

load1 <- rnorm(n_feat)
load2 <- rnorm(n_feat)
mat <- outer(axis1, load1) + outer(axis2, load2) + matrix(rnorm(n * n_feat, 0, 0.8), n, n_feat)
colnames(mat) <- sprintf("feature%02d", seq_len(n_feat))
mat <- round(mat, 3)

df <- data.frame(sample = sprintf("S%02d", seq_len(n)), group = groups, mat,
                 check.names = FALSE, stringsAsFactors = FALSE)
write.csv(df, "example/pca_samples.csv", row.names = FALSE)
cat(sprintf("wrote example/pca_samples.csv (%d samples x %d features, 3 groups)\n", n, n_feat))
