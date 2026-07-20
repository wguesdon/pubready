#!/usr/bin/env Rscript
# Generate the heatmap example data deterministically: a features-by-samples
# expression matrix plus a sample annotation table. Run once inside the pubplot
# container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.0 \
#     Rscript tools/make_examples/make_heatmap_data.R
set.seed(11)

n_genes <- 30
n_per   <- 10
samples <- c(sprintf("ctrl_%02d", 1:n_per), sprintf("trt_%02d", 1:n_per))
groups  <- rep(c("control", "treated"), each = n_per)

mat <- matrix(rnorm(n_genes * length(samples), mean = 8, sd = 1.5),
              nrow = n_genes, dimnames = list(sprintf("gene_%02d", 1:n_genes), samples))

trt <- groups == "treated"
mat[1:10, trt]  <- mat[1:10, trt]  + 3   # up in treated
mat[11:20, trt] <- mat[11:20, trt] - 3   # down in treated
# genes 21-30 unchanged

mat <- round(mat, 2)

# Matrix CSV: first column = feature id, remaining columns = samples.
out <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
write.csv(out, "example/expression_matrix.csv", row.names = FALSE)

# Annotation CSV: first column = sample id, then annotation columns.
ann <- data.frame(sample = samples, group = groups)
write.csv(ann, "example/expression_annotation.csv", row.names = FALSE)

cat(sprintf("wrote example/expression_matrix.csv (%d x %d) and annotation\n",
            n_genes, length(samples)))
