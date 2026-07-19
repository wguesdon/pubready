#!/usr/bin/env Rscript
# Generate a DESeq2-style differential-expression table for the volcano recipe.
# Run once inside the pubplot container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.0 \
#     Rscript example/make_volcano_data.R
set.seed(42)

n_null <- 1800
n_up   <- 120
n_down <- 100
n <- n_null + n_up + n_down
gene <- sprintf("GENE%04d", seq_len(n))

log2fc <- c(rnorm(n_null, 0, 0.4),
            rnorm(n_up,  2.0, 0.6),
            rnorm(n_down, -2.0, 0.6))
# p small for the DE genes, roughly uniform for the null set.
pvalue <- c(runif(n_null, 0, 1),
            10^(-runif(n_up,  3, 12)),
            10^(-runif(n_down, 3, 12)))
padj <- p.adjust(pvalue, method = "BH")

de <- data.frame(gene = gene,
                 log2FoldChange = round(log2fc, 3),
                 pvalue = signif(pvalue, 4),
                 padj = signif(padj, 4))
write.csv(de, "example/de_results.csv", row.names = FALSE)
cat(sprintf("wrote example/de_results.csv (%d genes: %d up, %d down)\n", n, n_up, n_down))
