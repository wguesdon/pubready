#!/usr/bin/env Rscript
# Generate the committed volcano example: a small, synthetic differential-
# expression table shaped like a DESeq2 result. Deterministic (seeded) and base
# R only, so it is cheap to regenerate and keeps the repo light. The columns
# match what the volcano recipe auto-detects: gene, log2FoldChange, pvalue, padj.
# Run once inside the pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.4.1 \
#     Rscript tools/make_examples/make_volcano_data.R
#
# For the genuine airway benchmark instead, see make_volcano_airway.R.
set.seed(1)

n_genes <- 1200
n_de    <- 150                                   # truly differential (up + down)
n_null  <- n_genes - n_de

# Null genes: fold change scattered around zero, p-values ~ uniform.
null_lfc <- rnorm(n_null, mean = 0, sd = 0.45)
null_p   <- runif(n_null, min = 0, max = 1)

# DE genes: half up, half down, large fold change, very small p-values.
sign_de <- rep(c(1, -1), length.out = n_de)
de_lfc  <- sign_de * (1.2 + rexp(n_de, rate = 0.9))
de_p    <- 10^(-runif(n_de, min = 4, max = 30))

lfc    <- c(null_lfc, de_lfc)
pvalue <- c(null_p, de_p)
padj   <- p.adjust(pvalue, method = "BH")        # BH FDR, as DESeq2 reports

out <- data.frame(
  gene           = sprintf("GENE%04d", seq_len(n_genes)),
  log2FoldChange = round(lfc, 4),
  pvalue         = signif(pvalue, 4),
  padj           = signif(padj, 4),
  stringsAsFactors = FALSE)
out <- out[order(out$padj), ]

write.csv(out, "example/de_results.csv", row.names = FALSE)
n_sig <- sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)
cat(sprintf("wrote example/de_results.csv (%d genes, %d significant at padj<0.05 & |log2FC|>1)\n",
            nrow(out), n_sig))
