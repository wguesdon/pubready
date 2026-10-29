#!/usr/bin/env Rscript
# Generate the committed ma_plot example: a synthetic differential-expression
# table shaped like a DESeq2 result, with the baseMean column that an MA plot
# needs. The standard error of a fold change falls with the mean count, so the
# low-count genes spread wide and seldom reach significance, as in real data.
# Genes with a mean count below 5 get no adjusted p, as DESeq2's independent
# filtering does. Deterministic (seeded) and base R only. Run once inside the
# pubready container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:1.0.0 \
#     Rscript tools/make_examples/make_de_basemean_data.R
set.seed(20261029)

n_genes <- 6000
base_mean <- 10^stats::rnorm(n_genes, mean = 2.4, sd = 0.75)
se <- 0.08 + 1.6 / sqrt(base_mean)

# 7% of the genes are truly differential, the rest have a null effect.
is_de <- stats::runif(n_genes) < 0.07
true_lfc <- ifelse(is_de, sample(c(-1, 1), n_genes, replace = TRUE) * (0.8 + stats::rexp(n_genes, 0.8)), 0)
lfc <- true_lfc + stats::rnorm(n_genes, 0, se)
stat <- lfc / se
pvalue <- 2 * stats::pnorm(-abs(stat))
padj <- rep(NA_real_, n_genes)
tested <- base_mean >= 5
padj[tested] <- stats::p.adjust(pvalue[tested], method = "BH")

out <- data.frame(
  gene = sprintf("GENE%04d", seq_len(n_genes)),
  baseMean = signif(base_mean, 6),
  log2FoldChange = signif(lfc, 5),
  lfcSE = signif(se, 4),
  stat = signif(stat, 5),
  pvalue = signif(pvalue, 4),
  padj = signif(padj, 4),
  stringsAsFactors = FALSE)
out <- out[order(out$padj, out$gene, na.last = TRUE), ]

utils::write.csv(out, "example/de_basemean.csv", row.names = FALSE, na = "NA")
n_sig <- sum(out$padj <= 0.05 & abs(out$log2FoldChange) >= 1, na.rm = TRUE)
cat(sprintf("wrote example/de_basemean.csv (%d genes, %d at padj <= 0.05 and |log2FC| >= 1)\n",
            nrow(out), n_sig))
