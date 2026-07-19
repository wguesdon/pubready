#!/usr/bin/env Rscript
# Build the volcano example from the airway benchmark (Himes et al. 2014,
# GSE52778): primary human airway smooth muscle cells, dexamethasone-treated vs
# untreated, four donors. This is the standard Bioconductor DE reference. We run
# DESeq2 with the donor (cell line) blocked, contrast treated vs untreated, and
# label genes by symbol. Run once inside the pubplot container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.1 \
#     Rscript example/make_volcano_data.R
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

data("airway")
airway$dex <- relevel(airway$dex, "untrt")            # untreated as the reference

dds <- DESeqDataSet(airway, design = ~ cell + dex)    # block on donor (cell line)
dds <- dds[rowSums(counts(dds)) > 1, ]                # drop all-zero genes
dds <- DESeq(dds)
res <- results(dds, contrast = c("dex", "trt", "untrt"))
res <- as.data.frame(res)
res$ensembl <- rownames(res)

# Map Ensembl gene IDs to HGNC symbols; fall back to the Ensembl ID when unmapped.
res$gene <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = res$ensembl,
                                  column = "SYMBOL", keytype = "ENSEMBL",
                                  multiVals = "first")
res$gene[is.na(res$gene)] <- res$ensembl[is.na(res$gene)]

out <- data.frame(
  gene           = res$gene,
  log2FoldChange = round(res$log2FoldChange, 4),
  pvalue         = signif(res$pvalue, 4),
  padj           = signif(res$padj, 4),
  stringsAsFactors = FALSE)
out <- out[!is.na(out$pvalue), ]                      # keep genes that were tested
out <- out[order(out$padj), ]

write.csv(out, "example/de_results.csv", row.names = FALSE)
n_sig <- sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)
cat(sprintf("wrote example/de_results.csv (%d genes tested, %d significant at padj<0.05 & |log2FC|>1)\n",
            nrow(out), n_sig))
