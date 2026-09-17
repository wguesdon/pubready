#!/usr/bin/env Rscript
# OPTIONAL: build a real-data volcano input from the airway benchmark (Himes et
# al. 2014, GSE52778): primary human airway smooth muscle cells, dexamethasone-
# treated vs untreated, four donors. The standard Bioconductor DE reference. Run
# DESeq2 with the donor (cell line) blocked, contrast treated vs untreated, and
# label genes by symbol.
#
# The committed volcano example uses the small synthetic table from
# make_volcano_data.R instead, so the repo stays light. This script is here for
# anyone who wants the genuine benchmark; it writes de_results_airway.csv, which
# is gitignored. Point the recipe at it with:
#   figkit plot --recipe volcano --data example/de_results_airway.csv --y padj
#
# Needs DESeq2 + airway + org.Hs.eg.db (already in the image). Run from the repo
# root inside the container:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubready:0.4.1 \
#     Rscript tools/make_examples/make_volcano_airway.R
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

write.csv(out, "example/de_results_airway.csv", row.names = FALSE)
n_sig <- sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)
cat(sprintf("wrote example/de_results_airway.csv (%d genes tested, %d significant at padj<0.05 & |log2FC|>1)\n",
            nrow(out), n_sig))
