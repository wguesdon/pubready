#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

raw <- read.csv("input_de_results.csv", check.names = FALSE)
fc <- as.numeric(raw[["log2FoldChange"]]); p <- as.numeric(raw[["pvalue"]])
lab <- as.character(raw[["gene"]])
keep <- is.finite(fc) & is.finite(p) & p > 0; fc <- fc[keep]; p <- p[keep]; lab <- lab[keep]
fcc <- 1; pcc <- 0.05
cat <- ifelse(abs(fc) >= fcc & p <= pcc, "FC and p", ifelse(p <= pcc, "p", ifelse(abs(fc) >= fcc, "FC", "NS")))
cat <- factor(cat, levels = c("NS", "FC", "p", "FC and p"))
pal <- c(NS = "grey70", FC = "#2E8B57", p = "#3B6DB3", `FC and p` = "#C1432B")
d <- data.frame(fc = fc, neglog = -log10(p), cat = cat, lab = lab)
sig <- d[d$cat == "FC and p", ]; sig <- sig[order(-(abs(sig$fc) * sig$neglog)), ]
lab_df <- head(sig, 15)
p <- ggplot(d, aes(fc, neglog, colour = cat)) + geom_point(size = 1.4, alpha = 0.7) +
  scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
  geom_vline(xintercept = c(-fcc, fcc), linetype = "dashed") +
  geom_hline(yintercept = -log10(pcc), linetype = "dashed") +
  geom_text_repel(data = lab_df, aes(label = lab), size = 3, colour = "black", seed = 1) +
  labs(x = "log2 fold change", y = "-log10 p") + theme_classic() + theme(legend.position = "bottom")
ggsave("figure_de_results_volcano.pdf", p, width = 6, height = 5.6)
cat("Reproduced figure.\n")
