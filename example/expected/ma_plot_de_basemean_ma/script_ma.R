#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

raw <- read.csv("input_de_basemean.csv", check.names = FALSE)
m <- as.numeric(raw[["log2FoldChange"]]); p <- as.numeric(raw[["padj"]])
a <- log10(as.numeric(raw[["baseMean"]])); lab <- as.character(raw[["gene"]])
keep <- is.finite(m) & is.finite(a); m <- m[keep]; p <- p[keep]; a <- a[keep]; lab <- lab[keep]
fcc <- 1; pcc <- 0.05
sig <- !is.na(p) & p <= pcc
cat <- factor(ifelse(sig & m >= fcc, "Up", ifelse(sig & m <= -fcc, "Down", "NS")), levels = c("NS", "Down", "Up"))
d <- data.frame(A = a, M = m, cat = cat, lab = lab)
hits <- d[d$cat != "NS", ]; hits <- hits[order(p[d$cat != "NS"], -abs(hits$M)), ]
top <- head(hits, 15)
lo <- min(m, -1.5 * fcc); hi <- max(m, 1.5 * fcc); lim <- c(lo, hi) + c(-0.06, 0.06) * (hi - lo)
g <- ggplot(d, aes(A, M, colour = cat)) +
  geom_point(data = d[d$cat == "NS", ], size = 0.45, shape = 16) +
  geom_hline(yintercept = 0) + geom_hline(yintercept = c(-fcc, fcc), linetype = "dashed", colour = "grey40") +
  geom_point(data = hits, size = 1.2, alpha = 0.85, shape = 16) +
  geom_text_repel(data = top, aes(label = lab), colour = "black", size = 3, seed = 1, max.overlaps = Inf) +
  scale_colour_manual(values = c(NS = "grey70", Down = "#3B6DB3", Up = "#C1432B"), drop = FALSE, name = NULL) +
  scale_y_continuous(limits = lim) + labs(x = "log10 mean of normalized counts", y = "log2 fold change") +
  theme_classic(base_size = 13) + theme(legend.position = "bottom")
ggsave("figure_de_basemean_ma.pdf", g, width = 6, height = 5.4)
cat("Reproduced figure.\n")
