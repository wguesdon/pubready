#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })

raw <- read.csv("input_gene_b_vs_gene_a.csv", check.names = FALSE)
x <- suppressWarnings(as.numeric(raw[["gene_a"]]))
y <- suppressWarnings(as.numeric(raw[["gene_b"]]))
keep <- is.finite(x) & is.finite(y)
df <- data.frame(x = x[keep], y = y[keep])
ct <- cor.test(df$x, df$y, method = "pearson")
print(ct)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(colour = "#3B6DB3", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#333333", fill = "grey75") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  labs(x = "Gene A (a.u.)", y = "Gene B (a.u.)") + theme_classic()
ggsave("figure_gene_b_vs_gene_a_corr.pdf", p, width = 4.6, height = 4.2)
cat("Reproduced figure.\n")
