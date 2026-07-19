#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_pca_samples.csv", check.names = FALSE)
num <- raw[, sapply(raw, is.numeric), drop = FALSE]
grp <- factor(raw[["group"]])
num <- num[, apply(num, 2, function(v) sd(v) > 0), drop = FALSE]
pc <- prcomp(num, center = TRUE, scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2) * 100
d <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], group = grp)
p <- ggplot(d, aes(PC1, PC2, colour = group)) + geom_point(size = 2.4) +
  stat_ellipse(level = 0.95) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic()
ggsave("figure_pca_samples_pca.pdf", p, width = 5.2, height = 4.6)
cat("Reproduced figure.\n")
