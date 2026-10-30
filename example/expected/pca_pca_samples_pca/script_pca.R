#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_pca_samples.csv", check.names = FALSE)
# The features of the figure; a cell that does not read as a number is NA,
# and a sample with a missing feature or group is dropped, as in the figure.
feats <- c("feature01", "feature02", "feature03", "feature04", "feature05", "feature06", "feature07", "feature08", "feature09", "feature10", "feature11", "feature12", "feature13", "feature14", "feature15", "feature16", "feature17", "feature18", "feature19", "feature20", "feature21", "feature22", "feature23", "feature24", "feature25", "feature26", "feature27", "feature28", "feature29", "feature30", "feature31", "feature32", "feature33", "feature34", "feature35", "feature36", "feature37", "feature38", "feature39", "feature40")
num <- raw[feats]
for (f in feats) num[[f]] <- suppressWarnings(as.numeric(as.character(num[[f]])))
g <- as.character(raw[["group"]])
keep <- complete.cases(num) & !is.na(g) & trimws(g) != ""
num <- num[keep, , drop = FALSE]; grp <- factor(g[keep])
pc <- prcomp(num, center = TRUE, scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2) * 100
d <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], group = grp)
p <- ggplot(d, aes(PC1, PC2, colour = group)) + geom_point(size = 2.4) +
  stat_ellipse(type = "norm", level = 0.95) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic()
ggsave("figure_pca_samples_pca.pdf", p, width = 5.2, height = 4.6)
cat("Reproduced figure.\n")
