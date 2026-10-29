#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(cluster) })

raw <- read.csv("input_umap_clusters.csv", check.names = FALSE)
d <- data.frame(x = as.numeric(raw[["UMAP_1"]]), y = as.numeric(raw[["UMAP_2"]]), g = raw[["cell_type"]])
d <- d[is.finite(d$x) & is.finite(d$y) & !is.na(d$g), ]
d$g <- factor(d$g)
step <- max(1, ceiling(nrow(d) / 5000)); i <- seq(1, nrow(d), by = step)
s <- silhouette(as.integer(d$g[i]), dist(d[i, c("x", "y")]))
print(summary(s)$clus.avg.widths)
cent <- aggregate(cbind(x, y) ~ g, data = d, FUN = median)
p <- ggplot(d, aes(x, y, colour = g)) + geom_point(size = 0.6, shape = 16) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  geom_label(data = cent, aes(label = g), colour = "black", fill = "white", label.size = 0) +
  coord_equal() + labs(x = "UMAP_1", y = "UMAP_2") + theme_classic(base_size = 13) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave("figure_umap_clusters_embedding.pdf", p, width = 6.6, height = 5.2)
cat("Reproduced figure.\n")
