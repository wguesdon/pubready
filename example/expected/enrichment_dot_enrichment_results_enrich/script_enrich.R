#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubplot container.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_enrichment_results.csv", check.names = FALSE)
d <- data.frame(term = as.character(raw[["Description"]]), padj = as.numeric(raw[["p.adjust"]]))
d$count <- as.numeric(raw[["Count"]])
d$x <- sapply(strsplit(as.character(raw[["GeneRatio"]]), "/"), function(a) as.numeric(a[1])/as.numeric(a[2])); xlab <- "gene ratio"
d <- d[order(d$padj), ][seq_len(min(15, nrow(d))), ]
d$term <- factor(d$term, levels = d$term[order(d$x)])
p <- ggplot(d, aes(x, term)) + geom_point(aes(size = count, colour = padj)) +
  scale_colour_gradient(low = "#C1432B", high = "#3B6DB3") + labs(x = xlab, y = NULL) + theme_classic()
p <- p
ggsave("figure_enrichment_results_enrich.pdf", p, width = 7.2, height = 6)
cat("Reproduced figure.\n")
