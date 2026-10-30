#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_enrichment_results.csv", check.names = FALSE)
d <- data.frame(term = as.character(raw[["Description"]]), padj = as.numeric(raw[["p.adjust"]]))
# An adjusted p of exactly 0 is set to the smallest positive double, as in the figure.
d$padj[!is.na(d$padj) & d$padj == 0] <- .Machine$double.xmin
d$count <- as.numeric(raw[["Count"]])
d$x <- sapply(strsplit(as.character(raw[["GeneRatio"]]), "/"), function(a) as.numeric(a[1])/as.numeric(a[2])); xlab <- "gene ratio"
d <- d[order(d$padj), ][seq_len(min(15, nrow(d))), ]
# Wrap each term at about 45 characters, as the figure does.
d$label <- d$term
for (i in seq_len(nrow(d))) d$label[i] <- paste(strwrap(d$term[i], 46), collapse = "\n")
d$label <- factor(d$label, levels = d$label[order(d$x)])
p <- ggplot(d, aes(x, label)) + geom_point(aes(size = count, colour = -log10(padj))) +
  scale_colour_gradient(low = "#3B6DB3", high = "#C1432B", name = "-log10 p.adjust") + labs(x = xlab, y = NULL) + theme_classic()
p <- p
ggsave("figure_enrichment_results_enrich.pdf", p, width = 7.72, height = 7.00)
cat("Reproduced figure.\n")
