#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })

raw <- read.csv("input_gwas_summary.csv", check.names = FALSE)
chr <- sub("^chr", "", as.character(raw[["chr"]]), ignore.case = TRUE)
pos <- as.numeric(raw[["pos"]]); p <- as.numeric(raw[["p"]])
keep <- !is.na(chr) & is.finite(pos) & is.finite(p); chr <- chr[keep]; pos <- pos[keep]; p <- p[keep]
p[p == 0] <- .Machine$double.xmin
num <- grepl("^[0-9]+$", unique(chr)); sp <- match(toupper(unique(chr)), c("X", "Y", "XY", "MT", "M"))
lv <- unique(chr)[order(ifelse(num, 1, ifelse(is.na(sp), 3, 2)), ifelse(num, suppressWarnings(as.numeric(unique(chr))), ifelse(is.na(sp), 0, sp)))]
ci <- match(chr, lv); lo <- tapply(pos, ci, min); span <- tapply(pos, ci, max) - lo
offset <- c(0, cumsum(span + 0.004 * sum(span)))[seq_along(lv)]
d <- data.frame(x = pos - lo[ci] + offset[ci], y = -log10(p), shade = factor(ci %% 2))
lambda <- median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
man <- ggplot(d, aes(x, y, colour = shade)) + geom_point(size = 0.7, shape = 16) +
  scale_colour_manual(values = c("#3B6DB3", "#1F3A5F"), guide = "none") +
  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "#C1432B") +
  scale_x_continuous(breaks = offset + span / 2, labels = lv) +
  labs(x = "Chromosome", y = "-log10 p") + theme_classic(base_size = 13)
n <- length(p); q <- data.frame(exp = -log10((seq_len(n) - 0.5) / n), obs = -log10(sort(p)))
q$lo <- -log10(qbeta(0.975, seq_len(n), n - seq_len(n) + 1)); q$hi <- -log10(qbeta(0.025, seq_len(n), n - seq_len(n) + 1))
qqp <- ggplot(q, aes(exp, obs)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") +
  geom_abline(colour = "#C1432B") + geom_point(size = 0.7, colour = "#3B6DB3") +
  labs(x = "Expected -log10 p", y = "Observed -log10 p", subtitle = sprintf("lambda_GC = %.3f", lambda)) +
  theme_classic(base_size = 13)
g <- ggarrange(man, qqp, ncol = 2, widths = c(2.7, 1), align = "h")
ggsave("figure_gwas_summary_manhattan.pdf", g, width = 10.5, height = 4.3)
cat("Reproduced figure.\n")
