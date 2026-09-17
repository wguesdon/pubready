#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(rstatix); library(readr)
})

df <- as.data.frame(read_csv("input_expression.csv", show_col_types = FALSE))
df[["genotype"]] <- factor(df[["genotype"]])
df[["expression"]] <- as.numeric(df[["expression"]])

ph <- rstatix::tukey_hsd(df, expression ~ genotype)
ph$group1 <- as.character(ph$group1); ph$group2 <- as.character(ph$group2)
sig <- ph[which(ph$p.adj < 0.05), , drop = FALSE]
if (nrow(sig) > 0) {
  ymax <- max(df[["expression"]], na.rm = TRUE); rng <- diff(range(df[["expression"]], na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) rng <- max(abs(ymax), 1)
  sig$y.position <- ymax + rng * (0.07 + 0.09 * (seq_len(nrow(sig)) - 1))
  sig$lab <- sig$p.adj.signif
}

omni <- as.data.frame(rstatix::anova_test(df, expression ~ genotype))
subtitle_txt <- paste0("One-way ANOVA, P = ", formatC(omni$p[1], format = "g", digits = 2))

p <- ggplot(df, aes(x = .data[["genotype"]], y = .data[["expression"]])) +
  geom_boxplot(aes(fill = .data[["genotype"]]), width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(fill = .data[["genotype"]]), width = 0.12, size = 1.4, alpha = 0.7, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("#3B6DB3", "#C1432B", "#2E8B57", "#7A5195")) +
  labs(x = "Genotype", y = "Expression (a.u.)", subtitle = subtitle_txt) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey30"))
if (nrow(sig) > 0) p <- p + stat_pvalue_manual(sig, label = "lab", tip.length = 0.01, bracket.size = 0.4)

ggsave("figure_expression_anova.pdf", p, width = 4.2, height = 4.2)
cat("Reproduced figure.\n")
