#!/usr/bin/env Rscript
# Standalone reproduction of this figure. Run inside the pinned pubready
# container (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(rstatix); library(readr)
})

df <- as.data.frame(read_csv("input_concentration.csv", show_col_types = FALSE))
df[["group"]] <- factor(df[["group"]])
df[["concentration"]] <- as.numeric(df[["concentration"]])

stat <- rstatix::wilcox_test(df, concentration ~ group, paired = FALSE)
stat <- rstatix::add_significance(stat)
stat <- rstatix::add_xy_position(stat, x = "group")

p <- ggplot(df, aes(x = .data[["group"]], y = .data[["concentration"]])) +
  geom_boxplot(aes(fill = .data[["group"]]), width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(aes(fill = .data[["group"]]), width = 0.12, size = 1.6, alpha = 0.75, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("#3B6DB3", "#C1432B")) +
  stat_pvalue_manual(stat, label = "p.signif", tip.length = 0.01, bracket.size = 0.4) +
  labs(x = "Group", y = "IL-6 (pg/mL)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave("figure_concentration_mwu.pdf", p, width = 3.8, height = 4.0)
cat("Reproduced figure.\n")
