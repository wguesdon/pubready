#!/usr/bin/env Rscript
# Standalone reproduction of this figure. Run inside the pinned pubready
# container (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(rstatix); library(readr)
})

df <- as.data.frame(read_csv("input_volume.csv", show_col_types = FALSE))
df[["group"]] <- factor(df[["group"]])
df[["volume"]] <- as.numeric(df[["volume"]])

stat <- rstatix::t_test(df, volume ~ group, paired = FALSE, var.equal = TRUE)
stat <- rstatix::add_significance(stat)
stat <- rstatix::add_xy_position(stat, x = "group")
stat$plabel <- paste0("p = ", formatC(stat$p, format = "g", digits = 2))

p <- ggplot(df, aes(x = .data[["group"]], y = .data[["volume"]])) +
  geom_violin(aes(fill = .data[["group"]]), trim = FALSE, width = 0.7, alpha = 0.9) +
  geom_jitter(aes(fill = .data[["group"]]), width = 0.12, size = 1.6, alpha = 0.75, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("#2E8B57", "#7A5195")) +
  stat_pvalue_manual(stat, label = "plabel", tip.length = 0.01, bracket.size = 0.4) +
  labs(x = "Group", y = "Tumor volume (mm^3)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave("figure_volume_violin.pdf", p, width = 3.8, height = 4.0)
cat("Reproduced figure.\n")
