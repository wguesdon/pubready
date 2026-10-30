#!/usr/bin/env Rscript
# Standalone reproduction of this figure. Run inside the pinned pubready
# container (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(ggpubr); library(rstatix); library(readr)
})

df <- as.data.frame(read_csv("input_volume.csv", show_col_types = FALSE))
df[["group"]] <- factor(df[["group"]], levels = c("control", "treated"))
df[["volume"]] <- as.numeric(df[["volume"]])
# The test runs on a copy with syntactic column names.
sdf <- data.frame(.x = df[["group"]], .y = df[["volume"]])

stat <- rstatix::t_test(sdf, .y ~ .x, paired = FALSE, var.equal = FALSE)
stat <- rstatix::add_significance(stat)
stat <- rstatix::add_xy_position(stat, x = ".x")
# The bracket sits 6% of the data span above the largest value.
stat$y.position <- max(sdf$.y) + 0.06 * diff(range(sdf$.y))

point_k <- ave(seq_len(nrow(df)), df[["group"]], FUN = seq_along)
df$.point_x <- as.integer(df[["group"]]) + -0.12 + 0.24 * ((point_k * 0.6180339887498949) %% 1)
p <- ggplot(df, aes(x = .data[["group"]], y = .data[["volume"]])) +
  geom_boxplot(aes(fill = .data[["group"]]), width = 0.6, outlier.shape = NA, alpha = 0.9) +
  geom_point(aes(x = .point_x, fill = .data[["group"]]), size = 1.5, alpha = 0.7, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("#3B6DB3", "#C1432B")) +
  stat_pvalue_manual(stat, label = "p.signif", tip.length = 0.01, bracket.size = 0.4) +
  labs(x = "Group", y = "Tumor volume (mm^3)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

ggsave("figure_volume_prism.pdf", p, width = 3.8, height = 4.0)
cat("Reproduced figure.\n")
