#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(rstatix); library(readr); library(ARTool)
})

df <- as.data.frame(read_csv("input_response.csv", show_col_types = FALSE))
df[["genotype"]] <- factor(df[["genotype"]]); df[["treatment"]] <- factor(df[["treatment"]]); df[["sex"]] <- factor(df[["sex"]]); 
df[["response"]] <- as.numeric(df[["response"]])

at <- as.data.frame(rstatix::anova_test(df, response ~ `genotype` * `treatment` * `sex`, type = 2))
eff <- data.frame(effect = at$Effect, p = at$p)
plab <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))
sig <- eff[which(eff$p < 0.05), , drop = FALSE]
subtitle_txt <- if (nrow(sig) > 0) paste(sprintf("%s: %s", gsub(":", " × ", sig$effect), vapply(sig$p, plab, character(1))), collapse = "\n") else "No significant effects (P ≥ 0.05)"

p <- ggplot(df, aes(x = .data[["genotype"]], y = .data[["response"]], fill = .data[["treatment"]])) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.3, alpha = 0.7, shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("#3B6DB3", "#C1432B")) +
  labs(x = "genotype", y = "Response (a.u.)", fill = "treatment", subtitle = subtitle_txt) +
  theme_classic(base_size = 13) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey30"))
p <- p + facet_wrap(vars(.data[["sex"]]))

ggsave("figure_response_threeway.pdf", p, width = 6.8, height = 4.4)
cat("Reproduced figure.\n")
