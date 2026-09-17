#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_target_lesion_mm.csv", check.names = FALSE)
raw[["week"]] <- as.numeric(raw[["week"]]); raw[["target_lesion_mm"]] <- as.numeric(raw[["target_lesion_mm"]])
raw <- raw[complete.cases(raw[, c("patient","week","target_lesion_mm")]), ]
raw <- raw[order(raw[["patient"]], raw[["week"]]), ]
base <- tapply(raw[["target_lesion_mm"]], raw[["patient"]], function(v) v[1])
raw$change <- 100 * (raw[["target_lesion_mm"]] - base[raw[["patient"]]]) / base[raw[["patient"]]]

p <- ggplot(raw, aes(.data[["week"]], change, group = .data[["patient"]])) +
  geom_hline(yintercept = 0, colour = "grey30") +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = -30, linetype = "dashed", colour = "grey55") +
  geom_line() + geom_point(size = 1.8) +
  labs(x = "week", y = "Change from baseline (%)") + theme_classic()
ggsave("figure_target_lesion_mm_spider.pdf", p, width = 5.8, height = 4.2)
cat("Reproduced figure.\n")
