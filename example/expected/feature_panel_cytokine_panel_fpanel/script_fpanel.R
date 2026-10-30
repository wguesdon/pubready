#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_cytokine_panel.csv", check.names = FALSE)
d <- data.frame(feature = as.character(raw[["cytokine"]]), group = as.character(raw[["group"]]), value = as.numeric(raw[["concentration_pg_ml"]]))
d <- d[is.finite(d$value) & !is.na(d$group), ]
d$group <- factor(d$group, levels = c("Healthy control", "Mild disease", "Severe disease"))
feats <- unique(d$feature)
test_one <- function(f) kruskal.test(value ~ group, data = d[d$feature == f, ])$p.value
p <- vapply(feats, test_one, numeric(1))
q <- p.adjust(p, method = "BH")
titles <- paste0(feats, "\nq = ", formatC(q, format = "g", digits = 2))
d$panel <- factor(titles[match(d$feature, feats)], levels = titles)
g <- ggplot(d, aes(group, value, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(shape = 21, width = 0.12, height = 0) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  theme_classic(base_size = 13) + theme(legend.position = "none", strip.background = element_blank())
ggsave("figure_cytokine_panel_fpanel.pdf", g, width = 9, height = 8)
cat("Reproduced figure.\n")
