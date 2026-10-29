#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

raw <- read.csv("input_cell_type.csv", check.names = FALSE)
d <- data.frame(id = as.character(raw[["sample"]]), category = as.character(raw[["cell_type"]]))
d$count <- as.numeric(raw[["n_cells"]])
d$group <- as.character(raw[["condition"]])
d <- d[stats::complete.cases(d), ]
cats <- unique(d$category)
tab <- xtabs(count ~ id + factor(category, levels = cats), data = d)
P <- prop.table(tab, 1)
P <- P[rowSums(tab) > 0, , drop = FALSE]
g <- d$group[match(rownames(P), d$id)]
lv <- sort(unique(g), method = "radix")
if (length(lv) >= 2) for (j in seq_along(cats)) {
  v <- as.numeric(P[, j])
  if (length(unique(v)) < 2) next
  if (length(lv) == 2) print(suppressWarnings(wilcox.test(v[g == lv[1]], v[g == lv[2]])))
  else print(kruskal.test(v, factor(g, levels = lv)))
}

long <- data.frame(sample = rep(rownames(P), times = length(cats)), group = rep(g, times = length(cats)),
                   category = factor(rep(cats, each = nrow(P)), levels = cats), proportion = as.vector(P))
p <- ggplot(long, aes(x = sample, y = proportion, fill = category)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  scale_y_continuous(labels = label_percent(), expand = expansion(mult = c(0, 0.01))) +
  facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
  labs(x = "sample", y = "Proportion", fill = "cell_type") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = element_blank())

ggsave("figure_cell_type_comp.pdf", p, width = max(4.8, 2.5 + 0.3 * nrow(P)), height = 4.8)
cat("Reproduced figure.\n")
