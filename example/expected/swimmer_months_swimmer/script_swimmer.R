#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_months.csv", check.names = FALSE)
d <- data.frame(patient = as.character(raw[["patient"]]), duration = as.numeric(raw[["months"]]),
                group = as.character(raw[["arm"]]), ongoing = raw[["ended"]] == 0)
d <- d[order(d$duration), ]; d$rank <- seq_len(nrow(d))
p <- ggplot(d) + geom_segment(aes(x = 0, xend = duration, y = rank, yend = rank, colour = group),
                              linewidth = 3) +
  geom_point(data = d[d$ongoing, ], aes(x = duration, y = rank), shape = 62, size = 4) +
  scale_y_continuous(breaks = d$rank, labels = d$patient) +
  labs(x = "Months on treatment", y = "Patient") + theme_classic()
ggsave("figure_months_swimmer.pdf", p, width = 6.4, height = max(3.6, min(10, 1.6 + 0.18 * nrow(d))))
cat("Reproduced figure.\n")
