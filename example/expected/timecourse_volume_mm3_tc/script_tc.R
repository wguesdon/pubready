#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder. The ANOVA uses complete subjects.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_volume_mm3.csv", check.names = FALSE)
d <- data.frame(id = factor(raw[["mouse"]]), time = raw[["day"]], y = as.numeric(raw[["volume_mm3"]]),
                group = factor(raw[["arm"]]))
d <- d[complete.cases(d), ]
n_times <- length(unique(d$time))
d <- d[d$id %in% names(which(table(d$id) == n_times)), ]
d$id <- droplevels(d$id); d$time_f <- factor(d$time)
print(summary(aov(y ~ group * time_f + Error(id / time_f), data = d)))
m <- aggregate(y ~ group + time, data = d, FUN = mean)
p <- ggplot(m, aes(time, y, colour = group)) + geom_line() + geom_point() +
  labs(x = "Day", y = "Tumour volume (mm^3)") + theme_classic()
ggsave("figure_volume_mm3_tc.pdf", p, width = 6.0, height = 4.4)
cat("Reproduced figure.\n")
