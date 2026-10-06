#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_target_lesion_mm.csv", check.names = FALSE)
raw <- raw[complete.cases(raw[, c("patient", "week", "target_lesion_mm")]), ]
raw <- raw[order(raw[["patient"]], raw[["week"]]), ]
base <- tapply(raw[["target_lesion_mm"]], raw[["patient"]], function(v) v[1])
raw$change <- 100 * (raw[["target_lesion_mm"]] - base[raw[["patient"]]]) / base[raw[["patient"]]]
post <- raw[duplicated(raw[["patient"]]), ]
best <- tapply(post$change, post[["patient"]], min)
d <- data.frame(patient = names(best), best = as.numeric(best))
d <- d[order(-d$best), ]; d$rank <- seq_len(nrow(d))
p <- ggplot(d, aes(rank, best)) + geom_col() + geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(20, -30), linetype = "dashed") +
  labs(x = "Patient", y = "Best change from baseline (%)") + theme_classic()
ggsave("figure_target_lesion_mm_waterfall.pdf", p, width = 6, height = 4.4)
cat("Reproduced figure.\n")
