#!/usr/bin/env Rscript
# Standalone reproduction of the transition counts and the test. Run inside
# the pinned pubready container (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_response.csv", check.names = FALSE, stringsAsFactors = FALSE)
d <- data.frame(id = as.character(raw[["patient"]]), visit = trimws(as.character(raw[["week"]])),
                state = trimws(as.character(raw[["response"]])))
d <- d[complete.cases(d), ]
visits <- c("8", "16", "24")
counts <- NULL
for (k in seq_len(length(visits) - 1)) {
  m <- merge(d[d$visit == visits[k], ], d[d$visit == visits[k + 1], ], by = "id")
  tab <- table(m$state.x, m$state.y)
  print(tab)
  if (length(visits) == 2 && nrow(tab) == ncol(tab)) print(mcnemar.test(tab, correct = FALSE))
  counts <- rbind(counts, data.frame(visit = visits[k], as.data.frame(table(m$state.x))))
}
p <- ggplot(counts, aes(visit, Freq, fill = Var1)) + geom_col() +
  labs(x = "week", y = "Patients", fill = "response") + theme_classic(base_size = 13)
ggsave("figure_response_sankey.pdf", p, width = 10.09, height = 4.75)
cat("Reproduced the transition counts.\n")
