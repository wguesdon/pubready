#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container.
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

raw <- read.csv("input_response.csv", check.names = FALSE)
raw <- raw[complete.cases(raw[, c("arm","response")]), ]
tab <- table(raw[["arm"]], raw[["response"]])
print(chisq.test(tab, correct = FALSE))
pd <- as.data.frame(prop.table(tab, 1)); names(pd) <- c("group", "outcome", "prop")
p <- ggplot(pd, aes(group, prop, fill = outcome)) + geom_col(width = 0.7, colour = "white") +
  geom_text(aes(label = percent(prop, accuracy = 1)), position = position_stack(vjust = 0.5), colour = "white") +
  scale_y_continuous(labels = percent) + labs(x = "arm", y = "proportion") + theme_classic()
ggsave("figure_response_prop.pdf", p, width = 4.2, height = 4.2)
cat("Reproduced figure.\n")
