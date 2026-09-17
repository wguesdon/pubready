#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container.
suppressPackageStartupMessages({ library(ggplot2); library(ggpubr) })

raw <- read.csv("input_value.csv", check.names = FALSE)
raw[["value"]] <- as.numeric(raw[["value"]]); raw <- raw[complete.cases(raw[, c("subject","condition","value")]), ]
a <- raw[raw[["condition"]] == "pre", c("subject","value")]; names(a) <- c("id","v1")
b <- raw[raw[["condition"]] == "post", c("subject","value")]; names(b) <- c("id","v2")
m <- merge(a[!duplicated(a$id), ], b[!duplicated(b$id), ], by = "id")
res <- t.test(m$v2, m$v1, paired = TRUE)
print(res)
pl <- rbind(data.frame(id=m$id, cond="pre", v=m$v1), data.frame(id=m$id, cond="post", v=m$v2))
pl$cond <- factor(pl$cond, levels = c("pre","post"))
p <- ggplot(pl, aes(cond, v)) + geom_line(aes(group = id), colour = "grey70") +
  geom_point(aes(colour = cond), size = 2.4) +
  labs(x = "condition", y = "value") + theme_classic() + theme(legend.position = "none")
ggsave("figure_value_paired.pdf", p, width = 3.8, height = 4.2)
cat("Reproduced figure.\n")
