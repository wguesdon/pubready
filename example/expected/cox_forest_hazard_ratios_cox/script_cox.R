#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(survival); library(survminer); library(readr)
})

df <- as.data.frame(read_csv("input_hazard_ratios.csv", show_col_types = FALSE))
for (c in c("arm", "age", "sex", "stage")) if (is.character(df[[c]])) df[[c]] <- factor(df[[c]])

cox <- coxph(Surv(`time`, `event`) ~ `arm` + `age` + `sex` + `stage`, data = df)
p <- ggforest(cox, data = df)

ggsave("figure_hazard_ratios_cox.pdf", p, width = 7.4, height = 5.0)
cat("Reproduced figure.\n")
