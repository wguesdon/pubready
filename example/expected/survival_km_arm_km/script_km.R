#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggpubr); library(readr)
})

df <- as.data.frame(read_csv("input_arm.csv", show_col_types = FALSE))
df[["arm"]] <- factor(df[["arm"]])
fit <- surv_fit(Surv(`time`, `event`) ~ `arm`, data = df)

g <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                censor = TRUE, palette = c("#3B6DB3", "#C1432B"), legend.title = "arm",
                legend.labs = c("control", "treated"), xlab = "Months", ylab = "Survival probability",
                ggtheme = theme_classic(base_size = 13),
                tables.theme = theme_cleantable())
p <- ggarrange(g$plot, g$table, ncol = 1, nrow = 2, heights = c(0.74, 0.26), align = "v")

ggsave("figure_arm_km.pdf", p, width = 5.6, height = 5.8)
cat("Reproduced figure.\n")
