#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggpubr); library(readr)
})

df <- as.data.frame(read_csv("input_arm.csv", show_col_types = FALSE))
df[["arm"]] <- factor(df[["arm"]])
fit <- surv_fit(Surv(`time`, `event`) ~ `arm`, data = df)

sd <- survdiff(Surv(`time`, `event`) ~ `arm`, data = df)
p_lr <- pchisq(sd$chisq, length(sd$n) - 1, lower.tail = FALSE)
tmax <- max(df[["time"]])
# The risk table of survminer counts the subjects with a time at or after each break.
g <- ggsurvplot(fit, data = df, pval = FALSE, risk.table = TRUE, conf.int = FALSE,
                censor = TRUE, palette = c("#3B6DB3", "#C1432B"), legend = "right", legend.title = "arm",
                legend.labs = c("control", "treated"), xlab = "Months", ylab = "Survival probability",
                xlim = c(0, tmax), break.time.by = 5, risk.table.y.text = FALSE,
                ggtheme = theme_classic(base_size = 13),
                tables.theme = theme_cleantable())
g$plot <- g$plot + labs(subtitle = if (p_lr < 0.001) "Log-rank P < 0.001" else sprintf("Log-rank P = %.2g", p_lr))
p <- ggarrange(g$plot, g$table, ncol = 1, nrow = 2, heights = c(0.74, 0.26), align = "v")

ggsave("figure_arm_km.pdf", p, width = 5.67, height = 5.8)
cat("Reproduced figure.\n")
