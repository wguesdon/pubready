#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(readr)
})

df <- as.data.frame(read_csv("input_hazard_ratios.csv", show_col_types = FALSE))
covs <- c("arm", "age", "sex", "stage")
# The first level of each text covariate, in code-point order, is the reference.
for (c in covs) if (is.character(df[[c]])) df[[c]] <- factor(df[[c]], levels = sort(unique(df[[c]]), method = "radix"))

cox <- coxph(Surv(`time`, `event`) ~ `arm` + `age` + `sex` + `stage`, data = df)
ci <- exp(confint(cox)); pv <- summary(cox)$coefficients[, "Pr(>|z|)"]

# One row per model term; the first level of each factor is the reference.
rows <- NULL
for (i in seq_along(cox$assign)) {
  idx <- cox$assign[[i]]; lev <- cox$xlevels[[names(cox$assign)[i]]]
  if (is.null(lev)) {
    rows <- rbind(rows, data.frame(label = covs[i], k = idx[1]))
  } else {
    rows <- rbind(rows, data.frame(label = paste0(covs[i], ": ", lev[seq_len(length(idx) + 1)]),
                                   k = c(NA, idx)))
  }
}
rows$HR <- exp(coef(cox))[rows$k]; rows$lo <- ci[rows$k, 1]; rows$hi <- ci[rows$k, 2]
rows$p <- pv[rows$k]; rows$y <- rev(seq_len(nrow(rows)))
stars <- ifelse(is.na(rows$p), "", ifelse(rows$p < 0.001, "***",
                ifelse(rows$p < 0.01, "**", ifelse(rows$p < 0.05, "*", ""))))
p_txt <- ifelse(rows$p < 0.001, "< 0.001", formatC(rows$p, format = "g", digits = 2))
rows$hr_txt <- ifelse(is.na(rows$k), "reference",
                      sprintf("%.2f (%.2f\u2013%.2f)", rows$HR, rows$lo, rows$hi))
rows$p_txt <- ifelse(is.na(rows$k), "", trimws(paste(p_txt, stars)))
est <- rows[!is.na(rows$k), ]
ticks <- c(0.5, 1.0, 2.0)
ylim <- c(0.4, nrow(rows) + 1.6)

forest <- ggplot(est) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey40") +
  geom_errorbarh(aes(y = y, xmin = lo, xmax = hi), height = 0, colour = "grey15") +
  geom_point(aes(x = HR, y = y), shape = 15, size = 3, colour = "#3B6DB3") +
  scale_x_log10(breaks = ticks, labels = trimws(formatC(ticks, format = "fg", digits = 3)),
                limits = range(c(est$lo, est$hi, ticks), finite = TRUE)) +
  scale_y_continuous(breaks = rows$y, labels = rows$label, limits = ylim, expand = c(0, 0)) +
  labs(x = "Hazard ratio (95% CI, log scale)", y = NULL,
       subtitle = sprintf("%d subjects, %d events", cox$n, cox$nevent)) +
  theme_classic(base_size = 13) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.y = element_text(hjust = 0, size = 10.5))
table <- ggplot(rows, aes(y = y)) +
  geom_text(aes(x = 0, label = hr_txt), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.66, label = p_txt), hjust = 0, size = 3.5) +
  annotate("text", x = c(0, 0.66), y = nrow(rows) + 1, label = c("HR (95% CI)", "P"),
           hjust = 0, fontface = "bold", size = 3.5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = ylim, expand = c(0, 0)) +
  labs(subtitle = " ", x = " ") + theme_void() +
  theme(axis.title.x = element_text(size = 13, colour = "white"),
        axis.text.x = element_text(size = 10.4, colour = "white"))
p <- ggpubr::ggarrange(forest, table, ncol = 2, align = "h", widths = c(4.32, 2.5))

ggsave("figure_hazard_ratios_cox.pdf", p, width = 6.82, height = 4.30, device = cairo_pdf)
cat("Reproduced figure.\n")
