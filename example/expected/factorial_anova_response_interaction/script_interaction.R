#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({
  library(ggplot2); library(readr); library(ARTool)
})

df <- as.data.frame(read_csv("input_response.csv", show_col_types = FALSE))
df[["genotype"]] <- factor(df[["genotype"]], levels = c("ko", "wt"))
df[["treatment"]] <- factor(df[["treatment"]], levels = c("high", "low", "vehicle"))
df[["response"]] <- as.numeric(df[["response"]])

# The models run on a copy with syntactic column names; the effects take the
# names of the columns back.
safe <- c(".f1" = "genotype", ".f2" = "treatment")
sdf <- setNames(data.frame(df[["response"]], df[c("genotype", "treatment")]), c(".y", names(safe)))
at <- as.data.frame(car::Anova(lm(.y ~ .f1 * .f2, data = sdf), type = 2))
at <- at[rownames(at) != "Residuals", ]
eff <- data.frame(effect = rownames(at), p = at[["Pr(>F)"]])
term_name <- function(t) paste(safe[strsplit(t, ":", fixed = TRUE)[[1]]], collapse = ":")
eff$effect <- vapply(eff$effect, term_name, character(1), USE.NAMES = FALSE)
plab <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else paste0("P = ", formatC(p, format = "g", digits = 2))
sig <- eff[which(eff$p < 0.05), , drop = FALSE]
subtitle_txt <- if (nrow(sig) > 0) paste(sprintf("%s: %s", gsub(":", " × ", sig$effect), vapply(sig$p, plab, character(1))), collapse = "\n") else "No significant effects (P ≥ 0.05)"

# The mean and SEM of each cell.
by <- df[c("genotype", "treatment")]
summ <- aggregate(df["response"], by = by, FUN = mean); names(summ)[ncol(summ)] <- "mean"
summ$sem <- aggregate(df["response"], by = by, FUN = sd)[["response"]] / sqrt(aggregate(df["response"], by = by, FUN = length)[["response"]])

pd <- position_dodge(width = 0.25)
p <- ggplot(summ, aes(x = .data[["genotype"]], y = mean, colour = .data[["treatment"]], group = .data[["treatment"]])) +
  geom_line(position = pd, linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.12, position = pd, linewidth = 0.6) +
  geom_point(position = pd, size = 2.4) +
  scale_colour_manual(values = c("#3B6DB3", "#C1432B", "#2E8B57")) +
  labs(x = "genotype", y = "Response (a.u.)", colour = "treatment", subtitle = subtitle_txt) +
  theme_classic(base_size = 13) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9.5, colour = "grey30"))

ggsave("figure_response_interaction.pdf", p, width = 5.2, height = 4.4)
cat("Reproduced figure.\n")
