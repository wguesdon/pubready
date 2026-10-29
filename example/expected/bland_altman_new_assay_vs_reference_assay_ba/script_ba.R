#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_new_assay_vs_reference_assay.csv", check.names = FALSE)
a <- as.numeric(raw[["new_assay"]]); b <- as.numeric(raw[["reference_assay"]])
keep <- is.finite(a) & is.finite(b); a <- a[keep]; b <- b[keep]
d <- data.frame(mean = (a + b) / 2, diff = a - b)
bias <- mean(d$diff); s <- sd(d$diff)
loa <- bias + c(-1, 1) * qnorm(0.975) * s
cat(sprintf("bias %.4f, SD %.4f, limits %.4f to %.4f\n", bias, s, loa[1], loa[2]))
print(summary(lm(diff ~ mean, data = d)))
p <- ggplot(d, aes(mean, diff)) + geom_point() +
  geom_hline(yintercept = bias) + geom_hline(yintercept = loa, linetype = "dashed") +
  labs(x = "Mean of new_assay and reference_assay", y = "Difference (new_assay - reference_assay)") + theme_classic()
ggsave("figure_new_assay_vs_reference_assay_ba.pdf", p, width = 6.4, height = 4.4)
cat("Reproduced figure.\n")
