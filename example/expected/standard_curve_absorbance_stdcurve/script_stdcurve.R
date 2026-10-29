#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_absorbance.csv", check.names = FALSE)
d <- data.frame(conc = suppressWarnings(as.numeric(as.character(raw[["conc_ng_ml"]]))), sig = as.numeric(raw[["absorbance"]]))
d$row <- seq_len(nrow(d))
std <- d[is.finite(d$conc) & is.finite(d$sig), ]
unk <- d[is.na(d$conc) & is.finite(d$sig), ]
st <- list(bottom = 0.1565, top = 2.7130000000000001, log_ec50 = 0.3979400086720376, hill = 1)
pos <- std[std$conc > 0, ]; pos$lx <- log10(pos$conc)
fit <- nls(sig ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill)), data = pos, start = st,
           control = nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12), algorithm = "default")
cf <- coef(fit); if (cf[['bottom']] > cf[['top']]) cf <- c(bottom = cf[['top']], top = cf[['bottom']], log_ec50 = cf[['log_ec50']], hill = -cf[['hill']])
f <- function(c) cf[['bottom']] + (cf[['top']] - cf[['bottom']]) / (1 + 10^((cf[['log_ec50']] - log10(c)) * cf[['hill']]))
finv <- function(s) 10^(cf[['log_ec50']] - log10((cf[['top']] - cf[['bottom']]) / (s - cf[['bottom']]) - 1) / cf[['hill']])
rng <- range(pos$conc)
print(summary(fit))
ends <- f(rng)
inside <- unk$sig >= min(ends) & unk$sig <= max(ends)
unk$concentration <- ifelse(inside, finv(unk$sig), NA)
print(unk)
grid <- data.frame(conc = seq(rng[1], rng[2], length.out = 200)); grid$sig <- f(grid$conc)
p <- ggplot() + geom_point(data = std, aes(conc, sig)) + geom_line(data = grid, aes(conc, sig)) +
  geom_point(data = unk[inside, ], aes(concentration, sig), shape = 23) +
  scale_x_log10() + labs(x = "IL-6 (ng/mL)", y = "Absorbance at 450 nm") + theme_classic()
ggsave("figure_absorbance_stdcurve.pdf", p, width = 5.2, height = 4.8)
cat("Reproduced figure.\n")
