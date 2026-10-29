#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_viability_pct.csv", check.names = FALSE)
d <- data.frame(dose = as.numeric(raw[["concentration_um"]]), y = as.numeric(raw[["viability_pct"]]),
                group = as.character(raw[["compound"]]))
d <- d[is.finite(d$dose) & is.finite(d$y) & d$dose > 0, ]
d$lx <- log10(d$dose)
# Only the curves with a determined midpoint are refitted; the others are
# drawn as points, as in the figure.
starts <- list(
  "cmpd_A" = list(bottom = 2.966666667, top = 98.86666667, log_ec50 = -1.522878745, hill = -1),
  "cmpd_B" = list(bottom = 11.4, top = 98.9, log_ec50 = -0.5228787453, hill = -1),
  "cmpd_C" = list(bottom = 23.13333333, top = 102.3, log_ec50 = 0, hill = -1)
)
algos <- c("cmpd_A" = "default", "cmpd_B" = "port", "cmpd_C" = "port")
fml <- y ~ bottom + (top - bottom) / (1 + 10^((log_ec50 - lx) * hill))
ctrl <- nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)
curves <- NULL
for (g in names(starts)) {
  fit <- nls(fml, data = d[d$group == g, ], start = starts[[g]], control = ctrl,
             algorithm = algos[[g]])
  print(summary(fit))
  grid <- seq(min(d$lx[d$group == g]), max(d$lx[d$group == g]), length.out = 200)
  curves <- rbind(curves, data.frame(group = g, dose = 10^grid,
                                     y = predict(fit, data.frame(lx = grid))))
}
p <- ggplot() + geom_point(data = d, aes(dose, y, colour = group), alpha = 0.6) +
  scale_x_log10() + labs(x = "Concentration (uM)", y = "Viability (%)") + theme_classic()
if (!is.null(curves)) p <- p + geom_line(data = curves, aes(dose, y, colour = group))
ggsave("figure_viability_pct_drc.pdf", p, width = 6.2, height = 4.2)
cat("Reproduced figure.\n")
