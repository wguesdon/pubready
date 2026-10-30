#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages(library(ggplot2))

raw <- read.csv("input_rate.csv", check.names = FALSE)
d <- data.frame(s = as.numeric(raw[["substrate_um"]]), v = as.numeric(raw[["rate"]]),
                group = as.character(raw[["inhibitor"]]))
d <- d[is.finite(d$s) & is.finite(d$v) & !is.na(d$group) & d$s >= 0, ]
starts <- list(
  "0 uM" = list(vmax = 114.31, km = 20),
  "2 uM" = list(vmax = 111.75, km = 20),
  "10 uM" = list(vmax = 105.73, km = 40)
)
algos <- c("0 uM" = "default", "2 uM" = "default", "10 uM" = "default")
ctrl <- nls.control(maxiter = 1000, tol = 1e-8, minFactor = 1e-12)
curves <- NULL
for (g in names(starts)) {
  fit <- nls(v ~ vmax * s / (km + s), data = d[d$group == g, ], start = starts[[g]],
             control = ctrl, algorithm = algos[[g]])
  print(summary(fit))
  grid <- seq(0, max(d$s), length.out = 200)
  curves <- rbind(curves, data.frame(group = g, s = grid, v = predict(fit, data.frame(s = grid))))
}
# A curve that is not determined is drawn as points only, as in the figure.
p <- ggplot() +
  geom_point(data = d, aes(s, v, colour = group), alpha = 0.6) +
  expand_limits(x = 0, y = 0) + labs(x = "Substrate (uM)", y = "Rate (nmol/min/mg)") + theme_classic()
if (!is.null(curves)) p <- p + geom_line(data = curves, aes(s, v, colour = group))
ggsave("figure_rate_enzyme.pdf", p, width = 6.6, height = 4.4)
cat("Reproduced figure.\n")
