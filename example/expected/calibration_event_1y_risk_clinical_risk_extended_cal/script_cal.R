#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ResourceSelection); library(ggplot2) })

raw <- read.csv("input_event_1y_risk_clinical_risk_extended.csv", check.names = FALSE, stringsAsFactors = FALSE)
models <- c("risk_clinical", "risk_extended")
raw <- raw[!is.na(raw[["event_1y"]]) & complete.cases(raw[, models, drop = FALSE]), ]
yb <- as.integer(as.character(raw[["event_1y"]]) == "yes")
curves <- NULL
for (m in models) {
  p <- raw[[m]]
  lp <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
  cat(m, "Brier", mean((p - yb)^2), "\n")
  print(coef(summary(glm(yb ~ 1 + offset(lp), family = binomial))))
  print(coef(summary(glm(yb ~ lp, family = binomial))))
  print(hoslem.test(yb, p, g = 10))
  sm <- lowess(p, yb, f = 2 / 3, iter = 0, delta = 0.01 * diff(range(p)))
  curves <- rbind(curves, data.frame(model = m, x = sm$x, y = sm$y))
}
g <- ggplot(curves, aes(x, y, colour = model)) + geom_abline(linetype = "dashed") + geom_path() +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability", y = "Observed proportion") + theme_classic(base_size = 13)
ggsave("figure_event_1y_risk_clinical_risk_extended_cal.pdf", g, width = 5.20, height = 8.56)
cat("Reproduced figure.\n")
