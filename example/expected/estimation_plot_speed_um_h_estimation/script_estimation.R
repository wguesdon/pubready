#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder. It recomputes the BCa bootstrap
# intervals with the same Park-Miller resamples and draws a simple version of
# the estimation plot.
suppressPackageStartupMessages(library(ggplot2))

M <- 2147483647; A <- 48271; SEED <- 20261029; B <- 5000
mulmod <- function(z, a) {
  hi <- floor(a / 65536); lo <- a - hi * 65536
  ((((z * hi) %% M) * 65536) %% M + (z * lo) %% M) %% M
}
uniforms <- function(n) {
  s <- SEED; out <- numeric(n)
  for (i in seq_len(n)) { s <- mulmod(s, A); out[i] <- s }
  out / M
}
bca <- function(r, g, paired) {
  per <- if (paired) length(r) else length(r) + length(g)
  u <- matrix(uniforms(B * per), nrow = per)
  if (paired) {
    i <- floor(u * length(r)) + 1
    boot <- colSums(matrix(g[i] - r[i], nrow = length(r))) / length(r)
    d <- g - r; jack <- (sum(d) - d) / (length(d) - 1)
  } else {
    ir <- floor(u[seq_along(r), , drop = FALSE] * length(r)) + 1
    ig <- floor(u[length(r) + seq_along(g), , drop = FALSE] * length(g)) + 1
    boot <- colSums(matrix(g[ig], nrow = length(g))) / length(g) -
      colSums(matrix(r[ir], nrow = length(r))) / length(r)
    jack <- c(mean(g) - (sum(r) - r) / (length(r) - 1), (sum(g) - g) / (length(g) - 1) - mean(r))
  }
  est <- sum(g) / length(g) - sum(r) / length(r)
  p0 <- sum(boot < est - 1e-11 * max(abs(c(r, g)))) / B
  z0 <- qnorm(min(max(p0, 0.5 / B), 1 - 0.5 / B))
  dj <- mean(jack) - jack; acc <- sum(dj^3) / (6 * sum(dj^2)^1.5)
  zq <- qnorm(c(0.025, 0.975))
  c(est, quantile(boot, pnorm(z0 + (z0 + zq) / (1 - acc * (z0 + zq))), type = 7))
}

raw <- read.csv("input_speed_um_h.csv", check.names = FALSE)
d <- data.frame(group = as.character(raw[["treatment"]]), y = as.numeric(raw[["speed_um_h"]]))
d <- d[!is.na(d$group) & is.finite(d$y), ]
lv <- c("vehicle", "compound_A", "compound_B", "compound_C")
res <- NULL
for (g in lv[-1]) {
  r <- d$y[d$group == lv[1]]; v <- d$y[d$group == g]
  ci <- bca(r, v, FALSE)
  tt <- t.test(v, r, paired = FALSE)
  res <- rbind(res, data.frame(group = g, mean_diff = ci[1], ci_low = ci[2], ci_high = ci[3],
                               p_value = tt$p.value))
}
print(res, digits = 6)
d$group <- factor(d$group, levels = lv)
p1 <- ggplot(d, aes(group, y, colour = group)) + geom_jitter(width = 0.15) +
  labs(x = NULL, y = "Migration speed (um/h)") + theme_classic() + theme(legend.position = "none")
res$group <- factor(res$group, levels = lv)
p2 <- ggplot(res, aes(group, mean_diff)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = "Mean difference vs vehicle") + theme_classic()
ggsave("figure_speed_um_h_estimation.pdf", ggpubr::ggarrange(p1, p2, ncol = 1, align = "v"), width = 5, height = 6)
cat("Reproduced figure.\n")
