#!/usr/bin/env Rscript
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
suppressPackageStartupMessages({ library(ggplot2) })

raw <- read.csv("input_cell_area_um2.csv", check.names = FALSE)
d <- data.frame(group = as.character(raw[["condition"]]), value = as.numeric(raw[["cell_area_um2"]]))
d <- d[is.finite(d$value) & !is.na(d$group), ]

# Freedman-Diaconis bins on the pooled values (Sturges when the IQR is zero).
v <- d$value; n <- length(v)
iqr <- diff(quantile(v, c(0.25, 0.75), type = 7, names = FALSE))
w <- if (iqr > 0) 2 * iqr * n^(-1/3) else if (max(v) > min(v)) (max(v) - min(v)) / ceiling(log2(n) + 1) else 1
if (max(v) > min(v) && (max(v) - min(v)) / w > 200) w <- (max(v) - min(v)) / 200
start <- floor(min(v) / w) * w; nb <- floor((max(v) - start) / w) + 1

bars <- NULL; curves <- NULL
for (g in sort(unique(d$group), method = "radix")) {
  vg <- d$value[d$group == g]
  idx <- pmin(pmax(floor((vg - start) / w) + 1, 1), nb)
  bars <- rbind(bars, data.frame(group = g, xmin = start + (seq_len(nb) - 1) * w,
                                 xmax = start + seq_len(nb) * w, count = tabulate(idx, nb)))
  if (length(vg) < 2) next
  bw <- bw.nrd0(vg)
  grid <- seq(min(vg) - 3 * bw, max(vg) + 3 * bw, length.out = 512)
  dens <- rowSums(dnorm(outer(grid, vg, "-") / bw)) / (length(vg) * bw)
  curves <- rbind(curves,
    data.frame(group = g, x = grid, y = dens * length(vg) * w, curve = "Kernel density"),
    data.frame(group = g, x = grid, y = dnorm(grid, mean(vg), sd(vg)) * length(vg) * w, curve = "Normal fit"))
  if (length(vg) >= 3 && length(vg) <= 5000) print(shapiro.test(vg))
}

p <- ggplot() +
  geom_rect(data = bars, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = group),
            colour = "white", linewidth = 0.25, alpha = 0.65) +
  geom_line(data = curves, aes(x = x, y = y, linetype = curve), colour = "grey15", linewidth = 0.7) +
  scale_linetype_manual(values = c("Kernel density" = "solid", "Normal fit" = "22"), name = NULL) +
  labs(x = "Cell area (um^2)", y = "Count") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom", strip.background = element_blank()) +
  guides(fill = "none")
p <- p + facet_wrap(~ group, ncol = 1, scales = "free_y")

ggsave("figure_cell_area_um2_dist.pdf", p, width = 5.6, height = 1.5 + 1.6 * length(unique(d$group)))
cat("Reproduced figure.\n")
