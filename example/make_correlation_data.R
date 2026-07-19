#!/usr/bin/env Rscript
# Generate the correlation example data deterministically. Run once inside the
# pubplot container from the repo root:
#   podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.3.0 \
#     Rscript example/make_correlation_data.R
#
# Produces two files:
#   correlation_xy.csv   two co-expressed genes (for the `correlation` recipe)
#   correlation_vars.csv a panel of cytokines (for `correlation_heatmap`)
set.seed(23)

# Two variables with a clear positive correlation.
n1 <- 60
gene_a <- rnorm(n1, mean = 8, sd = 1.5)
gene_b <- 0.8 * gene_a + rnorm(n1, sd = 1.0) + 2
xy <- data.frame(gene_a = round(gene_a, 2), gene_b = round(gene_b, 2))
write.csv(xy, "example/correlation_xy.csv", row.names = FALSE)

# A panel of numeric variables with two correlated blocks (pro- and
# anti-inflammatory), so the clustered correlation heatmap has structure.
n2 <- 80
f_pro  <- rnorm(n2)
f_anti <- rnorm(n2)
vars <- data.frame(
  il6  = 2.0 * f_pro  + rnorm(n2, sd = 0.6),
  tnf  = 1.8 * f_pro  + rnorm(n2, sd = 0.7),
  crp  = 1.5 * f_pro  + rnorm(n2, sd = 0.8),
  il1b = 1.7 * f_pro  + rnorm(n2, sd = 0.7),
  il10 = 1.9 * f_anti - 0.4 * f_pro + rnorm(n2, sd = 0.6),
  tgfb = 1.6 * f_anti + rnorm(n2, sd = 0.7),
  il4  = 1.4 * f_anti + rnorm(n2, sd = 0.8),
  ifng = 1.3 * f_pro  - 0.6 * f_anti + rnorm(n2, sd = 0.8)
)
vars <- round(vars, 2)
write.csv(vars, "example/correlation_vars.csv", row.names = FALSE)

cat(sprintf("wrote example/correlation_xy.csv (%d rows) and correlation_vars.csv (%d x %d)\n",
            n1, n2, ncol(vars)))
