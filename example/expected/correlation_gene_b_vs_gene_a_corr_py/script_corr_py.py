#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats


def lm_band(x, y, n=80):
    """The least-squares line and its 95% confidence band over the range of x.

    The band is the t-based interval of the fitted mean, the band that
    geom_smooth(method = "lm") draws in the R engine on 80 points. It has no
    random part, so the figure is the same in every run.

    Args:
        x: The x values.
        y: The y values.
        n: The count of grid points.

    Returns:
        tuple: The grid, the fitted line, the lower bound and the upper bound.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    grid = np.linspace(x.min(), x.max(), n)
    fit = intercept + slope * grid
    dof = len(x) - 2
    s = np.sqrt(np.sum((y - (intercept + slope * x)) ** 2) / dof)
    se = s * np.sqrt(1 / len(x) + (grid - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    t = stats.t.ppf(0.975, dof)
    return grid, fit, fit - t * se, fit + t * se


raw = pd.read_csv("input_gene_b_vs_gene_a.csv")
x = pd.to_numeric(raw["gene_a"], errors="coerce")
y = pd.to_numeric(raw["gene_b"], errors="coerce")
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask].to_numpy(), y[mask].to_numpy()
res = pg.corr(x, y, method="pearson"); print(res)
sub = "Pearson correlation, r = %.2f, P = %.2g" % (res["r"].iloc[0], res["p_val"].iloc[0])

fig, ax = plt.subplots(figsize=(4.6, 4.2))
ax.scatter(x, y, s=18.22, color="#3B6DB3", alpha=0.8, edgecolors="none")
grid, fit, lo, hi = lm_band(x, y)
ax.fill_between(grid, lo, hi, color="#BFBFBF", alpha=0.4, linewidth=0)
ax.plot(grid, fit, color="#333333", linewidth=1.4)
ax.set_xlabel("Gene A (a.u.)"); ax.set_ylabel("Gene B (a.u.)")
ax.set_title(sub, fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figure_gene_b_vs_gene_a_corr_py.pdf")
print("Reproduced figure.")
