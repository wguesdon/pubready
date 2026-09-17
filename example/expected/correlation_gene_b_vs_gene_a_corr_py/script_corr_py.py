#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, seaborn as sns, pingouin as pg

raw = pd.read_csv("input_gene_b_vs_gene_a.csv")
x = pd.to_numeric(raw["gene_a"], errors="coerce")
y = pd.to_numeric(raw["gene_b"], errors="coerce")
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask].to_numpy(), y[mask].to_numpy()
print(pg.corr(x, y, method="pearson"))

fig, ax = plt.subplots(figsize=(4.6, 4.2))
sns.regplot(x=x, y=y, ax=ax, color="#3B6DB3",
            scatter_kws={"s": 24, "alpha": 0.8, "edgecolor": "none"},
            line_kws={"color": "#333333", "linewidth": 1.4})
ax.set_xlabel("Gene A (a.u.)"); ax.set_ylabel("Gene B (a.u.)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figure_gene_b_vs_gene_a_corr_py.pdf")
print("Reproduced figure.")
