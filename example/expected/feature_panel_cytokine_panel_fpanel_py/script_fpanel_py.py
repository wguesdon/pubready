#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def label(v):
    """Write a group value as the figure does: a whole number with no decimal point."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


raw = pd.read_csv("input_cytokine_panel.csv")
d = pd.DataFrame({"feature": raw["cytokine"].astype(str), "group": raw["group"].map(label), "value": pd.to_numeric(raw["concentration_pg_ml"], errors="coerce")})
d = d[np.isfinite(d["value"])]
feats = list(pd.unique(d["feature"]))
levels = ['Healthy control', 'Mild disease', 'Severe disease']
p = []
for f in feats:
    g = [d.loc[(d["feature"] == f) & (d["group"] == lv), "value"].to_numpy() for lv in levels]
    p.append(stats.kruskal(*g).pvalue)
q = multipletests(p, method="fdr_bh")[1]
nrow = math.ceil(len(feats) / 3)
fig, axes = plt.subplots(nrow, 3, figsize=(3 * 3, 2.6 * nrow), squeeze=False)
for i, f in enumerate(feats):
    ax = axes[i // 3][i % 3]
    g = [d.loc[(d["feature"] == f) & (d["group"] == lv), "value"].to_numpy() for lv in levels]
    ax.boxplot(g, showfliers=False)
    ax.set_xticks(range(1, len(levels) + 1)); ax.set_xticklabels(levels, fontsize=8)
    ax.set_title(f"{f}\nq = {q[i]:.2g}", fontsize=10)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_cytokine_panel_fpanel_py.pdf")
print("Reproduced figure.")
