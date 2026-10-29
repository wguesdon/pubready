#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("input_cell_area_um2.csv")
d = pd.DataFrame({"group": raw["condition"].astype(str), "value": pd.to_numeric(raw["cell_area_um2"], errors="coerce")})
d = d[np.isfinite(d["value"]) & d["group"].notna()]
v = d["value"].to_numpy(); n = len(v)
iqr = np.quantile(v, 0.75) - np.quantile(v, 0.25)
w = 2 * iqr * n ** (-1 / 3) if iqr > 0 else ((v.max() - v.min()) / math.ceil(math.log2(n) + 1) if v.max() > v.min() else 1.0)
if v.max() > v.min() and (v.max() - v.min()) / w > 200:
    w = (v.max() - v.min()) / 200
start = math.floor(v.min() / w) * w; nb = int((v.max() - start) // w) + 1
groups = sorted(d["group"].unique())
fig, axes = plt.subplots(len(groups), 1, figsize=(5.6, 1.5 + 1.6 * len(groups)), sharex=True, squeeze=False)
for ax, g in zip(axes[:, 0], groups):
    vg = d.loc[d["group"] == g, "value"].to_numpy()
    idx = np.clip(np.floor((vg - start) / w).astype(int), 0, nb - 1)
    ax.bar(start + w * np.arange(nb), np.bincount(idx, minlength=nb), width=w, align="edge", alpha=0.65, edgecolor="white")
    if len(vg) >= 2:
        sd = vg.std(ddof=1); iq = np.quantile(vg, 0.75) - np.quantile(vg, 0.25)
        bw = 0.9 * (min(sd, iq / 1.34) or sd or 1.0) * len(vg) ** -0.2
        grid = np.linspace(vg.min() - 3 * bw, vg.max() + 3 * bw, 512)
        dens = np.exp(-0.5 * ((grid[:, None] - vg[None, :]) / bw) ** 2).sum(axis=1) / (len(vg) * bw * math.sqrt(2 * math.pi))
        ax.plot(grid, dens * len(vg) * w, color="#262626")
        ax.plot(grid, stats.norm.pdf(grid, vg.mean(), sd) * len(vg) * w, color="#262626", ls="--")
    if 3 <= len(vg) <= 5000:
        print(g, stats.shapiro(vg))
    ax.set_title(str(g) if True else "", loc="left", fontsize=10); ax.set_ylabel("Count")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1, 0].set_xlabel("Cell area (um^2)")
fig.tight_layout(); fig.savefig("figure_cell_area_um2_dist_py.pdf")
print("Reproduced figure.")
