#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_de_basemean.csv")
m = pd.to_numeric(raw["log2FoldChange"], errors="coerce").to_numpy(float)
p = pd.to_numeric(raw["padj"], errors="coerce").to_numpy(float)
a = np.log10(pd.to_numeric(raw["baseMean"], errors="coerce").to_numpy(float))
lab = raw["gene"].astype(str).to_numpy()
keep = np.isfinite(m) & np.isfinite(a)
m, p, a, lab = m[keep], p[keep], a[keep], lab[keep]
fcc, pcc = 1.0, 0.05
sig = ~np.isnan(p) & (p <= pcc)
up, down = sig & (m >= fcc), sig & (m <= -fcc)
hit = np.where(up | down)[0]
hit = hit[np.lexsort((hit, -np.abs(m[hit]), p[hit]))]
fig, ax = plt.subplots(figsize=(6, 5.4))
ns = ~(up | down)
ax.scatter(a[ns], m[ns], s=1.6, c="#B3B3B3", linewidths=0)
ax.axhline(0, color="black", lw=0.9)
ax.axhline(fcc, color="grey", ls="--", lw=0.7); ax.axhline(-fcc, color="grey", ls="--", lw=0.7)
ax.scatter(a[up], m[up], s=9, c="#C1432B", linewidths=0, label="Up")
ax.scatter(a[down], m[down], s=9, c="#3B6DB3", linewidths=0, label="Down")
for i in hit[:15]:
    ax.annotate(lab[i], (a[i], m[i]), fontsize=8)
lo, hi = min(m.min(), -1.5 * fcc), max(m.max(), 1.5 * fcc)
ax.set_ylim(lo - 0.06 * (hi - lo), hi + 0.06 * (hi - lo))
ax.set_xlabel("log10 mean of normalized counts"); ax.set_ylabel("log2 fold change")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figure_de_basemean_ma_py.pdf")
print("Reproduced figure.")
