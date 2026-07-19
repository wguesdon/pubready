#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_de_results.csv")
fc = pd.to_numeric(raw["log2FoldChange"], errors="coerce").to_numpy()
p = pd.to_numeric(raw["padj"], errors="coerce").to_numpy()
lab = raw["gene"].astype(str).to_numpy()
m = np.isfinite(fc) & np.isfinite(p) & (p > 0)
fc, p, lab = fc[m], p[m], lab[m]
neglog = -np.log10(p)
fcc, pcc = 1.0, 0.05
cat = np.where((np.abs(fc) >= fcc) & (p <= pcc), "FC and p",
      np.where(p <= pcc, "p", np.where(np.abs(fc) >= fcc, "FC", "NS")))
colors = {"NS": "grey", "FC": "#2E8B57", "p": "#3B6DB3", "FC and p": "#C1432B"}
fig, ax = plt.subplots(figsize=(6, 5.6))
for c in ("NS", "FC", "p", "FC and p"):
    mm = cat == c
    if mm.any():
        ax.scatter(fc[mm], neglog[mm], s=10, alpha=0.7, c=colors[c], label=c, edgecolors="none")
ax.axvline(-fcc, ls="--", lw=0.6, color="grey"); ax.axvline(fcc, ls="--", lw=0.6, color="grey")
ax.axhline(-np.log10(pcc), ls="--", lw=0.6, color="grey")
ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 p")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False)
fig.tight_layout(); fig.savefig("figure_de_results_volcano_py.pdf")
print("Reproduced figure.")
