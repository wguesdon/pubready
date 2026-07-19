#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_enrichment_results.csv")
term = raw["Description"].astype(str).to_numpy()
padj = pd.to_numeric(raw["p.adjust"], errors="coerce").to_numpy()
ct = pd.to_numeric(raw["Count"], errors="coerce").to_numpy()
x = raw["GeneRatio"].astype(str).map(lambda s: float(s.split("/")[0])/float(s.split("/")[1]) if "/" in s else float(s)).to_numpy(); xlab = "gene ratio"
idx = np.argsort(padj)[:15]
t, x, pj, ct = term[idx], x[idx], padj[idx], ct[idx]
yo = np.argsort(x); t, x, pj, ct = t[yo], x[yo], pj[yo], ct[yo]
fig, ax = plt.subplots(figsize=(7.2, max(4.0, 0.34*len(t)+1.6)))
sc = ax.scatter(x, np.arange(len(t)), c=pj, s=70, cmap="coolwarm_r", edgecolors="none")
ax.set_yticks(np.arange(len(t))); ax.set_yticklabels(t)
ax.set_xlabel(xlab)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.colorbar(sc, ax=ax, label="p.adjust", fraction=0.046, pad=0.04)
fig.tight_layout(); fig.savefig("figure_enrichment_results_enrich_py.pdf")
print("Reproduced figure.")
