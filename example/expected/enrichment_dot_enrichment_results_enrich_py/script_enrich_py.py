#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap
import numpy as np, pandas as pd

raw = pd.read_csv("input_enrichment_results.csv")
term = raw["Description"].astype(str).to_numpy()
padj = pd.to_numeric(raw["p.adjust"], errors="coerce").to_numpy()
ct = pd.to_numeric(raw["Count"], errors="coerce").to_numpy()
x = np.array([float(s.split("/")[0]) / float(s.split("/")[1]) if "/" in s else float(s) for s in raw["GeneRatio"].astype(str)]); xlab = "gene ratio"
idx = np.argsort(padj, kind="stable")[:15]
t, x, pj, ct = term[idx], x[idx], padj[idx], ct[idx]
yo = np.argsort(x, kind="stable"); t, x, pj, ct = t[yo], x[yo], pj[yo], ct[yo]
# Wrap each term at about 45 characters, as the figure does.
t = [textwrap.fill(s, 45, break_long_words=False, break_on_hyphens=False) for s in t]
fig, ax = plt.subplots(figsize=(7.72, 7.00))
sc = ax.scatter(x, np.arange(len(t)), c=-np.log10(pj), s=70, cmap="coolwarm", edgecolors="none")
ax.set_yticks(np.arange(len(t))); ax.set_yticklabels(t)
ax.set_xlabel(xlab)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.colorbar(sc, ax=ax, label="-log10 p.adjust", fraction=0.046, pad=0.04)
fig.tight_layout(); fig.savefig("figure_enrichment_results_enrich_py.pdf")
print("Reproduced figure.")
