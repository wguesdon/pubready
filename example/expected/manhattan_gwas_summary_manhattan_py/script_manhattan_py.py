#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("input_gwas_summary.csv")
chrs = raw["chr"].astype(str).str.replace(r"^(?i:chr)", "", regex=True).to_numpy()
pos = pd.to_numeric(raw["pos"], errors="coerce").to_numpy(float)
p = pd.to_numeric(raw["p"], errors="coerce").to_numpy(float)
keep = np.isfinite(pos) & np.isfinite(p)
chrs, pos, p = chrs[keep], pos[keep], np.where(p[keep] == 0, np.finfo(float).tiny, p[keep])
special = ["X", "Y", "XY", "MT", "M"]
key = {c: ((1, float(c), "") if c.isdigit() else (2, float(special.index(c.upper()) + 1), "")
           if c.upper() in special else (3, 0.0, c)) for c in set(chrs)}
lv = sorted(set(chrs), key=key.get)
ci = np.array([lv.index(c) for c in chrs])
lo = np.array([pos[ci == k].min() for k in range(len(lv))])
span = np.array([pos[ci == k].max() for k in range(len(lv))]) - lo
offset = np.concatenate([[0.0], np.cumsum(span + 0.004 * span.sum())])[:len(lv)]
x = pos - lo[ci] + offset[ci]
lam = np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)
fig, (ax, aq) = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw={"width_ratios": [2.7, 1]})
ax.scatter(x, -np.log10(p), s=4, c=np.where(ci % 2 == 0, "#3B6DB3", "#1F3A5F"), linewidths=0)
ax.axhline(-np.log10(1e-5), ls="--", color="grey"); ax.axhline(-np.log10(5e-8), ls="--", color="#C1432B")
ax.set_xticks(offset + span / 2); ax.set_xticklabels(lv, fontsize=8)
ax.set_xlabel("Chromosome"); ax.set_ylabel("-log10 p")
n = len(p); i = np.arange(1, n + 1)
e = -np.log10((i - 0.5) / n)
aq.fill_between(e, -np.log10(stats.beta.ppf(0.975, i, n - i + 1)),
                -np.log10(stats.beta.ppf(0.025, i, n - i + 1)), color="#D9D9D9")
aq.plot([0, e.max()], [0, e.max()], color="#C1432B")
aq.scatter(e, -np.log10(np.sort(p)), s=4, c="#3B6DB3", linewidths=0)
aq.set_xlabel("Expected -log10 p"); aq.set_ylabel("Observed -log10 p")
aq.set_title(f"lambda_GC = {lam:.3f}")
for a in (ax, aq):
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_gwas_summary_manhattan_py.pdf")
print("Reproduced figure.")
