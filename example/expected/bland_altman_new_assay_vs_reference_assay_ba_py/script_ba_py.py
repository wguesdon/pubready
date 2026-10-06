#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("input_new_assay_vs_reference_assay.csv")
a = pd.to_numeric(raw["new_assay"], errors="coerce"); b = pd.to_numeric(raw["reference_assay"], errors="coerce")
keep = np.isfinite(a) & np.isfinite(b)
a, b = a[keep].to_numpy(), b[keep].to_numpy()
mean, diff = (a + b) / 2, a - b
bias, s = diff.mean(), diff.std(ddof=1)
loa = bias + np.array([-1, 1]) * stats.norm.ppf(0.975) * s
print(f"bias {bias:.4f}, SD {s:.4f}, limits {loa[0]:.4f} to {loa[1]:.4f}")
print(stats.linregress(mean, diff))
fig, ax = plt.subplots(figsize=(5.4, 4.4))
ax.scatter(mean, diff)
ax.axhline(bias, color="k"); [ax.axhline(v, color="k", ls="--") for v in loa]
ax.set_xlabel("Mean of new_assay and reference_assay"); ax.set_ylabel("Difference (new_assay - reference_assay)")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("figure_new_assay_vs_reference_assay_ba_py.pdf")
print("Reproduced figure.")
