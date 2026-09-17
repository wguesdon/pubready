#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_pca_samples.csv")
grp = raw["group"].astype(str).to_numpy()
num = raw.select_dtypes(include="number")
X = num.to_numpy(float); X = X[:, X.std(0, ddof=1) > 0]
Xs = (X - X.mean(0)) / X.std(0, ddof=1)
U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
sc = U * S; ve = S**2 / np.sum(S**2) * 100
fig, ax = plt.subplots(figsize=(5.2, 4.6))
for lv in pd.unique(grp):
    m = grp == lv
    ax.scatter(sc[m, 0], sc[m, 1], s=40, alpha=0.9, label=lv)
ax.set_xlabel(f"PC1 ({ve[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({ve[1]:.1f}%)")
ax.legend(frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_pca_samples_pca_py.pdf")
print("Reproduced figure.")
