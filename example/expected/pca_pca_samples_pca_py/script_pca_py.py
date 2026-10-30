#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_pca_samples.csv")
# The features of the figure; a cell that does not read as a number is NaN,
# and a sample with a missing feature or group is dropped, as in the figure.
feats = ['feature01', 'feature02', 'feature03', 'feature04', 'feature05', 'feature06', 'feature07', 'feature08', 'feature09', 'feature10', 'feature11', 'feature12', 'feature13', 'feature14', 'feature15', 'feature16', 'feature17', 'feature18', 'feature19', 'feature20', 'feature21', 'feature22', 'feature23', 'feature24', 'feature25', 'feature26', 'feature27', 'feature28', 'feature29', 'feature30', 'feature31', 'feature32', 'feature33', 'feature34', 'feature35', 'feature36', 'feature37', 'feature38', 'feature39', 'feature40']
num = raw[feats].apply(pd.to_numeric, errors="coerce")
g = raw["group"]
keep = (num.notna().all(axis=1) & g.notna() & (g.astype(str).str.strip() != "")).to_numpy()
grp = g.astype(str).to_numpy()[keep]
X = num.to_numpy(float)[keep]
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
