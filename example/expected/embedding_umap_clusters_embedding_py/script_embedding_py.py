#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

raw = pd.read_csv("input_umap_clusters.csv").dropna(subset=["UMAP_1", "UMAP_2", "cell_type"])
fig, ax = plt.subplots(figsize=(6.6, 5.2))
ax.scatter(raw["UMAP_1"], raw["UMAP_2"], s=1.64, linewidths=0, c=pd.factorize(raw["cell_type"].astype(str))[0], cmap="tab20")
if not False:
    for name, sub in raw.groupby(raw["cell_type"].astype(str)):
        ax.text(sub["UMAP_1"].median(), sub["UMAP_2"].median(), name, ha="center",
                bbox={"boxstyle": "round", "fc": "white", "ec": "none"})
ax.set_aspect("equal", adjustable="datalim"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP_1"); ax.set_ylabel("UMAP_2")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_umap_clusters_embedding_py.pdf")
print("Reproduced figure.")
