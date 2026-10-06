#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_target_lesion_mm.csv")
raw = raw.dropna(subset=["patient", "week", "target_lesion_mm"]).sort_values(["patient", "week"])
base = raw.groupby("patient")["target_lesion_mm"].first()
raw["change"] = 100 * (raw["target_lesion_mm"] - raw["patient"].map(base)) / raw["patient"].map(base)
best = raw[raw["patient"].duplicated()].groupby("patient")["change"].min()
best = best.sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(6, 4.4))
ax.bar(np.arange(1, len(best) + 1), best.to_numpy())
ax.axhline(0, color="k"); ax.axhline(20.0, ls="--", color="grey"); ax.axhline(-30.0, ls="--", color="grey")
ax.set_xlabel("Patient"); ax.set_ylabel("Best change from baseline (%)")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("figure_target_lesion_mm_waterfall_py.pdf")
print(best)
print("Reproduced figure.")
