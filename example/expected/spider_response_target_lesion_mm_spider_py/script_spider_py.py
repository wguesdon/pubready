#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

raw = pd.read_csv("input_target_lesion_mm.csv")
raw["week"] = pd.to_numeric(raw["week"], errors="coerce")
raw["target_lesion_mm"] = pd.to_numeric(raw["target_lesion_mm"], errors="coerce")
raw = raw.dropna(subset=["patient", "week", "target_lesion_mm"]).sort_values(["patient", "week"])
base = raw.groupby("patient")["target_lesion_mm"].first()
raw["change"] = 100 * (raw["target_lesion_mm"] - raw["patient"].map(base)) / raw["patient"].map(base)

fig, ax = plt.subplots(figsize=(5.8, 4.2))
ax.axhline(0, color="grey", lw=0.8)
ax.axhline(20.0, color="grey", lw=0.8, ls="--")
ax.axhline(-30.0, color="grey", lw=0.8, ls="--")
for _, rows in raw.groupby("patient"):
    ax.plot(rows["week"], rows["change"], lw=1.1)
    ax.scatter(rows["week"], rows["change"], s=14)
ax.set_xlabel("week"); ax.set_ylabel("Change from baseline (%)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_target_lesion_mm_spider_py.pdf")
print("Reproduced figure.")
