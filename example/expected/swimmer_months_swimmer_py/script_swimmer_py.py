#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_months.csv").sort_values("months").reset_index(drop=True)
raw["ongoing"] = raw["ended"] == 0
raw["colour"] = raw["arm"].astype("category").cat.codes
fig, ax = plt.subplots(figsize=(6.40, max(3.6, min(10, 1.6 + 0.18 * len(raw)))))
ax.barh(np.arange(1, len(raw) + 1), raw["months"], color=plt.cm.tab10(raw["colour"]))
ongoing = raw[raw["ongoing"]]
ax.scatter(ongoing["months"], ongoing.index + 1, marker=">", color="black")
ax.set_yticks(np.arange(1, len(raw) + 1)); ax.set_yticklabels(raw["patient"])
ax.set_xlabel("months"); ax.set_ylabel("Patient")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("figure_months_swimmer_py.pdf")
print("Reproduced figure.")
