#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats as ss

raw = pd.read_csv("input_value.csv")
raw["value"] = pd.to_numeric(raw["value"], errors="coerce")
raw = raw.dropna(subset=["subject", "condition", "value"])
raw["condition"] = raw["condition"].astype(str)
lv = ['pre', 'post']
a = raw[raw["condition"] == lv[0]][["subject", "value"]].drop_duplicates("subject").rename(columns={"value": "v1"})
b = raw[raw["condition"] == lv[1]][["subject", "value"]].drop_duplicates("subject").rename(columns={"value": "v2"})
m = a.merge(b, on="subject")
res = ss.ttest_rel(m["v2"], m["v1"])
print(res)
fig, ax = plt.subplots(figsize=(3.8, 4.2))
for _, r in m.iterrows(): ax.plot([0, 1], [r["v1"], r["v2"]], color="grey", lw=0.6, alpha=0.7)
ax.scatter([0]*len(m), m["v1"], s=45); ax.scatter([1]*len(m), m["v2"], s=45)
ax.set_xticks([0, 1]); ax.set_xticklabels(lv); ax.set_xlabel("condition"); ax.set_ylabel("value")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_value_paired_py.pdf")
print("Reproduced figure.")
