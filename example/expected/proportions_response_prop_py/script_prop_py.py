#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np, pandas as pd
from scipy import stats as ss

raw = pd.read_csv("input_response.csv").dropna(subset=["arm", "response"])
tab = pd.crosstab(raw["arm"].astype(str), raw["response"].astype(str))
print(ss.chi2_contingency(tab.to_numpy(), correction=False))
prop = tab.div(tab.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(4.2, 4.2))
xpos = np.arange(tab.shape[0]); bottom = np.zeros(tab.shape[0])
for oc in prop.columns:
    ax.bar(xpos, prop[oc].to_numpy(), bottom=bottom, width=0.7, edgecolor="white", label=oc)
    bottom += prop[oc].to_numpy()
ax.set_xticks(xpos); ax.set_xticklabels(list(prop.index))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0)); ax.set_ylabel("proportion")
ax.legend(title="response", frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_response_prop_py.pdf")
print("Reproduced figure.")
