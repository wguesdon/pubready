#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("input_cell_type.csv")
d = pd.DataFrame({"id": raw["sample"].astype(str), "category": raw["cell_type"].astype(str),
                  "group": raw["condition"].astype(str), "count": pd.to_numeric(raw["n_cells"], errors="coerce")}).dropna()
cats = list(pd.unique(d["category"]))
tab = d.pivot_table(index="id", columns="category", values="count", aggfunc="sum", fill_value=0)[cats]
tab = tab[tab.sum(axis=1) > 0]
prop = tab.div(tab.sum(axis=1), axis=0)
g = d.drop_duplicates("id").set_index("id").loc[prop.index, "group"].to_numpy()
lv = sorted(set(g))
for c in cats if len(lv) >= 2 else []:
    v = prop[c].to_numpy()
    if len(lv) == 2:
        # The p-value rule of R's wilcox.test: exact below 50 per group without ties.
        a, b = v[g == lv[0]], v[g == lv[1]]
        exact = len(a) < 50 and len(b) < 50 and len(np.unique(v)) == len(v)
        print(c, stats.mannwhitneyu(a, b, alternative="two-sided",
                                    method="exact" if exact else "asymptotic", use_continuity=True))
    else:
        print(c, stats.kruskal(*[v[g == x] for x in lv]))
order = np.argsort(g, kind="stable")
prop = prop.iloc[order]
fig, ax = plt.subplots(figsize=(max(4.8, 2.5 + 0.3 * len(prop)), 4.8))
bottom = np.zeros(len(prop))
for c in reversed(cats):
    ax.bar(np.arange(len(prop)), prop[c].to_numpy(), bottom=bottom, width=0.85, edgecolor="white", label=c)
    bottom += prop[c].to_numpy()
ax.set_xticks(np.arange(len(prop))); ax.set_xticklabels(list(prop.index), rotation=90)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0)); ax.set_ylabel("Proportion")
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1], title="cell_type", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("figure_cell_type_comp_py.pdf")
print("Reproduced figure.")
