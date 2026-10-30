#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

df = pd.read_csv("input_arm.csv")
# The group labels as the recipe writes them: a whole number with no ".0".
df["arm"] = df["arm"].astype(str).str.replace(r"\.0$", "", regex=True)
order, pal = ['control', 'treated'], ['#3B6DB3', '#C1432B']
tmax = df["time"].max()
breaks = [5.0 * k for k in range(int(np.floor(tmax / 5.0 + 1e-9)) + 1)]
fig, (ax, axr) = plt.subplots(2, 1, figsize=(5.67, 5.8),
                              gridspec_kw={"height_ratios": [0.74, 0.26]})
for i, grp in enumerate(order):
    sub = df[df["arm"] == grp]
    kmf = KaplanMeierFitter().fit(sub["time"], sub["event"], label=grp)
    kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i], show_censors=True,
                               censor_styles={"marker": "|", "ms": 6})
    # At risk at t: the subjects whose time is t or later, as in R.
    for b in breaks:
        axr.text(b, len(order) - i, str(int((sub["time"] >= b).sum())), ha="center",
                 va="center", fontsize=9, color=pal[i])
res = multivariate_logrank_test(df["time"], df["arm"], df["event"])
p = res.p_value
ax.set_title("Log-rank P < 0.001" if p < 0.001 else "Log-rank P = %.2g" % p, fontsize=10)
ax.set_xlabel("Months")
ax.set_ylabel("Survival probability")
ax.set_ylim(0, 1.02); ax.set_xlim(-0.05 * tmax, 1.05 * tmax); ax.set_xticks(breaks)
ax.legend(title="arm", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
for s in ("top", "right"): ax.spines[s].set_visible(False)
axr.set_xlim(ax.get_xlim()); axr.set_ylim(0.5, len(order) + 0.5); axr.axis("off")
axr.set_title("Number at risk", loc="left", fontsize=10, fontweight="bold")
fig.tight_layout(); fig.savefig("figure_arm_km_py.pdf", bbox_inches="tight")
print("Reproduced figure.")
