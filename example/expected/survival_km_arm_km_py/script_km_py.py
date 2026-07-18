#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test

df = pd.read_csv("input_arm.csv"); df["arm"] = df["arm"].astype(str)
order, pal = ['control', 'treated'], ['#3B6DB3', '#C1432B']
fig, ax = plt.subplots(figsize=(5.6, 5.8))
kmfs = []
for i, grp in enumerate(order):
    sub = df[df["arm"] == grp]
    kmf = KaplanMeierFitter().fit(sub["time"], sub["event"], label=grp)
    kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i], show_censors=True, censor_styles={"marker": "|", "ms": 6})
    kmfs.append(kmf)
res = multivariate_logrank_test(df["time"], df["arm"], df["event"])
p = res.p_value
ax.text(0.04, 0.08, "log-rank p " + ("< 0.001" if p < 0.001 else "= %.2g" % p), transform=ax.transAxes, fontsize=11)
ax.set_xlabel("Months"); ax.set_ylabel("Survival probability")
ax.set_ylim(0, 1.02); ax.legend(title="arm", frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
add_at_risk_counts(*kmfs, ax=ax)
fig.tight_layout(); fig.savefig("figure_arm_km_py.pdf")
print("Reproduced figure.")
