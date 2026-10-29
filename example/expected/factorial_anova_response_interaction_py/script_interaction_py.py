#!/usr/bin/env python3
# Standalone reproduction (parametric ANOVA). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import numpy as np
import pandas as pd, seaborn as sns
import statsmodels.api as sm, statsmodels.formula.api as smf

df = pd.read_csv("input_response.csv")
factors = ['genotype', 'treatment']
for f in factors: df[f] = df[f].astype(str)
ren = {"response": "_y_"}; ren.update({f: "_f%d_" % i for i, f in enumerate(factors)})
d2 = df.rename(columns=ren)
rhs = "*".join("C(_f%d_)" % i for i in range(len(factors)))
aov = sm.stats.anova_lm(smf.ols("_y_ ~ " + rhs, data=d2).fit(), typ=2)
def clean(t):
    for i, f in enumerate(factors): t = t.replace("C(_f%d_)" % i, f)
    return t
def plab(p): return "P < 0.001" if p < 0.001 else "P = %.2g" % p
sig = [(clean(t), r["PR(>F)"]) for t, r in aov.iterrows() if t != "Residual" and r["PR(>F)"] < 0.05]
subtitle = "\n".join("%s: %s" % (t.replace(":", " × "), plab(p)) for t, p in sig) or "No significant effects (P ≥ 0.05)"

f1_lv, f2_lv = ['ko', 'wt'], ['high', 'low', 'vehicle']
facet_lv = [None]
fig, axes = plt.subplots(1, len(facet_lv), figsize=(5.2, 4.4), sharey=True, squeeze=False)
axes = axes[0]
for ax, fl in zip(axes, facet_lv):
    sub = df if fl is None else df[df["None"] == fl]
    k = len(f2_lv)
    for i, lv in enumerate(f2_lv):
        g = sub[sub["treatment"] == lv].groupby("genotype")["response"]
        m = g.mean().reindex(f1_lv); se = (g.std() / g.size() ** 0.5).reindex(f1_lv)
        x = [j + (i + 0.5 - k / 2) * 0.25 / k for j in range(len(f1_lv))]
        ax.errorbar(x, m, yerr=se, color=['#3B6DB3', '#C1432B', '#2E8B57'][i], marker="o", capsize=3, label=lv)
    ax.set_xticks(range(len(f1_lv))); ax.set_xticklabels(f1_lv)
    if fl is not None: ax.set_title(str(fl))
    ax.set_xlabel("genotype"); ax.set_ylabel("Response (a.u.)" if ax is axes[0] else "")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1].legend(title="treatment", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
fig.suptitle(subtitle, fontsize=9.5, color="#4d4d4d")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig("figure_response_interaction_py.pdf")
print("Reproduced figure.")
