#!/usr/bin/env python3
# Standalone reproduction (parametric ANOVA). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import numpy as np
import pandas as pd, seaborn as sns
import statsmodels.api as sm, statsmodels.formula.api as smf


def label(v):
    """Write a group value as the figure does: a whole number with no decimal point."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


df = pd.read_csv("input_response.csv")
factors = ['genotype', 'treatment', 'sex']
for f in factors: df[f] = df[f].map(label)
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

f1_lv, f2_lv = ['ko', 'wt'], ['drug', 'vehicle']
facet_lv = ['female', 'male']
fig, axes = plt.subplots(1, len(facet_lv), figsize=(6.8, 4.4), sharey=True, squeeze=False)
axes = axes[0]
for ax, fl in zip(axes, facet_lv):
    sub = df if fl is None else df[df["sex"] == fl]
    sns.boxplot(data=sub, x="genotype", y="response", hue="treatment", order=f1_lv, hue_order=f2_lv, palette=['#3B6DB3', '#C1432B'], ax=ax, width=0.7, fliersize=0, legend=False)
    # Each point sits at its dodged cell centre plus a fixed golden-ratio offset.
    k = len(f2_lv); a = 0.15 / (k + 2)
    for i, lv in enumerate(f2_lv):
        for j, v in enumerate(f1_lv):
            vals = sub.loc[(sub["treatment"] == lv) & (sub["genotype"] == v), "response"].to_numpy(dtype=float)
            n = np.arange(1, len(vals) + 1)
            x = j + (i + 0.5 - k / 2) * 0.7 / k - a + 2 * a * np.mod(n * 0.6180339887498949, 1.0)
            ax.scatter(x, vals, s=10, color=['#3B6DB3', '#C1432B'][i], alpha=0.7, edgecolors="black", linewidths=0.3, zorder=3)
    if fl is not None: ax.set_title(str(fl))
    ax.set_xlabel("genotype"); ax.set_ylabel("Response (a.u.)" if ax is axes[0] else "")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1].legend(handles=[mpatches.Patch(color=['#3B6DB3', '#C1432B'][i], label=l) for i, l in enumerate(f2_lv)], title="treatment", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
fig.suptitle(subtitle, fontsize=9.5, color="#4d4d4d")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig("figure_response_threeway_py.pdf")
print("Reproduced figure.")
