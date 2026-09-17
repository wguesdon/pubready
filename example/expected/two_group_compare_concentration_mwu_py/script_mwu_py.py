#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg
from scipy import stats as ss
from statannotations.Annotator import Annotator

df = pd.read_csv("input_concentration.csv")
df["group"] = df["group"].astype(str)
order = ['control', 'treated']
a = df.loc[df["group"] == order[0], "concentration"].to_numpy()
b = df.loc[df["group"] == order[1], "concentration"].to_numpy()
exact = len(a) < 50 and len(b) < 50 and len(set(list(a) + list(b))) == len(a) + len(b)
res = (ss.mannwhitneyu(a, b, alternative="two-sided", method="exact") if exact
       else ss.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic",
                            use_continuity=True))
p = float(res.pvalue)

fig, ax = plt.subplots(figsize=(3.8, 4.0))
sns.boxplot(data=df, x="group", y="concentration", hue="group", order=order, hue_order=order, palette=['#3B6DB3', '#C1432B'], ax=ax, legend=False, width=0.6, fliersize=0)
sns.stripplot(data=df, x="group", y="concentration", hue="group", order=order, hue_order=order, palette=['#3B6DB3', '#C1432B'],
              ax=ax, legend=False, size=4, alpha=0.75, edgecolor="black", linewidth=0.3, jitter=0.12)
annot = Annotator(ax, [(order[0], order[1])], data=df, x="group", y="concentration", order=order)
annot.configure(line_width=1.0, fontsize=12)
annot.set_custom_annotations(["***"])
annot.annotate()
ax.set_xlabel("Group"); ax.set_ylabel("IL-6 (pg/mL)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figure_concentration_mwu_py.pdf")
print("Reproduced figure.")
