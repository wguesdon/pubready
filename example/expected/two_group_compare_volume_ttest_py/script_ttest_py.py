#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg
from statannotations.Annotator import Annotator

df = pd.read_csv("input_volume.csv")
df["group"] = df["group"].astype(str)
order = ['control', 'treated']
a = df.loc[df["group"] == order[0], "volume"].to_numpy()
b = df.loc[df["group"] == order[1], "volume"].to_numpy()
res = pg.ttest(a, b, paired=False, correction=False); p = float(res["p_val"].iloc[0])

fig, ax = plt.subplots(figsize=(3.8, 4.0))
sns.boxplot(data=df, x="group", y="volume", hue="group", order=order, palette=['#3B6DB3', '#C1432B'], ax=ax, legend=False, width=0.6, fliersize=0)
sns.stripplot(data=df, x="group", y="volume", hue="group", order=order, palette=['#3B6DB3', '#C1432B'],
              ax=ax, legend=False, size=4, alpha=0.75, edgecolor="black", linewidth=0.3, jitter=0.12)
annot = Annotator(ax, [(order[0], order[1])], data=df, x="group", y="volume", order=order)
annot.configure(line_width=1.0, fontsize=12)
annot.set_custom_annotations(["****"])
annot.annotate()
ax.set_xlabel("Group"); ax.set_ylabel("Tumor volume (mm^3)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figure_volume_ttest_py.pdf")
print("Reproduced figure.")
