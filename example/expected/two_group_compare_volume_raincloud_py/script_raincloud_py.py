#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd, seaborn as sns, pingouin as pg
from scipy import stats as ss
from statannotations.Annotator import Annotator


def label(v):
    """Write a group value as the figure does: a whole number with no decimal point."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


df = pd.read_csv("input_volume.csv")
df["group"] = df["group"].map(label)
order = ['control', 'treated']
a = df.loc[df["group"] == order[0], "volume"].to_numpy()
b = df.loc[df["group"] == order[1], "volume"].to_numpy()
res = pg.ttest(a, b, paired=False, correction=True); p = float(res["p_val"].iloc[0])

fig, ax = plt.subplots(figsize=(3.8, 4.0))
# A half violin, a narrow box and the rain to the left of the box.
pal = ['#3B6DB3', '#C1432B']
for i, g in enumerate(order):
    vals = df.loc[df["group"] == g, "volume"].dropna().to_numpy()
    if len(vals) < 2:
        continue
    v = ax.violinplot([vals], positions=[i], showextrema=False, widths=0.6)
    for body in v["bodies"]:
        verts = body.get_paths()[0].vertices
        verts[:, 0] = np.clip(verts[:, 0], i, np.inf)  # keep the right half
        body.set_facecolor(pal[i]); body.set_edgecolor("none"); body.set_alpha(0.85)
    bp = ax.boxplot([vals], positions=[i], widths=0.1, patch_artist=True,
                    showfliers=False, manage_ticks=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(pal[i]); patch.set_alpha(0.9)
    for med in bp["medians"]:
        med.set_color("black")
    k = np.arange(1, len(vals) + 1)
    jit = -0.18 + 0.12 * np.mod(k * 0.6180339887498949, 1.0)
    ax.scatter(np.full(len(vals), i) + jit, vals, s=10, color=pal[i], alpha=0.6,
               edgecolors="black", linewidths=0.2, zorder=3)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
annot = Annotator(ax, [(order[0], order[1])], data=df, x="group", y="volume", order=order)
annot.configure(line_width=1.0, fontsize=12)
annot.set_custom_annotations(["****"])
annot.annotate()
ax.set_xlabel("group"); ax.set_ylabel("Tumor volume (mm^3)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("figure_volume_raincloud_py.pdf")
print("Reproduced figure.")
