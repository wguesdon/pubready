#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The ANOVA
# uses the subjects seen at every time point.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, pingouin as pg, seaborn as sns

raw = pd.read_csv("input_volume_mm3.csv").dropna(subset=["mouse", "day", "volume_mm3"])
n_times = raw["day"].nunique()
raw = raw[raw.groupby("mouse")["day"].transform("count") == n_times]
print(pg.mixed_anova(data=raw, dv="volume_mm3", within="day", subject="mouse", between="arm"))
fig, ax = plt.subplots(figsize=(6.0, 4.4))
sns.lineplot(data=raw, x="day", y="volume_mm3", hue="arm", errorbar="se", marker="o", ax=ax)
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("figure_volume_mm3_tc_py.pdf")
print("Reproduced figure.")
