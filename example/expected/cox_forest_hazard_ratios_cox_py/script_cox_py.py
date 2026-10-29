#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from matplotlib import transforms
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter

df = pd.read_csv("input_hazard_ratios.csv")
covs = ['arm', 'age', 'sex', 'stage']
design = df[["time", "event"]].copy()
rows = []  # (label, term); term None marks the reference level of a factor
for c in covs:
    s = df[c]
    if s.dtype == object or not np.issubdtype(s.dtype, np.number):
        cats = sorted(s.astype(str).unique())
        rows.append(("%s: %s" % (c, cats[0]), None))
        for lvl in cats[1:]:
            design["%s%s" % (c, lvl)] = (s.astype(str) == lvl).astype(int)
            rows.append(("%s: %s" % (c, lvl), "%s%s" % (c, lvl)))
    else:
        design[c] = pd.to_numeric(s, errors="coerce")
        rows.append((c, c))
cph = CoxPHFitter().fit(design, "time", "event", fit_options={"precision": 1e-11})
s = cph.summary
print(s)

k = len(rows); yy = np.arange(k, 0, -1).astype(float)
ticks = [0.5, 1.0, 2.0]
w, h, lab_in = 6.8200, 4.3000, 1.3200
fig = plt.figure(figsize=(w, h))
ax = fig.add_axes([lab_in / w, 0.85 / h, 2.85 / w, (h - 1.3) / h])
axt = fig.add_axes([(lab_in + 3.0) / w, 0.85 / h, 2.3 / w, (h - 1.3) / h])
ax.axvline(1.0, ls=":", color="0.4", lw=0.9)
trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
for y_, (label, term) in zip(yy, rows):
    fig.text(0.15 / w, y_, label, transform=trans, ha="left", va="center", fontsize=10.5)
    if term is None:
        axt.text(0, y_, "reference", ha="left", va="center", fontsize=10)
        continue
    hr_, lo_, hi_, p_ = s.loc[term, ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    ax.hlines(y_, lo_, hi_, color="0.15", lw=1.0)
    ax.scatter([hr_], [y_], marker="s", s=48, color="#3B6DB3", zorder=3)
    star = "***" if p_ < 0.001 else "**" if p_ < 0.01 else "*" if p_ < 0.05 else ""
    p_txt = ("< 0.001" if p_ < 0.001 else "%.2g" % p_) + (" " + star if star else "")
    axt.text(0, y_, "%.2f (%.2f\u2013%.2f)" % (hr_, lo_, hi_), ha="left", va="center", fontsize=10)
    axt.text(0.66, y_, p_txt, ha="left", va="center", fontsize=10)
axt.text(0, k + 1, "HR (95% CI)", ha="left", va="center", fontsize=10, fontweight="bold")
axt.text(0.66, k + 1, "P", ha="left", va="center", fontsize=10, fontweight="bold")
ax.set_xscale("log")
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FuncFormatter("{:g}".format))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_ylim(0.4, k + 1.6); axt.set_ylim(0.4, k + 1.6); axt.set_xlim(0, 1); axt.axis("off")
ax.set_yticks([]); ax.set_xlabel('Hazard ratio (95% CI, log scale)', fontsize=13)
ax.set_title("%d subjects, %d events" % (len(df), int(df["event"].sum())), fontsize=9.5,
             color="#4d4d4d", loc="left")
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
fig.savefig("figure_hazard_ratios_cox_py.pdf", bbox_inches="tight")
print("Reproduced figure.")
