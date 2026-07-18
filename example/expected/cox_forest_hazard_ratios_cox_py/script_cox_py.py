#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv("input_hazard_ratios.csv")
covs = ['arm', 'age', 'sex', 'stage']
design = df[["time", "event"]].copy()
for c in covs:
    s = df[c]
    if s.dtype == object or not np.issubdtype(s.dtype, np.number):
        for lvl in sorted(s.astype(str).unique())[1:]:
            design["%s[%s]" % (c, lvl)] = (s.astype(str) == lvl).astype(int)
    else:
        design[c] = pd.to_numeric(s, errors="coerce")
cph = CoxPHFitter().fit(design, "time", "event")
s = cph.summary
terms = list(s.index); hr = s["exp(coef)"].to_numpy()
lo = s["exp(coef) lower 95%"].to_numpy(); hi = s["exp(coef) upper 95%"].to_numpy(); pv = s["p"].to_numpy()
n = len(terms); yy = np.arange(n)[::-1]
fig, ax = plt.subplots(figsize=(7.4, max(3.6, 0.5 * n + 2.0)))
ax.errorbar(hr, yy, xerr=[hr - lo, hi - hr], fmt="s", color="black", capsize=3, ms=6, lw=1.0)
ax.axvline(1.0, ls=":", color="grey"); ax.set_xscale("log")
ax.set_yticks(yy); ax.set_yticklabels(terms); ax.set_ylim(-0.6, n - 0.4)
ax.set_xlabel("Hazard ratio (95% CI)"); ax.set_title("Hazard ratio", fontweight="bold")
for y_, h_, l_, u_, p_ in zip(yy, hr, lo, hi, pv):
    ax.annotate("%.2f (%.2f-%.2f)  P = %.2g" % (h_, l_, u_, p_), xy=(1.02, y_), xycoords=("axes fraction", "data"), va="center", fontsize=8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("figure_hazard_ratios_cox_py.pdf")
print("Reproduced figure.")
