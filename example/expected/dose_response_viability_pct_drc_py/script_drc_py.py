#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import curve_fit


def fpl(lx, bottom, top, log_ec50, hill):
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - lx) * hill))


raw = pd.read_csv("input_viability_pct.csv")
d = pd.DataFrame({"dose": pd.to_numeric(raw["concentration_um"], errors="coerce"),
                   "y": pd.to_numeric(raw["viability_pct"], errors="coerce"),
                   "group": raw["compound"].astype(str)})
d = d[np.isfinite(d["dose"]) & np.isfinite(d["y"]) & (d["dose"] > 0)]
d["lx"] = np.log10(d["dose"])
methods = {'cmpd_A': 'lm', 'cmpd_B': 'lm', 'cmpd_C': 'lm'}
starts = {
    "cmpd_A": [2.966666666666667, 98.86666666666667, -1.5228787452803376, -1.0],
    "cmpd_B": [11.4, 98.89999999999999, -0.5228787452803376, -1.0],
    "cmpd_C": [23.133333333333336, 102.3, 0.0, -1.0]
}
fig, ax = plt.subplots(figsize=(6.2, 4.2))
for g, p0 in starts.items():
    s = d[d["group"] == g]
    popt, pcov = curve_fit(fpl, s["lx"], s["y"], p0=p0, method=methods[g], maxfev=20000,
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    print(g, "bottom, top, log10 EC50, hill:", popt)
    grid = np.linspace(s["lx"].min(), s["lx"].max(), 200)
    ax.plot(10 ** grid, fpl(grid, *popt), label=g)
    ax.scatter(s["dose"], s["y"], s=14, alpha=0.6)
# A curve whose midpoint is not determined is drawn as points only, as in the figure.
for g in d["group"].unique():
    if g not in starts:
        s = d[d["group"] == g]
        ax.scatter(s["dose"], s["y"], s=14, alpha=0.6, label=f"{g} (not determined)")
ax.set_xscale("log")
ax.set_xlabel("Concentration (uM)"); ax.set_ylabel("Viability (%)")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figure_viability_pct_drc_py.pdf")
print("Reproduced figure.")
