#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import curve_fit


def mm(s, vmax, km):
    return vmax * s / (km + s)


raw = pd.read_csv("input_rate.csv")
d = pd.DataFrame({"s": pd.to_numeric(raw["substrate_um"], errors="coerce"),
                   "v": pd.to_numeric(raw["rate"], errors="coerce"),
                   "group": raw["inhibitor"].astype(str)})
d = d[np.isfinite(d["s"]) & np.isfinite(d["v"]) & (d["s"] >= 0)]
methods = {'0 uM': 'lm', '2 uM': 'lm', '10 uM': 'lm'}
starts = {
    "0 uM": [114.31, 20.0],
    "2 uM": [111.75, 20.0],
    "10 uM": [105.73, 40.0]
}
fig, ax = plt.subplots(figsize=(6.6, 4.4))
grid = np.linspace(0, d["s"].max(), 200)
for g, p0 in starts.items():
    s = d[d["group"] == g]
    popt, pcov = curve_fit(mm, s["s"], s["v"], p0=p0, method=methods[g], maxfev=20000,
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    print(g, "Vmax, Km:", popt, "SE:", np.sqrt(np.diag(pcov)))
    ax.plot(grid, mm(grid, *popt), label=g)
    ax.scatter(s["s"], s["v"], s=14, alpha=0.6)
ax.set_xlim(left=0); ax.set_ylim(bottom=0)
ax.set_xlabel("Substrate (uM)"); ax.set_ylabel("Rate (nmol/min/mg)")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figure_rate_enzyme_py.pdf")
print("Reproduced figure.")
