#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import curve_fit


def fpl(lx, bottom, top, log_ec50, hill):
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - lx) * hill))


raw = pd.read_csv("input_absorbance.csv")
d = pd.DataFrame({"conc": pd.to_numeric(raw["conc_ng_ml"], errors="coerce"),
                   "sig": pd.to_numeric(raw["absorbance"], errors="coerce")})
std = d[np.isfinite(d["conc"]) & np.isfinite(d["sig"])]
unk = d[d["conc"].isna() & np.isfinite(d["sig"])].copy()
pos = std[std["conc"] > 0]
p0 = [0.1565, 2.713, 0.3979400086720376, 1.0]
popt, _ = curve_fit(fpl, np.log10(pos["conc"]), pos["sig"], p0=p0, method="lm",
                    maxfev=20000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
b, t, le, h = popt
if b > t:
    b, t, h = t, b, -h


def f(c):
    return fpl(np.log10(c), b, t, le, h)


def finv(s):
    return 10 ** (le - np.log10((t - b) / (s - b) - 1) / h)


rng = (pos["conc"].min(), pos["conc"].max())
print("bottom, top, EC50, hill:", b, t, 10 ** le, h)
ends = [f(rng[0]), f(rng[1])]
inside = (unk["sig"] >= min(ends)) & (unk["sig"] <= max(ends))
unk["concentration"] = np.where(inside, finv(unk["sig"].where(inside)), np.nan)
print(unk)
fig, ax = plt.subplots(figsize=(5.2, 4.8))
grid = np.logspace(np.log10(rng[0]), np.log10(rng[1]), 200)
ax.plot(grid, f(grid)); ax.scatter(std["conc"], std["sig"])
ax.scatter(unk["concentration"], unk["sig"], marker="D", facecolors="white", edgecolors="red")
ax.set_xscale('log')
ax.set_xlabel("IL-6 (ng/mL)"); ax.set_ylabel("Absorbance at 450 nm")
fig.tight_layout(); fig.savefig("figure_absorbance_stdcurve_py.pdf")
print("Reproduced figure.")
