"""Normality QC for the Python engine. Mirrors core/r/R/qc.R: a QQ panel per
quantity a recipe's normality check runs on, labelled with the Shapiro-Wilk p,
skewness, and excess kurtosis. Used by the bundle writer and by figkit diagnose,
so the test choice is shown to the scientist, not just asserted. The QQ plot is
the tiebreaker the Shapiro p-value cannot be alone: Shapiro is underpowered at
small n and over-rejects at large n."""
import math

import matplotlib
matplotlib.use("Agg")  # headless container
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from .theme import PALETTE, apply_pub_style
from .util import shapiro_normal


def qc_moment_stats(v):
    """Sample skewness and excess kurtosis on the finite values, plus Shapiro p."""
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    n = int(len(v))
    out = {"n": n, "shapiro_p": float("nan"), "skewness": float("nan"), "kurtosis": float("nan")}
    if n < 3:
        return out
    out["shapiro_p"] = shapiro_normal(v)[0]
    s = float(np.std(v, ddof=1))
    if math.isfinite(s) and s > 0:
        z = (v - v.mean()) / s
        out["skewness"] = float(np.sum(z ** 3) / n)
        out["kurtosis"] = float(np.sum(z ** 4) / n - 3)
    return out


def qc_dims(n_panels):
    ncol = min(max(n_panels, 1), 3)
    nrow = math.ceil(n_panels / ncol)
    return (2.9 * ncol + 0.4, 2.9 * nrow + 1.0)


def _fmt(x, g=False):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{x:.2g}" if g else f"{x:.2f}"


def qc_normality(qc):
    """qc = {"quantity": label, "panels": {name: array}}.
    Returns (fig, panels) where panels is a list of per-panel stat dicts."""
    names = list(qc["panels"].keys())
    panels = []
    for k in names:
        s = qc_moment_stats(qc["panels"][k])
        s["panel"] = k
        panels.append(s)

    n = len(names)
    ncol = min(max(n, 1), 3)
    nrow = math.ceil(n / ncol)
    w, h = qc_dims(n)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), squeeze=False)
    flat = [ax for row in axes for ax in row]

    for ax, s in zip(flat, panels):
        v = np.asarray(qc["panels"][s["panel"]], dtype=float)
        v = v[np.isfinite(v)]
        if len(v) >= 2:
            (osm, osr), (slope, intercept, _) = stats.probplot(v, dist="norm")
            ax.scatter(osm, osr, s=14, color=PALETTE[0], alpha=0.75, edgecolors="none")
            ax.plot(osm, slope * osm + intercept, color="#333333", linewidth=1.0)
        ax.set_title(f"{s['panel']}  (n = {s['n']})\n"
                     f"p = {_fmt(s['shapiro_p'], True)}  skew = {_fmt(s['skewness'])}  "
                     f"kurt = {_fmt(s['kurtosis'])}",
                     fontsize=8)
        ax.set_xlabel("Theoretical quantiles", fontsize=9)
        ax.set_ylabel("Sample quantiles", fontsize=9)
        apply_pub_style(ax)
    for ax in flat[n:]:
        ax.set_visible(False)

    fig.suptitle(f"Normality QC: {qc['quantity']}", fontweight="bold", fontsize=12, y=1.0)
    fig.text(0.5, 0.925,
             "Points near the line support normality; mild, judged deviation is usually acceptable.",
             ha="center", fontsize=8, color="#4d4d4d")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    return fig, panels
