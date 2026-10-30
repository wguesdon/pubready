"""Normality QC for the Python engine. Mirrors core/r/R/qc.R: a QQ panel per
quantity a recipe's normality check runs on, labelled with the Shapiro-Wilk p,
skewness, and excess kurtosis. Used by the bundle writer and by figkit diagnose,
so the test choice is shown to the scientist, not just asserted. The QQ plot is
the tiebreaker the Shapiro p-value cannot be alone: Shapiro is underpowered at
small n and over-rejects at large n."""
import math
import textwrap

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
    lay = qc_layout(qc["quantity"], n)
    w, h = lay["width"], lay["height"]
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

    # Figure text, not a suptitle, so tight_layout keeps the header space of
    # the rect and does not reserve it a second time.
    n_title, n_sub = len(lay["title"]), len(lay["subtitle"])
    sub_top = 0.12 + 0.22 * n_title
    fig.text(0.5, 1 - 0.08 / h, "\n".join(lay["title"]), ha="center", va="top",
             fontweight="bold", fontsize=12)
    fig.text(0.5, 1 - sub_top / h, "\n".join(lay["subtitle"]),
             ha="center", va="top", fontsize=8, color="#4d4d4d")
    header_in = sub_top + 0.15 * n_sub + 0.1
    fig.tight_layout(rect=[0, 0, 1, 1 - header_in / h])
    return fig, panels


#: The subtitle of the QC figure, and the characters per inch of the title and
#: the subtitle. Both engines wrap the two lines with this rule to the width of
#: the figure less its axis margin, so a one-panel figure does not cut them.
QC_SUBTITLE = ("Points near the line support normality; mild, judged deviation is "
               "usually acceptable.")
QC_TITLE_CPI = 8.5
QC_SUBTITLE_CPI = 15
QC_AXIS_MARGIN_IN = 0.5


def qc_wrap(text, width_in, cpi):
    """Wrap one header line of the QC figure, as R strwrap() does.

    Args:
        text: The line.
        width_in: The width of the figure in inches.
        cpi: The characters per inch of the font.

    Returns:
        list: The lines, each of floor((width - margin) * cpi) characters or
        fewer.
    """
    width = int(math.floor((width_in - QC_AXIS_MARGIN_IN) * cpi))
    return textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False)


def qc_layout(quantity, n_panels):
    """The size of the QC figure and its wrapped title and subtitle.

    Each line past the first adds height, so the panels keep their size. The R
    engine (qc.R) has the same rule.

    Args:
        quantity: The checked quantity, for the title.
        n_panels: The number of QQ panels.

    Returns:
        dict: The width and the height in inches, and the title and subtitle
        lines.
    """
    w, h = qc_dims(n_panels)
    title = qc_wrap(f"Normality QC: {quantity}", w, QC_TITLE_CPI)
    subtitle = qc_wrap(QC_SUBTITLE, w, QC_SUBTITLE_CPI)
    h = h + 0.25 * (len(title) - 1) + 0.18 * (len(subtitle) - 1)
    return {"width": w, "height": h, "title": title, "subtitle": subtitle}
