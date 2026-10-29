"""Recipe: Bland-Altman agreement between two methods. Mirrors the R recipe.

Input is one row per subject with the reading of method A in --x and of method B
in --y. Each difference is A - B and each mean is (A + B) / 2. The figure plots
the difference against the mean with the bias and the 95% limits of agreement
(bias +/- 1.96 SD), each with a shaded 95% confidence interval. The bias CI is
t-based, and the limits-of-agreement CI uses the standard error of Bland and
Altman (1999), SD * sqrt(1/n + 1.96^2 / (2 (n - 1))). A regression of the
difference on the mean tests for a proportional bias.
"""
import functools

import numpy as np
import pandas as pd
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import shapiro_normal


def agreement_stats(a, b):
    """Compute the bias, the limits of agreement and their confidence intervals.

    Args:
        a: The readings of method A.
        b: The readings of method B, paired with ``a``.

    Returns:
        dict: The statistics, with the same keys as the R engine.
    """
    d = a - b
    m = (a + b) / 2
    n = len(d)
    bias, s = float(d.mean()), float(d.std(ddof=1))
    z, tq = stats.norm.ppf(0.975), stats.t.ppf(0.975, n - 1)
    se_bias = s / np.sqrt(n)
    se_loa = s * np.sqrt(1 / n + z ** 2 / (2 * (n - 1)))
    reg = stats.linregress(m, d)
    return {
        "n": n, "bias": bias,
        "bias_ci_low": bias - tq * se_bias, "bias_ci_high": bias + tq * se_bias,
        "sd_diff": s, "loa_low": bias - z * s, "loa_high": bias + z * s,
        "loa_low_ci_low": bias - z * s - tq * se_loa, "loa_low_ci_high": bias - z * s + tq * se_loa,
        "loa_high_ci_low": bias + z * s - tq * se_loa, "loa_high_ci_high": bias + z * s + tq * se_loa,
        "prop_bias_slope": float(reg.slope), "prop_bias_p": float(reg.pvalue),
        "shapiro_p_diff": shapiro_normal(d)[0],
    }


def recipe_bland_altman(df, spec):
    """Draw a Bland-Altman plot with the bias and the limits of agreement.

    Args:
        df: The tidy table, one row per subject.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing or fewer than 3 pairs are complete.
    """
    x, y = spec["data"].get("x"), spec["data"].get("y")
    if not x or not y:
        raise ValueError("bland_altman needs --x (method A) and --y (method B), one row per subject")
    for c in (x, y):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")
    n0 = len(df)
    a = pd.to_numeric(df[x], errors="coerce")
    b = pd.to_numeric(df[y], errors="coerce")
    keep = np.isfinite(a) & np.isfinite(b)
    a, b = a[keep].to_numpy(float), b[keep].to_numpy(float)
    if len(a) < 3:
        raise ValueError("bland_altman needs at least 3 complete pairs")
    dropped = int((~keep).sum())
    steps = ([f"Dropped {dropped} row(s) with a missing or non-numeric reading in '{x}' or "
              f"'{y}' ({n0} -> {len(a)})."] if dropped else ["No cleaning applied; input used as-is."])

    s = agreement_stats(a, b)
    mean, diff = (a + b) / 2, a - b
    ap = spec["appearance"]
    pal = palette(1, ap.get("palette"))
    w, h = 6.4, 4.4
    fig, ax = new_fig(w, h)
    for lo, hi in ((s["loa_low_ci_low"], s["loa_low_ci_high"]),
                   (s["bias_ci_low"], s["bias_ci_high"]),
                   (s["loa_high_ci_low"], s["loa_high_ci_high"])):
        ax.axhspan(lo, hi, color="0.8", alpha=0.5, lw=0)
    for val, ls, lab in ((s["loa_low"], "--", f"-1.96 SD: {s['loa_low']:.3g}"),
                         (s["bias"], "-", f"Bias: {s['bias']:.3g}"),
                         (s["loa_high"], "--", f"+1.96 SD: {s['loa_high']:.3g}")):
        ax.axhline(val, ls=ls, color="0.15", lw=1.0)
        # The label sits in the right margin, next to its line, so that it never
        # covers a point; tight_layout leaves room for it, as R does.
        ax.annotate(lab, xy=(1.0, val), xycoords=("axes fraction", "data"), xytext=(5, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=9,
                    color="0.15", annotation_clip=False)
    ax.scatter(mean, diff, s=26, color=pal[0], alpha=0.85, zorder=3)
    ax.set_title(f"n = {s['n']}; proportional bias slope {s['prop_bias_slope']:.3f}, "
                 f"P = {s['prop_bias_p']:.2g}", fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or f"Mean of {x} and {y}")
    ax.set_ylabel(ap.get("y_label") or f"Difference ({x} - {y})")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    stats_df = pd.DataFrame([{"method_a": x, "method_b": y, **s}])
    return {
        "fig": fig, "stats": stats_df,
        "test_meta": {"name": "Bland-Altman agreement analysis", **s},
        "methods": _methods(s, x, y),
        "df_used": pd.DataFrame({"a": a, "b": b, "mean": mean, "difference": diff}),
        "clean_steps": steps, "label": slugify(f"{x}_vs_{y}"), "width": w, "height": h,
        "qc": {"quantity": "differences between the methods", "panels": {f"{x} - {y}": diff}},
        "build_script": functools.partial(_script, spec),
    }


def _methods(s, x, y):
    """Write the methods paragraph.

    Args:
        s: The agreement statistics.
        x: The method A column.
        y: The method B column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    return (
        f"Agreement between {x} and {y} was assessed by the method of Bland and Altman "
        f"(n = {s['n']} pairs). The mean difference ({x} - {y}) was {s['bias']:.3g} "
        f"(95% CI {s['bias_ci_low']:.3g} to {s['bias_ci_high']:.3g}) with a standard deviation "
        f"of {s['sd_diff']:.3g}. The 95% limits of agreement were {s['loa_low']:.3g} (95% CI "
        f"{s['loa_low_ci_low']:.3g} to {s['loa_low_ci_high']:.3g}) and {s['loa_high']:.3g} "
        f"(95% CI {s['loa_high_ci_low']:.3g} to {s['loa_high_ci_high']:.3g}), computed as the "
        "mean difference plus or minus 1.96 standard deviations, with confidence intervals from "
        "the standard error of Bland and Altman (1999). A linear regression of the difference on "
        f"the mean gave a slope of {s['prop_bias_slope']:.3f} (P = {s['prop_bias_p']:.2g}) as a "
        "test of proportional bias. The normality of the differences, which the limits assume, "
        f"was checked with the Shapiro-Wilk test (P = {s['shapiro_p_diff']:.2g}). Computed in "
        f"Python {platform.python_version()} with scipy {v('scipy')} and drawn with matplotlib "
        f"{v('matplotlib')}."
    )


def _script(spec, in_name, fig_stub):
    """Write the standalone script that redraws the plot.

    Args:
        spec: The figure specification.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    x, y = spec["data"]["x"], spec["data"]["y"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("{in_name}")
a = pd.to_numeric(raw["{x}"], errors="coerce"); b = pd.to_numeric(raw["{y}"], errors="coerce")
keep = np.isfinite(a) & np.isfinite(b)
a, b = a[keep].to_numpy(), b[keep].to_numpy()
mean, diff = (a + b) / 2, a - b
bias, s = diff.mean(), diff.std(ddof=1)
loa = bias + np.array([-1, 1]) * stats.norm.ppf(0.975) * s
print(f"bias {{bias:.4f}}, SD {{s:.4f}}, limits {{loa[0]:.4f}} to {{loa[1]:.4f}}")
print(stats.linregress(mean, diff))
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.scatter(mean, diff)
ax.axhline(bias, color="k"); [ax.axhline(v, color="k", ls="--") for v in loa]
ax.set_xlabel("Mean of {x} and {y}"); ax.set_ylabel("Difference ({x} - {y})")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
