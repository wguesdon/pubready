"""Recipe: meta-analysis forest plot. Mirrors the R recipe, which uses metafor.

Input is one row per study with an effect estimate (--y), its 95% CI (--ci_low,
--ci_high) or its standard error on the analysis scale (--se), and a study label
(--label). A ratio measure is pooled on the log scale and a difference measure
on its own scale; --measure auto picks ratio when every value is positive and
the CIs are more symmetric on the log scale. The studies are pooled by inverse
variance with statsmodels combine_effects: a fixed-effect model and a
DerSimonian-Laird random-effects model, with Cochran's Q, I^2 and tau^2.
"""
import functools

import numpy as np
import pandas as pd
from matplotlib.patches import Polygon
from scipy import stats
from statsmodels.stats.meta_analysis import combine_effects

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette, plain_log_axis

Z975 = stats.norm.ppf(0.975)


def asymmetry(est, lo, hi):
    """Median relative asymmetry of the CIs around the estimates.

    Args:
        est: The estimates.
        lo: The lower bounds.
        hi: The upper bounds.

    Returns:
        float: The median of |(hi - est) - (est - lo)| / (hi - lo).
    """
    return float(np.median(np.abs((hi - est) - (est - lo)) / (hi - lo)))


def resolve_measure(measure, est, lo, hi):
    """Pick the ratio or the difference scale, as the R engine does.

    Args:
        measure: The requested measure: auto, ratio or difference.
        est: The estimates.
        lo: The lower bounds, or None when only standard errors are given.
        hi: The upper bounds, or None.

    Returns:
        str: "ratio" or "difference".

    Raises:
        ValueError: If the measure name is unknown.
    """
    measure = (measure or "auto").lower()
    if measure in ("ratio", "difference"):
        return measure
    if measure != "auto":
        raise ValueError(f"unknown --measure '{measure}' (use ratio or difference)")
    if lo is None or np.any(np.concatenate([est, lo, hi]) <= 0):
        return "difference"
    if asymmetry(np.log(est), np.log(lo), np.log(hi)) < asymmetry(est, lo, hi):
        return "ratio"
    return "difference"


def back(v, measure):
    """Return a value from the analysis scale to the reported scale.

    Args:
        v: The value on the analysis scale.
        measure: "ratio" or "difference".

    Returns:
        The value on the reported scale.
    """
    return np.exp(v) if measure == "ratio" else v


def recipe_forest_meta(df, spec):
    """Pool study estimates and draw a forest plot.

    Args:
        df: The tidy table, one row per study.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if fewer than 2 studies remain, or
            if a ratio measure holds a value at or below zero.
    """
    dat = spec["data"]
    y, lab = dat.get("y"), dat.get("label")
    lo_col, hi_col, se_col = dat.get("ci_low"), dat.get("ci_high"), dat.get("se")
    if not y:
        raise ValueError("forest_meta needs --y (the effect estimate of each study)")
    has_ci = bool(lo_col) and bool(hi_col)
    if not has_ci and not se_col:
        raise ValueError("forest_meta needs --ci_low and --ci_high (the 95% CI), or --se (the standard error)")
    cols = [y] + ([lo_col, hi_col] if has_ci else [se_col]) + ([lab] if lab else [])
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    est = pd.to_numeric(df[y], errors="coerce").to_numpy(float)
    labels = (df[lab].astype(str).to_numpy() if lab
              else np.array([f"Study {i + 1}" for i in range(n0)]))
    if has_ci:
        lo = pd.to_numeric(df[lo_col], errors="coerce").to_numpy(float)
        hi = pd.to_numeric(df[hi_col], errors="coerce").to_numpy(float)
        keep = np.isfinite(est) & np.isfinite(lo) & np.isfinite(hi)
        lo, hi = lo[keep], hi[keep]
        se_in = None
    else:
        se_in = pd.to_numeric(df[se_col], errors="coerce").to_numpy(float)
        keep = np.isfinite(est) & np.isfinite(se_in) & (se_in > 0)
        se_in = se_in[keep]
        lo = hi = None
    est, labels = est[keep], labels[keep]
    if len(est) < 2:
        raise ValueError("forest_meta needs at least 2 studies")

    measure = resolve_measure(dat.get("measure"), est, lo, hi)
    if measure == "ratio":
        vals = np.concatenate([est] + ([lo, hi] if has_ci else []))
        if np.any(vals <= 0):
            raise ValueError("a ratio measure needs every estimate and CI bound above zero")
    tr = np.log if measure == "ratio" else np.asarray
    yi = tr(est)
    sei = (tr(hi) - tr(lo)) / (2 * Z975) if has_ci else se_in
    if not has_ci:
        lo, hi = back(yi - Z975 * sei, measure), back(yi + Z975 * sei, measure)
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} study row(s) with a missing estimate or interval.")
    steps.append(f"Pooled {len(est)} of {n0} study row(s) as a {measure} measure"
                 f"{' on the log scale' if measure == 'ratio' else ''}.")
    if has_ci:
        steps.append("Derived each standard error from the 95% CI as (upper - lower) / (2 x 1.96) "
                     "on the analysis scale.")

    res = combine_effects(yi, sei ** 2, method_re="dl")
    sf = res.summary_frame()
    k = len(yi)
    q_p = float(stats.chi2.sf(res.q, k - 1))
    w_fe = 1 / sei ** 2
    w_re = 1 / (sei ** 2 + res.tau2)
    pooled_rows = []
    for name, key, wf, wr in (("Fixed effect", "fixed effect", 100.0, np.nan),
                              ("Random effects (DerSimonian-Laird)", "random effect", np.nan, 100.0)):
        eff, sd_eff = float(sf.loc[key, "eff"]), float(sf.loc[key, "sd_eff"])
        zval = eff / sd_eff
        pooled_rows.append({
            "row_type": "pooled", "label": name, "estimate": back(eff, measure),
            "ci_low": back(float(sf.loc[key, "ci_low"]), measure),
            "ci_high": back(float(sf.loc[key, "ci_upp"]), measure),
            "yi": eff, "se": sd_eff, "weight_fixed_pct": wf, "weight_random_pct": wr,
            "z": zval, "p_value": float(2 * stats.norm.sf(abs(zval))),
            "q": float(res.q), "q_df": k - 1, "q_p": q_p, "i2": 100 * float(res.i2),
            "tau2": float(res.tau2)})
    studies = pd.DataFrame({
        "row_type": "study", "label": labels, "estimate": est, "ci_low": lo, "ci_high": hi,
        "yi": yi, "se": sei, "weight_fixed_pct": 100 * w_fe / w_fe.sum(),
        "weight_random_pct": 100 * w_re / w_re.sum(), "z": np.nan, "p_value": np.nan,
        "q": np.nan, "q_df": np.nan, "q_p": np.nan, "i2": np.nan, "tau2": np.nan})
    pooled = pd.DataFrame(pooled_rows)
    stats_df = pd.concat([studies, pooled], ignore_index=True)

    fig, w, h = _plot(studies, pooled, measure, spec)
    test_meta = {
        "name": "Inverse-variance meta-analysis (fixed effect and DerSimonian-Laird random effects)",
        "measure": measure, "k": k,
        "fixed": {c: pooled_rows[0][c] for c in ("estimate", "ci_low", "ci_high", "p_value")},
        "random": {c: pooled_rows[1][c] for c in ("estimate", "ci_low", "ci_high", "p_value")},
        "heterogeneity": {"q": float(res.q), "df": k - 1, "p_value": q_p,
                          "i2": 100 * float(res.i2), "tau2": float(res.tau2)},
    }
    df_used = pd.DataFrame({"label": labels, "estimate": est, "ci_low": lo, "ci_high": hi,
                            "yi": yi, "se": sei})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(pooled, measure, k, y),
        "df_used": df_used, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, measure, has_ci),
    }


def _plot(studies, pooled, measure, spec):
    """Draw the forest plot with the estimate and weight columns on the right.

    Args:
        studies: The study rows.
        pooled: The two pooled rows.
        measure: "ratio" or "difference".
        spec: The figure specification.

    Returns:
        tuple: The figure, its width and its height.
    """
    ap = spec["appearance"]
    pal = palette(2, ap.get("palette"))
    k = len(studies)
    sy = np.arange(k, 0, -1) + 2.0
    py = np.array([1.5, 0.5])
    w, h = 8.2, 1.6 + 0.3 * (k + 3)
    fig, ax = new_fig(w, h)
    ax.axvline(1 if measure == "ratio" else 0, ls="--", color="0.55", lw=0.8)
    ax.axvline(pooled["estimate"].iloc[1], ls=":", color=pal[1], lw=0.8)
    ax.hlines(sy, studies["ci_low"], studies["ci_high"], color="0.15", lw=1.0)
    sizes = 12 + 140 * studies["weight_random_pct"] / studies["weight_random_pct"].max()
    ax.scatter(studies["estimate"], sy, marker="s", s=sizes, color=pal[0], zorder=3)
    for i in range(2):
        r = pooled.iloc[i]
        verts = [(r["ci_low"], py[i]), (r["estimate"], py[i] + 0.32),
                 (r["ci_high"], py[i]), (r["estimate"], py[i] - 0.32)]
        ax.add_patch(Polygon(verts, closed=True, facecolor=pal[1], edgecolor="0.15", lw=0.6))
    if measure == "ratio":
        ax.set_xscale("log")
        plain_log_axis(ax.xaxis)
    ax.set_yticks(list(sy) + list(py))
    ax.set_yticklabels(list(studies["label"]) + ["Fixed effect", "Random effects (DL)"], fontsize=10)
    ax.set_ylim(0, sy.max() + 0.9)
    ax.tick_params(axis="y", length=0)
    trans = ax.get_yaxis_transform()
    for yy, r in zip(list(sy) + list(py), list(studies.itertuples()) + list(pooled.itertuples())):
        ax.text(1.04, yy, f"{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]", transform=trans,
                va="center", ha="left", fontsize=9.5)
    for yy, wt in zip(sy, studies["weight_random_pct"]):
        ax.text(1.78, yy, f"{wt:.1f}%", transform=trans, va="center", ha="right", fontsize=9.5)
    ax.text(1.04, sy.max() + 0.8, "Estimate [95% CI]", transform=trans, va="center", ha="left",
            fontsize=9.5, fontweight="bold")
    ax.text(1.78, sy.max() + 0.8, "Weight", transform=trans, va="center", ha="right",
            fontsize=9.5, fontweight="bold")
    p = pooled.iloc[1]
    ax.set_title(f"Heterogeneity: Q({int(p['q_df'])}) = {p['q']:.1f}, P = {p['q_p']:.2g}; "
                 f"I2 = {p['i2']:.1f}%; tau2 = {p['tau2']:.3g}", fontsize=9.5, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or
                  f"{spec['data']['y']} (95% CI{', log scale' if measure == 'ratio' else ''})")
    apply_pub_style(ax, ap.get("theme"))
    ax.spines["left"].set_visible(False)
    fig.subplots_adjust(left=0.27, right=0.62, top=1 - 0.6 / h, bottom=0.9 / h)
    return fig, w, h


def _methods(pooled, measure, k, y):
    """Write the methods paragraph.

    Args:
        pooled: The two pooled rows.
        measure: "ratio" or "difference".
        k: The count of studies.
        y: The estimate column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    fe, re = pooled.iloc[0], pooled.iloc[1]
    scale = "on the log scale and back-transformed" if measure == "ratio" else "on the original scale"
    return (
        f"{k} studies were pooled by inverse-variance weighting {scale}, with the standard error "
        "of each study derived from its 95% confidence interval or given directly. The "
        f"fixed-effect estimate of {y} was {fe['estimate']:.3f} (95% CI {fe['ci_low']:.3f} to "
        f"{fe['ci_high']:.3f}, P = {fe['p_value']:.2g}), and the DerSimonian-Laird random-effects "
        f"estimate was {re['estimate']:.3f} (95% CI {re['ci_low']:.3f} to {re['ci_high']:.3f}, "
        f"P = {re['p_value']:.2g}). Heterogeneity was assessed with Cochran's Q (Q = {re['q']:.2f}, "
        f"df = {int(re['q_df'])}, P = {re['q_p']:.2g}), I2 = {re['i2']:.1f}% and tau2 = "
        f"{re['tau2']:.4f}. Squares are sized by the random-effects weight and the diamonds show "
        f"the pooled estimates. Computed in Python {platform.python_version()} with statsmodels "
        f"{v('statsmodels')} (combine_effects) and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, measure, has_ci, in_name, fig_stub):
    """Write the standalone script that repeats the pooling.

    Args:
        spec: The figure specification.
        measure: "ratio" or "difference".
        has_ci: Whether the input gives a CI rather than a standard error.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    d = spec["data"]
    tr = "np.log" if measure == "ratio" else "np.asarray"
    se_line = (f'sei = ({tr}(raw["{d["ci_high"]}"]) - {tr}(raw["{d["ci_low"]}"])) / (2 * 1.959963984540054)'
               if has_ci else f'sei = raw["{d["se"]}"].to_numpy()')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd
from statsmodels.stats.meta_analysis import combine_effects

raw = pd.read_csv("{in_name}")
yi = {tr}(raw["{d["y"]}"].to_numpy())
{se_line}
res = combine_effects(yi, np.asarray(sei) ** 2, method_re="dl")
print(res.summary_frame())
print("Q", res.q, "I2", res.i2, "tau2", res.tau2)
fig = res.plot_forest()
fig.savefig("{fig_stub}.pdf", bbox_inches="tight")
print("Reproduced figure.")
'''
