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
from matplotlib.ticker import NullLocator
from scipy import stats
from statsmodels.stats.meta_analysis import combine_effects

from ..bundle import slugify
from ..theme import apply_pub_style, log_ticks_ref, new_fig, palette, plain_log_axis
from .estimation_plot import wrap

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


def dersimonian_laird(yi, vi):
    """DerSimonian-Laird random-effects pooling, with the rules of metafor.

    tau2 is the moment estimate (Q - df) / C truncated at zero, and I2 is
    tau2 / (tau2 + s2) with the typical within-study variance s2, so I2 is zero
    when Q is below its degrees of freedom. statsmodels combine_effects does not
    truncate, which gives a negative tau2, a negative I2 and negative weights
    for a homogeneous set of studies.

    Args:
        yi: The study estimates on the analysis scale.
        vi: Their sampling variances.

    Returns:
        dict: q, tau2, i2 (in percent), the random-effects estimate and its
        standard error.
    """
    yi = np.asarray(yi, dtype=float)
    vi = np.asarray(vi, dtype=float)
    w = 1.0 / vi
    k = len(yi)
    fe = float(np.sum(w * yi) / np.sum(w))
    q = float(np.sum(w * (yi - fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    s2 = (k - 1) * np.sum(w) / (np.sum(w) ** 2 - np.sum(w ** 2)) if c > 0 else 0.0
    i2 = 100.0 * tau2 / (tau2 + s2) if (tau2 + s2) > 0 else 0.0
    w_re = 1.0 / (vi + tau2)
    est = float(np.sum(w_re * yi) / np.sum(w_re))
    return {"q": q, "tau2": tau2, "i2": i2, "estimate": est,
            "se": float(np.sqrt(1.0 / np.sum(w_re)))}


def back(v, measure):
    """Return a value from the analysis scale to the reported scale.

    Args:
        v: The value on the analysis scale.
        measure: "ratio" or "difference".

    Returns:
        The value on the reported scale.
    """
    return np.exp(v) if measure == "ratio" else v


def study_labels(v, n):
    """The study labels, with "Study <row>" for an empty cell, as .fm_labels() in R.

    Args:
        v: The --label column, or None.
        n: The count of rows.

    Returns:
        numpy.ndarray: The labels.
    """
    out = []
    for i in range(n):
        s = None if v is None or pd.isna(v.iloc[i]) else str(v.iloc[i]).strip()
        out.append(s if s else f"Study {i + 1}")
    return np.array(out, dtype=object)


def _num(v):
    """Format a number as R format() does, with 7 significant digits.

    Args:
        v: The number.

    Returns:
        str: The text.
    """
    return f"{float(v):.7g}"


def study_problem(i, est, lo, hi, se):
    """Say why one study row cannot be pooled, as .fm_study_problem() in R (F77).

    Args:
        i: The row.
        est: The estimates.
        lo: The lower CI bounds, or None with --se.
        hi: The upper CI bounds, or None with --se.
        se: The standard errors, or None with a CI.

    Returns:
        str or None: The reason, or None when the row can be pooled.
    """
    if se is not None:
        return f"its standard error is {_num(se[i])}" if se[i] <= 0 else None
    if lo[i] > hi[i]:
        return f"its CI bounds are reversed ({_num(lo[i])} to {_num(hi[i])})"
    if lo[i] == hi[i]:
        return f"its CI has zero width ({_num(lo[i])} to {_num(hi[i])})"
    if est[i] < lo[i] or est[i] > hi[i]:
        return (f"its CI ({_num(lo[i])} to {_num(hi[i])}) does not hold its estimate "
                f"{_num(est[i])}")
    return None


def check_studies(labels, est, lo, hi, se):
    """Stop and name each study that cannot be pooled, as .fm_check_studies() in R.

    Args:
        labels: The study labels.
        est: The estimates.
        lo: The lower CI bounds, or None.
        hi: The upper CI bounds, or None.
        se: The standard errors, or None.

    Raises:
        ValueError: If a study row cannot be pooled.
    """
    bad = [(labels[i], why) for i in range(len(est))
           if (why := study_problem(i, est, lo, hi, se)) is not None]
    if not bad:
        return
    msg = "; ".join(f"study '{lab}': {why}" for lab, why in bad[:5])
    if len(bad) > 5:
        msg += f"; and {len(bad) - 5} more"
    raise ValueError(f"forest_meta cannot pool {len(bad)} study row(s): {msg}")


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
    labels = study_labels(df[lab] if lab else None, n0)
    if has_ci:
        lo = pd.to_numeric(df[lo_col], errors="coerce").to_numpy(float)
        hi = pd.to_numeric(df[hi_col], errors="coerce").to_numpy(float)
        keep = np.isfinite(est) & np.isfinite(lo) & np.isfinite(hi)
        lo, hi = lo[keep], hi[keep]
        se_in = None
    else:
        se_in = pd.to_numeric(df[se_col], errors="coerce").to_numpy(float)
        keep = np.isfinite(est) & np.isfinite(se_in)
        se_in = se_in[keep]
        lo = hi = None
    est, labels = est[keep], labels[keep]
    check_studies(labels, est, lo, hi, se_in)
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

    # The fixed effect comes from statsmodels. The random-effects model is the
    # DerSimonian-Laird estimate with tau2 truncated at zero, as metafor does;
    # combine_effects lets tau2 go negative when Q is below its df.
    res = combine_effects(yi, sei ** 2, method_re="dl")
    sf = res.summary_frame()
    dl = dersimonian_laird(yi, sei ** 2)
    k = len(yi)
    q_p = float(stats.chi2.sf(dl["q"], k - 1))
    w_fe = 1 / sei ** 2
    w_re = 1 / (sei ** 2 + dl["tau2"])
    fe_eff = float(sf.loc["fixed effect", "eff"])
    fe_se = float(sf.loc["fixed effect", "sd_eff"])
    pooled_rows = []
    for name, eff, sd_eff, wf, wr in (
            ("Fixed effect", fe_eff, fe_se, 100.0, np.nan),
            ("Random effects (DerSimonian-Laird)", dl["estimate"], dl["se"], np.nan, 100.0)):
        zval = eff / sd_eff
        pooled_rows.append({
            "row_type": "pooled", "label": name, "estimate": back(eff, measure),
            "ci_low": back(eff - Z975 * sd_eff, measure),
            "ci_high": back(eff + Z975 * sd_eff, measure),
            "yi": eff, "se": sd_eff, "weight_fixed_pct": wf, "weight_random_pct": wr,
            "z": zval, "p_value": float(2 * stats.norm.sf(abs(zval))),
            "q": dl["q"], "q_df": k - 1, "q_p": q_p, "i2": dl["i2"], "tau2": dl["tau2"]})
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
        "heterogeneity": {"q": dl["q"], "df": k - 1, "p_value": q_p,
                          "i2": dl["i2"], "tau2": dl["tau2"]},
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


#: The width of the forest panel in inches, (0.62 - 0.27) of the 8.2 inch default.
PANEL_IN = (0.62 - 0.27) * 8.2


def figure_size(labels, k):
    """The figure size, as .fm_size() in R (F81).

    The label column grows with the longest line of the wrapped labels, at
    least the 2.21 inches of the default layout, and each extra line of a
    wrapped label adds 0.2 inch of height.

    Args:
        labels: The wrapped row labels.
        k: The count of studies.

    Returns:
        tuple: The width, the height and the label column width, in inches.
    """
    parts = [str(s).split("\n") for s in labels]
    label_in = max(0.27 * 8.2, 0.25 + 0.075 * max(len(p) for lines in parts for p in lines))
    w = 8.2 - 0.27 * 8.2 + label_in
    h = 1.6 + 0.3 * (k + 3) + 0.2 * sum(len(lines) - 1 for lines in parts)
    return w, h, label_in


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
    # A long study label wraps at 30 characters, and the figure grows with the
    # longest line, so the panel keeps its width (F81), as in R.
    ylabels = [wrap(s, 30) for s in studies["label"]] + ["Fixed effect", "Random effects (DL)"]
    w, h, label_in = figure_size(ylabels, k)
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
        # The same 1-2-5 ticks as the R engine, with at least one tick on each
        # side of the null line, and the axis widened to show every tick.
        lo = float(min(studies["ci_low"].min(), pooled["ci_low"].min()))
        hi = float(max(studies["ci_high"].max(), pooled["ci_high"].max()))
        ticks = log_ticks_ref(lo, hi, 1.0)
        llo, lhi = np.log10(min(lo, ticks[0])), np.log10(max(hi, ticks[-1]))
        pad = 0.05 * (lhi - llo)
        ax.set_xscale("log")
        ax.set_xlim(10 ** (llo - pad), 10 ** (lhi + pad))
        ax.set_xticks(ticks)
        ax.xaxis.set_minor_locator(NullLocator())
        plain_log_axis(ax.xaxis)
    ax.set_yticks(list(sy) + list(py))
    ax.set_yticklabels(ylabels, fontsize=10)
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
    fig.subplots_adjust(left=label_in / w, right=(label_in + PANEL_IN) / w, top=1 - 0.6 / h,
                        bottom=0.9 / h)
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
        f"{re['tau2']:.4f}{_truncation_text(re['q'], re['q_df'])}. Squares are sized by the "
        "random-effects weight and the "
        f"diamonds show the pooled estimates. Computed in Python {platform.python_version()} with "
        f"statsmodels {v('statsmodels')} (combine_effects for the fixed effect) and drawn with "
        f"matplotlib {v('matplotlib')}."
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


def _truncation_text(q, q_df):
    """Return the truncation clause of the methods paragraph.

    The DerSimonian-Laird tau2 is truncated at zero only when Q is below its
    degrees of freedom, so the paragraph names the truncation only then.

    Args:
        q: Cochran's Q.
        q_df: The degrees of freedom of Q.

    Returns:
        str: The clause, or an empty string when no truncation occurred.
    """
    if np.isfinite(q) and np.isfinite(q_df) and q < q_df:
        return ("; Q was below its degrees of freedom, so the DerSimonian-Laird tau2 was "
                "truncated at zero and I2 is 0%")
    return ""
