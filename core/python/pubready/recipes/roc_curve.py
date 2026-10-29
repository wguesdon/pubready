"""Recipe: ROC curves of one or more numeric markers against a binary outcome.

Mirrors the R recipe, which uses pROC. The rules of pROC are restated here so
that the two engines agree to the last digit:

* The case level is --positive, or the second level in code-point order.
* The direction is "<" (higher values indicate the case) when the control median
  is not above the case median, and ">" otherwise, so every AUC is at least 0.5.
* The thresholds are -Inf, the midpoints between sorted unique values, and +Inf,
  in the order that takes the sensitivity from 1 down to 0.
* The AUC confidence interval and the paired test follow DeLong et al. (1988),
  and the interval is clipped to [0, 1].
"""
import functools

import numpy as np
import pandas as pd
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette


def level_strings(s):
    """Turn an outcome column into strings that match the R engine.

    A whole-number float column (1.0, 0.0) becomes "1" and "0", as R prints it.

    Args:
        s: The outcome column.

    Returns:
        pandas.Series: The levels as strings, with missing values kept as NaN.
    """
    if pd.api.types.is_float_dtype(s):
        nonmiss = s.dropna()
        if len(nonmiss) and np.all(np.mod(nonmiss, 1) == 0):
            return s.map(whole_number_str)
    return s.astype(str).where(s.notna())


def whole_number_str(v):
    """Format a whole-number float as an integer string, or NaN for a missing value.

    Args:
        v: The value.

    Returns:
        str or float: The string, or NaN.
    """
    if pd.isna(v):
        return np.nan
    return str(int(v))


def direction_of(controls, cases):
    """Pick the direction with the pROC rule.

    Args:
        controls: The marker values of the controls.
        cases: The marker values of the cases.

    Returns:
        str: "<" or ">".
    """
    return "<" if np.median(controls) <= np.median(cases) else ">"


def placements(controls, cases, direction):
    """DeLong placement values of the cases and of the controls.

    Args:
        controls: The marker values of the controls.
        cases: The marker values of the cases.
        direction: "<" or ">".

    Returns:
        tuple: The case placements (V10) and the control placements (V01).
    """
    c_s, k_s = np.sort(controls), np.sort(cases)
    n, m = len(c_s), len(k_s)
    c_lt = np.searchsorted(c_s, cases, side="left")
    c_le = np.searchsorted(c_s, cases, side="right")
    k_lt = np.searchsorted(k_s, controls, side="left")
    k_le = np.searchsorted(k_s, controls, side="right")
    if direction == "<":
        v10 = (c_lt + 0.5 * (c_le - c_lt)) / n
        v01 = (m - k_le + 0.5 * (k_le - k_lt)) / m
    else:
        v10 = (n - c_le + 0.5 * (c_le - c_lt)) / n
        v01 = (k_lt + 0.5 * (k_le - k_lt)) / m
    return v10, v01


def roc_points(controls, cases, direction):
    """Thresholds, sensitivities and specificities in pROC order.

    Args:
        controls: The marker values of the controls.
        cases: The marker values of the cases.
        direction: "<" or ">".

    Returns:
        tuple: Three arrays, ordered so that the sensitivity falls from 1 to 0.
    """
    vals = np.unique(np.concatenate([controls, cases]))
    thr = np.concatenate([[-np.inf], (vals[:-1] + vals[1:]) / 2, [np.inf]])
    c_s, k_s = np.sort(controls), np.sort(cases)
    n, m = len(c_s), len(k_s)
    if direction == "<":
        sens = (m - np.searchsorted(k_s, thr, side="right")) / m
        spec = np.searchsorted(c_s, thr, side="left") / n
    else:
        thr = thr[::-1]
        sens = np.searchsorted(k_s, thr, side="left") / m
        spec = (n - np.searchsorted(c_s, thr, side="right")) / n
    return thr, sens, spec


def analyse_marker(controls, cases):
    """Compute the curve, the AUC with its DeLong CI and the Youden point.

    Args:
        controls: The marker values of the controls.
        cases: The marker values of the cases.

    Returns:
        dict: The direction, the placements, the AUC and its variance, the CI,
        the curve, and the Youden threshold with its sensitivity and specificity.
    """
    direction = direction_of(controls, cases)
    v10, v01 = placements(controls, cases, direction)
    m, n = len(cases), len(controls)
    auc = float(v10.mean())
    var = float(np.var(v10, ddof=1) / m + np.var(v01, ddof=1) / n)
    z = stats.norm.ppf(0.975)
    lo, hi = min(max(auc - z * np.sqrt(var), 0.0), 1.0), min(max(auc + z * np.sqrt(var), 0.0), 1.0)
    thr, sens, spec = roc_points(controls, cases, direction)
    i = int(np.argmax(sens + spec - 1))
    return {"direction": direction, "v10": v10, "v01": v01, "auc": auc, "var": var,
            "ci_low": lo, "ci_high": hi, "thr": thr, "sens": sens, "spec": spec,
            "youden_threshold": float(thr[i]), "sensitivity": float(sens[i]),
            "specificity": float(spec[i])}


def delong_paired(a, b):
    """DeLong's test for two correlated ROC curves on the same subjects.

    Args:
        a: The analysis of the first marker.
        b: The analysis of the second marker.

    Returns:
        tuple: The Z statistic of (AUC a - AUC b) and the two-sided p-value.
    """
    m, n = len(a["v10"]), len(a["v01"])
    cov = (np.cov(a["v10"], b["v10"], ddof=1)[0, 1] / m +
           np.cov(a["v01"], b["v01"], ddof=1)[0, 1] / n)
    zval = (a["auc"] - b["auc"]) / np.sqrt(a["var"] + b["var"] - 2 * cov)
    return float(zval), float(2 * stats.norm.sf(abs(zval)))


def legend_handles(ax, pal):
    """Legend keys that show each curve with its hollow Youden point, as R draws them.

    Args:
        ax: The ROC axes, with one labelled curve per marker.
        pal: The colour of each marker, in curve order.

    Returns:
        list: One Line2D key per marker, labelled as its curve.
    """
    from matplotlib.lines import Line2D
    curves = [ln for ln in ax.get_lines() if not ln.get_label().startswith("_")]
    return [Line2D([0], [0], color=pal[k], lw=1.6, marker="o", markersize=6,
                   markerfacecolor="white", markeredgecolor=pal[k], markeredgewidth=1.4,
                   label=ln.get_label()) for k, ln in enumerate(curves)]


def keep_right_margin(fig, ax, inches=0.2):
    """Keep white space right of the last x tick label, as the R plot margin does.

    The bundle writer crops the figure to its artists, so the right-most tick
    label "1.0" would end at the edge of the image. An invisible marker at the
    given distance right of the axes keeps a margin there.

    Args:
        fig: The figure, after its final layout.
        ax: The ROC axes.
        inches: The margin right of the axes in inches.
    """
    from matplotlib.patches import Rectangle
    x = ax.get_position().x1 + inches / fig.get_figwidth()
    y = ax.get_position().y0
    fig.add_artist(Rectangle((min(x, 1.0), y), 1e-3, 1e-3, transform=fig.transFigure,
                             fill=False, edgecolor="none"))


def recipe_roc_curve(df, spec):
    """Draw ROC curves with AUCs, DeLong intervals and the Youden points.

    Args:
        df: The tidy table, one row per subject.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if the outcome does not have exactly
            two levels, or if a group has fewer than two subjects.
    """
    y, x_raw = spec["data"].get("y"), spec["data"].get("x")
    if not y or not x_raw:
        raise ValueError("roc_curve needs --y (the binary outcome) and --x "
                         "(one marker, or several separated by commas)")
    markers = [m.strip() for m in str(x_raw).split(",") if m.strip()]
    for c in [y] + markers:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({"outcome": level_strings(df[y])})
    for m in markers:
        d[m] = pd.to_numeric(df[m], errors="coerce")
    keep = d["outcome"].notna() & d[markers].notna().all(axis=1)
    d = d[keep].reset_index(drop=True)
    lv = sorted(d["outcome"].unique())
    if len(lv) != 2:
        raise ValueError(f"roc_curve needs exactly 2 outcome levels in '{y}'; found {len(lv)} "
                         f"({', '.join(lv)})")
    case = str(spec["data"]["positive"]) if spec["data"].get("positive") is not None else lv[1]
    if case not in lv:
        raise ValueError(f"--positive '{case}' is not a level of '{y}' ({', '.join(lv)})")
    control = [v for v in lv if v != case][0]
    is_case = (d["outcome"] == case).to_numpy()
    n_case, n_ctrl = int(is_case.sum()), int((~is_case).sum())
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing outcome or marker value, "
                     "so every marker is read on the same subjects.")
    steps.append(f"Kept {len(d)} of {n0} row(s): {n_case} case(s) ({case}) and {n_ctrl} "
                 f"control(s) ({control}).")
    if n_case < 2 or n_ctrl < 2:
        raise ValueError("roc_curve needs at least 2 cases and 2 controls")

    res = [analyse_marker(d.loc[~is_case, m].to_numpy(float), d.loc[is_case, m].to_numpy(float))
           for m in markers]
    rows = []
    for k, (m, r) in enumerate(zip(markers, res)):
        zval, pval = delong_paired(res[0], r) if k > 0 else (np.nan, np.nan)
        rows.append({"marker": m, "case_level": case, "n_cases": n_case, "n_controls": n_ctrl,
                     "direction": r["direction"], "auc": r["auc"],
                     "auc_ci_low": r["ci_low"], "auc_ci_high": r["ci_high"],
                     "youden_threshold": r["youden_threshold"], "sensitivity": r["sensitivity"],
                     "specificity": r["specificity"], "delong_z_vs_first": zval,
                     "delong_p_vs_first": pval})
    stats_df = pd.DataFrame(rows)

    ap = spec["appearance"]
    pal = palette(len(markers), ap.get("palette"))
    w, h = 4.8, 5.0
    fig, ax = new_fig(w, h)
    ax.plot([0, 1], [0, 1], ls="--", color="0.6", lw=0.8)
    for k, (m, r) in enumerate(zip(markers, res)):
        lab = f"{m}: AUC {r['auc']:.2f} ({r['ci_low']:.2f} to {r['ci_high']:.2f})"
        ax.plot(1 - r["spec"], r["sens"], color=pal[k], lw=1.6, label=lab)
        ax.scatter([1 - r["specificity"]], [r["sensitivity"]], s=36, facecolor="white",
                   edgecolor=pal[k], linewidths=1.4, zorder=3)
    sub = f"{n_case} {case} vs {n_ctrl} {control}"
    if len(markers) == 2:
        sub += f"; DeLong P = {rows[1]['delong_p_vs_first']:.2g}"
    ax.set_title(sub, fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01); ax.set_aspect("equal")
    ax.set_xlabel(ap.get("x_label") or "1 - Specificity")
    ax.set_ylabel(ap.get("y_label") or "Sensitivity")
    # The legend keys carry the Youden point as well as the line, as in R.
    ax.legend(handles=legend_handles(ax, pal), loc="lower right", frameon=False, fontsize=9)
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()
    keep_right_margin(fig, ax)

    test_meta = {
        "name": "ROC curve analysis (DeLong)",
        "case_level": case, "control_level": control,
        "n_cases": n_case, "n_controls": n_ctrl,
        "markers": {r["marker"]: {k: r[k] for k in (
            "auc", "auc_ci_low", "auc_ci_high", "direction", "youden_threshold",
            "sensitivity", "specificity", "delong_p_vs_first")} for r in rows},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, case, control, n_case, n_ctrl, y),
        "df_used": d, "clean_steps": steps,
        "label": slugify("_".join([y] + markers)), "width": w, "height": h,
        "build_script": functools.partial(_script, spec, markers, case, control),
    }


def _methods(rows, case, control, n_case, n_ctrl, y):
    """Write the methods paragraph.

    Args:
        rows: The stats table.
        case: The case level.
        control: The control level.
        n_case: The count of cases.
        n_ctrl: The count of controls.
        y: The outcome column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    per = "; ".join(
        f"{r['marker']}, AUC = {r['auc']:.3f} (95% CI {r['auc_ci_low']:.3f} to "
        f"{r['auc_ci_high']:.3f}), optimal threshold {r['youden_threshold']:.4g} "
        f"(sensitivity {r['sensitivity']:.2f}, specificity {r['specificity']:.2f})"
        for _, r in rows.iterrows())
    cmp = ""
    if len(rows) >= 2:
        cmp = (f" The AUC of each marker was compared with that of {rows['marker'].iloc[0]} by "
               "DeLong's test for two correlated ROC curves: " + "; ".join(
                   f"{r['marker']}, Z = {r['delong_z_vs_first']:.2f}, P = {r['delong_p_vs_first']:.2g}"
                   for _, r in rows.iloc[1:].iterrows()) + ".")
    flipped = list(rows.loc[rows["direction"] == ">", "marker"])
    flip = (f" For {', '.join(flipped)}, lower values indicate the case, so the curve reads the "
            "marker in that direction.") if flipped else ""
    return (
        f"Receiver operating characteristic (ROC) curves were computed for the discrimination of "
        f"{y} = {case} (cases, n = {n_case}) from {control} (controls, n = {n_ctrl}). The area "
        "under each curve (AUC) is reported with a 95% confidence interval by the method of "
        "DeLong et al. (1988), and the optimal threshold maximises the Youden index "
        f"(sensitivity + specificity - 1). {per}.{flip}{cmp} Computed in Python "
        f"{platform.python_version()} with numpy {v('numpy')} and scipy {v('scipy')}, following "
        f"the rules of pROC, and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, markers, case, control, in_name, fig_stub):
    """Write the standalone script that redraws the curves.

    Args:
        spec: The figure specification.
        markers: The marker columns.
        case: The case level.
        control: The control level.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    y = spec["data"]["y"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The AUC is
# the Mann-Whitney probability that a case outranks a control.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
markers = {markers!r}
raw = raw.dropna(subset=["{y}"] + markers)
is_case = raw["{y}"].astype(str).str.replace(r"\\.0$", "", regex=True) == "{case}"
fig, ax = plt.subplots(figsize=(4.8, 5.0))
ax.plot([0, 1], [0, 1], ls="--", color="0.6")
for m in markers:
    cases, controls = raw.loc[is_case, m].to_numpy(), raw.loc[~is_case, m].to_numpy()
    sign = 1 if np.median(controls) <= np.median(cases) else -1
    gt = (sign * cases[:, None] > sign * controls[None, :]).mean()
    eq = (cases[:, None] == controls[None, :]).mean()
    auc = gt + 0.5 * eq
    thr = np.sort(np.unique(sign * raw[m].to_numpy()))
    sens = [(sign * cases >= t).mean() for t in thr]
    fpr = [(sign * controls >= t).mean() for t in thr]
    ax.plot(fpr, sens, label=f"{{m}}: AUC {{auc:.3f}}")
    print(m, "AUC", auc)
ax.set_xlabel("1 - Specificity"); ax.set_ylabel("Sensitivity"); ax.set_aspect("equal")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure ({control} = control, {case} = case).")
'''
