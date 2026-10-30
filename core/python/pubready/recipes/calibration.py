"""Recipe: calibration of a risk prediction model.

Mirrors the R recipe. The outcome is --y and the predicted probabilities are
--x, one column per model. For each model: the Brier score, the calibration
intercept (calibration-in-the-large) and slope from logistic models with a 95%
Wald CI, the Hosmer-Lemeshow test and the observed/expected ratio.

The Hosmer-Lemeshow groups restate ``ResourceSelection::hoslem.test``: the
breaks are the unique type-7 quantiles of p at 0, 1/g, ..., 1, and R's
``cut(include.lowest = TRUE)`` assigns the subjects. ``quantile7`` follows the
arithmetic of R's ``quantile.default`` step by step, so a subject that sits on a
break lands in the same group in both engines.
"""
import functools
import inspect
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

from ..bundle import slugify
from ..theme import apply_pub_style, palette
from .roc_curve import level_strings

#: The number of Hosmer-Lemeshow groups requested.
G = 10
#: Predictions are clipped to [CLIP, 1 - CLIP] before the logit.
CLIP = 1e-8
#: Below this many distinct predictions the slope and the HL test are not estimable.
MIN_DISTINCT = 3
#: A legend entry longer than this puts the model name on its own line.
WRAP = 62


def quantile7(x, probs):
    """Type-7 sample quantiles with the arithmetic of R's quantile.default.

    Args:
        x: The values.
        probs: The probabilities.

    Returns:
        numpy.ndarray: The quantiles.
    """
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    index = 1 + max(n - 1, 0) * np.asarray(probs, dtype=float)
    lo = np.floor(index).astype(int)
    hi = np.ceil(index).astype(int)
    qs = x[lo - 1].copy()
    sel = (index > lo) & (x[hi - 1] != qs)
    h = (index - lo)[sel]
    qs[sel] = (1 - h) * qs[sel] + h * x[hi[sel] - 1]
    return qs


def hl_groups(p, g):
    """The Hosmer-Lemeshow group of each subject, as hoslem.test forms it.

    Args:
        p: The predicted probabilities.
        g: The number of groups requested.

    Returns:
        numpy.ndarray: The group of each subject, from 1.
    """
    probs = np.minimum(0.0 + np.arange(g + 1) * (1 / g), 1.0)
    qq = np.unique(quantile7(p, probs))
    if len(qq) < 2:
        return np.ones(len(p), dtype=int)
    code = np.searchsorted(qq, p, side="left")
    code[code == 0] = 1
    return code


def group_table(p, y, g):
    """Per-group counts for the test and the figure.

    Args:
        p: The predicted probabilities.
        y: The outcome, 1 for an event.
        g: The number of groups requested.

    Returns:
        pandas.DataFrame: For each non-empty group: the size, the events, the
        expected events and non-events, the mean prediction, the observed
        proportion and its 95% Wilson interval.
    """
    grp = hl_groups(p, g)
    ids = np.unique(grp)
    n = np.array([np.sum(grp == i) for i in ids], dtype=float)
    o1 = np.array([np.sum(y[grp == i]) for i in ids], dtype=float)
    e1 = np.array([np.sum(p[grp == i]) for i in ids])
    e0 = np.array([np.sum(1 - p[grp == i]) for i in ids])
    obs = o1 / n
    z = stats.norm.ppf(0.975)
    centre = (obs + z ** 2 / (2 * n)) / (1 + z ** 2 / n)
    half = z * np.sqrt(obs * (1 - obs) / n + z ** 2 / (4 * n ** 2)) / (1 + z ** 2 / n)
    return pd.DataFrame({"group": ids, "n": n, "events": o1, "expected": e1,
                         "expected_non": e0, "mean_pred": e1 / n, "observed": obs,
                         "ci_low": centre - half, "ci_high": centre + half})


def hosmer_lemeshow(tab):
    """The Hosmer-Lemeshow statistic, its df and its p-value.

    Args:
        tab: The per-group table.

    Returns:
        dict: The statistic, df = groups - 2, the p-value and the group count.
    """
    stat = float(np.sum((tab["events"] - tab["expected"]) ** 2 / tab["expected"]) +
                 np.sum(((tab["n"] - tab["events"]) - tab["expected_non"]) ** 2 /
                        tab["expected_non"]))
    df = len(tab) - 2
    p = float(stats.chi2.sf(stat, df)) if df > 0 and np.isfinite(stat) else np.nan
    return {"statistic": stat, "df": df, "p_value": p, "groups": len(tab)}


def wald(fit, k):
    """Coefficient and 95% Wald CI of one term of a fitted GLM.

    Args:
        fit: The statsmodels GLM result.
        k: The index of the term.

    Returns:
        dict: The estimate, the standard error and the CI bounds.
    """
    est, se = float(fit.params[k]), float(fit.bse[k])
    z = stats.norm.ppf(0.975)
    return {"estimate": est, "se": se, "ci_low": est - z * se, "ci_high": est + z * se}


def fit_model(m, d, case):
    """The calibration statistics of one model.

    Args:
        m: The prediction column.
        d: The cleaned table.
        case: The event level.

    Returns:
        dict: The statistics, the per-group table and the lowess smooth.
    """
    p = d[m].to_numpy(float)
    y = (d["outcome"] == case).to_numpy().astype(int)
    lp = np.log(np.clip(p, CLIP, 1 - CLIP) / (1 - np.clip(p, CLIP, 1 - CLIP)))
    ones = np.ones((len(y), 1))
    fam = sm.families.Binomial()
    x_slope = np.column_stack([ones, lp])
    # Below 3 distinct predictions the data cannot identify a slope, and the
    # Hosmer-Lemeshow test has no degree of freedom: both are reported as not
    # estimable (F78), as in R.
    estimable = len(np.unique(p)) >= MIN_DISTINCT
    fit_int = sm.GLM(y, ones, family=fam, offset=lp).fit(tol=1e-12, maxiter=100)
    # A second fit that starts at the estimate reads the standard error at the
    # estimate itself, as the R engine does.
    fit_int = sm.GLM(y, ones, family=fam, offset=lp).fit(start_params=fit_int.params, tol=1e-12,
                                                         maxiter=100)
    nan = float("nan")
    slope = {"estimate": nan, "se": nan, "ci_low": nan, "ci_high": nan}
    if estimable:
        fit_slope = sm.GLM(y, x_slope, family=fam).fit(tol=1e-12, maxiter=100)
        fit_slope = sm.GLM(y, x_slope, family=fam).fit(start_params=fit_slope.params, tol=1e-12,
                                                        maxiter=100)
        slope = wald(fit_slope, 1)
    tab = group_table(p, y, G)
    hl = (hosmer_lemeshow(tab) if estimable else
          {"statistic": nan, "df": None, "p_value": nan, "groups": len(tab)})
    smooth = lowess(y, p, frac=2 / 3, it=0, delta=0.01 * (p.max() - p.min()),
                    return_sorted=True)
    return {"model": m, "n": len(y), "n_events": int(y.sum()),
            "brier": float(np.mean((p - y) ** 2)), "estimable": estimable,
            "intercept": wald(fit_int, 0), "slope": slope,
            "hl": hl, "oe_ratio": float(y.sum() / p.sum()),
            "groups": tab, "smooth_x": smooth[:, 0], "smooth_y": np.clip(smooth[:, 1], 0, 1),
            "p": p, "y": y}


def p_text(p):
    """Format a p-value with its relation sign.

    Args:
        p: The p-value.

    Returns:
        str: "< 0.001", "= 0.xx" or "NA".
    """
    if not np.isfinite(p):
        return "NA"
    return "< 0.001" if p < 0.001 else f"= {p:.2g}"


def legend_label(r):
    """The legend entry of one model.

    An entry longer than WRAP characters puts the model name on its own line, so
    the legend stays within the figure width.

    Args:
        r: The model statistics.

    Returns:
        str: The label.
    """
    slope = f"{r['slope']['estimate']:.2f}" if r["estimable"] else "n.e."
    st = (f"slope {slope}, intercept {r['intercept']['estimate']:.2f}, "
          f"Brier {r['brier']:.3f}, HL P {p_text(r['hl']['p_value'])}")
    one = f"{r['model']}: {st}"
    return f"{r['model']}\n{st}" if len(one) > WRAP else one


def band_base(k):
    """The baseline of the histogram band of model k (from 1).

    Args:
        k: The model index, from 1.

    Returns:
        float: The y position.
    """
    return -0.17 - 0.24 * (k - 1)


def figure_dims(n_models, labels):
    """The figure size, matched to the R engine.

    Args:
        n_models: The count of models.
        labels: The legend entries.

    Returns:
        dict: The width and the height in inches.
    """
    span = 1.03 - (band_base(n_models) - 0.11)
    lines = sum(1 + ("\n" in lab) for lab in labels)
    return {"width": 5.2, "height": 4.3 * span / 1.04 + 1.35 + 0.2 * lines}


def spike_counts(r):
    """Counts of events and non-events in 100 bins of width 0.01.

    Args:
        r: The model statistics.

    Returns:
        tuple: The event counts and the non-event counts.
    """
    b = np.minimum(np.floor(r["p"] * 100), 99).astype(int)
    return (np.bincount(b[r["y"] == 1], minlength=100),
            np.bincount(b[r["y"] == 0], minlength=100))


def recipe_calibration(df, spec):
    """Draw the calibration plot of one or more risk models.

    Args:
        df: The tidy table, one row per subject.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a prediction is outside [0, 1],
            or if the outcome does not have exactly two levels.
    """
    y, x_raw = spec["data"].get("y"), spec["data"].get("x")
    if not y or not x_raw:
        raise ValueError("calibration needs --y (the binary outcome) and --x (one column of "
                         "predicted probabilities, or several separated by commas)")
    models = [m.strip() for m in str(x_raw).split(",") if m.strip()]
    for c in [y] + models:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({"outcome": level_strings(df[y])})
    for m in models:
        d[m] = pd.to_numeric(df[m], errors="coerce")
    keep = d["outcome"].notna() & d[models].notna().all(axis=1)
    d = d[keep].reset_index(drop=True)
    for m in models:
        if ((d[m] < 0) | (d[m] > 1)).any():
            raise ValueError(f"column '{m}' holds values outside [0, 1]; calibration needs "
                             "predicted probabilities")
    lv = sorted(d["outcome"].unique())
    if len(lv) != 2:
        raise ValueError(f"calibration needs exactly 2 outcome levels in '{y}'; found {len(lv)} "
                         f"({', '.join(lv)})")
    pos = spec["data"].get("positive")
    case = str(pos) if pos is not None else lv[1]
    if case not in lv:
        raise ValueError(f"--positive '{case}' is not a level of '{y}' ({', '.join(lv)})")
    control = [v for v in lv if v != case][0]
    n_case = int((d["outcome"] == case).sum())
    if n_case < 2 or len(d) - n_case < 2:
        raise ValueError("calibration needs at least 2 events and 2 non-events")
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing outcome or prediction, "
                     "so every model is read on the same subjects.")
    steps.append(f"Kept {len(d)} of {n0} row(s): {n_case} event(s) ({case}) and "
                 f"{len(d) - n_case} non-event(s) ({control}).")
    steps.append(f"Predictions were clipped to [{CLIP:g}, 1 - {CLIP:g}] before the logit "
                 "transform.")

    fits = [fit_model(m, d, case) for m in models]
    stats_df = pd.DataFrame([{
        "model": r["model"], "case_level": case, "n": r["n"], "n_events": r["n_events"],
        "brier": r["brier"], "intercept": r["intercept"]["estimate"],
        "intercept_ci_low": r["intercept"]["ci_low"],
        "intercept_ci_high": r["intercept"]["ci_high"],
        "slope": r["slope"]["estimate"], "slope_ci_low": r["slope"]["ci_low"],
        "slope_ci_high": r["slope"]["ci_high"], "hl_statistic": r["hl"]["statistic"],
        "hl_df": r["hl"]["df"], "hl_p": r["hl"]["p_value"], "hl_groups": r["hl"]["groups"],
        "oe_ratio": r["oe_ratio"]} for r in fits])

    ap = spec["appearance"]
    dims = figure_dims(len(models), [legend_label(r) for r in fits])
    fig = draw_figure(fits, palette(len(models), ap.get("palette")), ap,
                      f"{len(d)} subjects, {n_case} events ({y} = {case})", dims)
    test_meta = {
        "name": "Calibration (Hosmer-Lemeshow test, calibration intercept and slope)",
        "case_level": case, "control_level": control, "n": len(d), "n_events": n_case,
        "groups_requested": G,
        "models": {r["model"]: {
            "n": r["n"], "n_events": r["n_events"], "brier": r["brier"],
            "slope_estimable": r["estimable"],
            "intercept": r["intercept"]["estimate"], "slope": r["slope"]["estimate"],
            "hosmer_lemeshow": r["hl"], "oe_ratio": r["oe_ratio"]} for r in fits},
    }
    groups = pd.concat([r["groups"].assign(model=r["model"]) for r in fits], ignore_index=True)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(fits, case, control, n_case, len(d), y),
        "df_used": d, "clean_steps": steps, "groups": groups,
        "label": slugify("_".join([y] + models)), "width": dims["width"],
        "height": dims["height"],
        "build_script": functools.partial(_script, spec, models, case, dims),
    }


def draw_figure(fits, pal, ap, sub, dims):
    """Draw the calibration panel with the spike histogram bands below it.

    Args:
        fits: The model statistics.
        pal: The colour of each model.
        ap: The appearance options.
        sub: The subtitle.
        dims: The figure size.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    n_models = len(fits)
    ylow = band_base(n_models) - 0.11
    fig = plt.figure(figsize=(dims["width"], dims["height"]), layout="constrained")
    ax = fig.add_subplot()
    ax.plot([0, 1], [0, 1], ls="--", color="0.6", lw=0.8, zorder=1)
    counts = [spike_counts(r) for r in fits]
    top = max(max(up.max(), down.max()) for up, down in counts)
    xs = (np.arange(100) + 0.5) / 100
    for k, (r, (up, down)) in enumerate(zip(fits, counts), start=1):
        base = band_base(k)
        ax.hlines(base, 0, 1, color="0.6", lw=0.6)
        sel_u, sel_d = up > 0, down > 0
        ax.vlines(xs[sel_u], base, base + 0.09 * up[sel_u] / top, color=pal[k - 1], lw=0.9)
        ax.vlines(xs[sel_d], base, base - 0.09 * down[sel_d] / top, color=pal[k - 1], lw=0.9)
        for dy, txt in ((0.045, "events"), (-0.045, "non-events")):
            ax.text(-0.025, base + dy, txt, ha="right", va="center", fontsize=7.4,
                    color="0.35", clip_on=False)
    handles = []
    for k, r in enumerate(fits):
        off = (k + 1 - (n_models + 1) / 2) * 0.008
        ax.plot(r["smooth_x"], r["smooth_y"], color=pal[k], lw=1.6, zorder=3)
        g = r["groups"]
        ax.errorbar(g["mean_pred"] + off, g["observed"],
                    yerr=[g["observed"] - g["ci_low"], g["ci_high"] - g["observed"]],
                    fmt="o", color=pal[k], ms=4.5, capsize=2.5, elinewidth=0.9, zorder=4)
        handles.append(Line2D([0], [0], color=pal[k], lw=1.6))
    ticks = [0, 0.25, 0.5, 0.75, 1.0]
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(ylow, 1.03)
    ax.set_aspect("equal")
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.2f}" for t in ticks])
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:.2f}" for t in ticks])
    ax.set_xlabel(ap.get("x_label") or "Predicted probability")
    ax.set_ylabel(ap.get("y_label") or "Observed proportion", y=(0.5 - ylow) / (1.03 - ylow))
    ax.set_title(sub, fontsize=10, color="#4d4d4d")
    apply_pub_style(ax, ap.get("theme"))
    ax.spines["left"].set_bounds(0, 1)
    ax.spines["bottom"].set_bounds(0, 1)
    fig.legend(handles, [legend_label(r) for r in fits], loc="outside lower center", ncol=1,
               frameon=False, fontsize=9)
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    return fig


def model_text(r):
    """The methods sentence of one model, as .cal_model_text() in R.

    Args:
        r: The model statistics.

    Returns:
        str: The sentence part.
    """
    head = (f"{r['model']}: Brier score {r['brier']:.3f}, calibration intercept "
            f"{r['intercept']['estimate']:.2f} (95% CI {r['intercept']['ci_low']:.2f} to "
            f"{r['intercept']['ci_high']:.2f}), ")
    if not r["estimable"]:
        return (f"{head}slope and Hosmer-Lemeshow test not estimable, because the predictions "
                f"take fewer than {MIN_DISTINCT} distinct values, observed/expected "
                f"{r['oe_ratio']:.2f}")
    return (f"{head}slope {r['slope']['estimate']:.2f} (95% CI "
            f"{r['slope']['ci_low']:.2f} to {r['slope']['ci_high']:.2f}), Hosmer-Lemeshow "
            f"chi-square({int(r['hl']['df'])}) = {r['hl']['statistic']:.2f}, P "
            f"{p_text(r['hl']['p_value'])}, observed/expected {r['oe_ratio']:.2f}")


def _methods(fits, case, control, n_case, n, y):
    """Write the methods paragraph.

    Args:
        fits: The model statistics.
        case: The event level.
        control: The non-event level.
        n_case: The count of events.
        n: The count of subjects.
        y: The outcome column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    per = "; ".join(model_text(r) for r in fits)
    return (
        f"The calibration of the predicted probabilities was assessed against the observed "
        f"outcome {y} = {case} (events, n = {n_case}) versus {control} (non-events, n = "
        f"{n - n_case}). For each model, the Brier score is the mean squared difference between "
        "the prediction and the outcome. The calibration intercept (calibration-in-the-large) is "
        "the intercept of a logistic regression of the outcome with logit(p) as an offset, and "
        "the calibration slope is the coefficient of logit(p) in a logistic regression, each "
        "with a 95% Wald CI; ideal values are 0 and 1. The Hosmer-Lemeshow test grouped the "
        f"subjects at the deciles of the prediction (g = {G}; tied deciles merge, df = groups - "
        "2). The figure shows the observed proportion against the mean prediction in each group "
        "with a 95% Wilson CI, a lowess smooth of the outcome on the prediction (f = 2/3, without "
        "robustness iterations, which would treat the events of a binary outcome as outliers)"
        f", and the spike histogram of the predictions. {per}. Computed in "
        f"Python {platform.python_version()} with statsmodels {v('statsmodels')}, grouping as "
        f"ResourceSelection::hoslem.test does, and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, models, case, dims, in_name, fig_stub):
    """Write the standalone script that recomputes the statistics.

    Args:
        spec: The figure specification.
        models: The prediction columns.
        case: The event level.
        dims: The figure size.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    y = spec["data"]["y"]
    helpers = "\n\n".join(inspect.getsource(fn) for fn in
                          (quantile7, hl_groups, group_table, hosmer_lemeshow))
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The helper
# functions below group the subjects as ResourceSelection::hoslem.test does.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

{helpers}

raw = pd.read_csv("{in_name}")
models = {models!r}
raw = raw.dropna(subset=["{y}"] + models)
yb = (raw["{y}"].astype(str).str.replace(r"\\.0$", "", regex=True) == "{case}").to_numpy().astype(int)
fig, ax = plt.subplots(figsize=({dims["width"]:.2f}, {dims["height"]:.2f}))
ax.plot([0, 1], [0, 1], ls="--", color="0.6")
for m in models:
    p = raw[m].to_numpy(float)
    lp = np.log(np.clip(p, 1e-8, 1 - 1e-8) / (1 - np.clip(p, 1e-8, 1 - 1e-8)))
    ones = np.ones((len(p), 1))
    print(m, "Brier", np.mean((p - yb) ** 2))
    print(sm.GLM(yb, ones, family=sm.families.Binomial(), offset=lp).fit().summary2().tables[1])
    print(sm.GLM(yb, np.column_stack([ones, lp]), family=sm.families.Binomial()).fit().summary2().tables[1])
    print("Hosmer-Lemeshow", hosmer_lemeshow(group_table(p, yb, 10)))
    s = lowess(yb, p, frac=2 / 3, it=0, delta=0.01 * (p.max() - p.min()))
    ax.plot(s[:, 0], s[:, 1], label=m)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed proportion")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
