"""Recipe: a repeated-measures time course. Mirrors the R recipe.

Each subject (--id) is measured at every time point (--x) and may belong to a
group (--group). The figure shows the group mean +/- SEM at each time with thin
lines for the subjects. The test is a mixed ANOVA (group between subjects, time
within subjects), or a one-way repeated-measures ANOVA when there is no group.

The sums of squares are the classical split-plot decomposition, which equals R's
aov with an Error(id/time) stratum. The Greenhouse-Geisser epsilon and Mauchly's
test come from the pooled within-group covariance of the time points, as
car::Anova computes them. With --test auto the within-subject effects take the
Greenhouse-Geisser correction when Mauchly's test rejects sphericity; --test gg
always applies it and --test uncorrected never does. A subject without a value
at every time point is dropped, and the drop is recorded.
"""
import functools

import numpy as np
import pandas as pd
from scipy import linalg, stats

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import cap_first


def p_label(p):
    """Format a P value for the figure and the methods.

    Args:
        p: The P value.

    Returns:
        str: "P < 0.001" or "P = <value>".
    """
    return "P < 0.001" if p < 0.001 else f"P = {p:.2g}"


def mauchly(m, p, df):
    """Mauchly's sphericity test, as stats::mauchly.test and car::Anova compute it.

    The correction term of R reads 3 * (p + 1), the number of time points, at one
    position, and it is kept so that the P value matches R exactly.

    Args:
        m: The contrast-space sums of squares and products matrix.
        p: The number of contrasts (time points - 1).
        df: The error degrees of freedom (subjects - groups).

    Returns:
        dict: W, the chi-square statistic, its df and the P value.
    """
    if p < 2:
        return {"W": 1.0, "statistic": 0.0, "df": 0, "p_value": 1.0}
    log_w = np.log(np.linalg.det(m)) - p * np.log(np.trace(m) / p)
    rho = 1 - (2 * p ** 2 + p + 2) / (6 * p * df)
    w2 = ((p + 2) * (p - 1) * (p - 2) * (2 * p ** 3 + 6 * p ** 2 + 3 * (p + 1) + 2) /
          (288 * (df * p * rho) ** 2))
    z = -df * rho * log_w
    f = p * (p + 1) / 2 - 1
    pr1, pr2 = stats.chi2.sf(z, f), stats.chi2.sf(z, f + 4)
    return {"W": float(np.exp(log_w)), "statistic": float(z), "df": int(f),
            "p_value": float(pr1 + w2 * (pr2 - pr1))}


def mixed_anova(y, g, levels):
    """Split-plot ANOVA of a subjects-by-times matrix.

    Args:
        y: The subjects-by-times matrix.
        g: The group of each row.
        levels: The group levels, in order.

    Returns:
        dict: The effect rows, the GG epsilon, Mauchly's test and the cell means.
    """
    n, t = y.shape
    grand = y.mean()
    n_g = np.array([(g == lv).sum() for lv in levels], dtype=float)
    cells = np.vstack([y[g == lv].mean(axis=0) for lv in levels])
    ss_subj = t * np.sum((y.mean(axis=1) - grand) ** 2)
    ss_group = t * np.sum(n_g * (cells.mean(axis=1) - grand) ** 2)
    ss_time = n * np.sum((y.mean(axis=0) - grand) ** 2)
    ss_cells = np.sum(n_g * np.sum((cells - grand) ** 2, axis=1))
    ss_inter = ss_cells - ss_group - ss_time
    ss_err_b = ss_subj - ss_group
    ss_err_w = np.sum((y - grand) ** 2) - ss_subj - ss_time - ss_inter
    k = len(levels)
    df_err_b, df_err_w = n - k, (n - k) * (t - 1)
    ssp = sum((y[g == lv] - y[g == lv].mean(axis=0)).T @ (y[g == lv] - y[g == lv].mean(axis=0))
              for lv in levels)
    c = linalg.null_space(np.ones((1, t)))
    m = c.T @ ssp @ c
    p = t - 1
    eps = 1.0 if p == 1 else float(np.trace(m) ** 2 / (p * np.trace(m @ m)))
    rows = [{"term": "time", "ss": ss_time, "df1": t - 1, "ss_err": ss_err_w, "df2": df_err_w,
             "within": True}]
    if k >= 2:
        rows = ([{"term": "group", "ss": ss_group, "df1": k - 1, "ss_err": ss_err_b,
                  "df2": df_err_b, "within": False}] + rows +
                [{"term": "interaction", "ss": ss_inter, "df1": (k - 1) * (t - 1),
                  "ss_err": ss_err_w, "df2": df_err_w, "within": True}])
    return {"rows": rows, "eps": eps, "mauchly": mauchly(m, p, df_err_b), "cells": cells}


def effect_row(r, eps, use_gg, labels):
    """Turn one effect into a stats row.

    Args:
        r: The effect with its sums of squares and degrees of freedom.
        eps: The Greenhouse-Geisser epsilon.
        use_gg: Whether the within-subject effects take the correction.
        labels: The display name of each term.

    Returns:
        dict: The stats row.
    """
    fval = (r["ss"] / r["df1"]) / (r["ss_err"] / r["df2"])
    p = float(stats.f.sf(fval, r["df1"], r["df2"]))
    p_gg = float(stats.f.sf(fval, eps * r["df1"], eps * r["df2"])) if r["within"] else np.nan
    corrected = r["within"] and use_gg
    return {"effect": labels[r["term"]], "ss": float(r["ss"]), "df1": r["df1"], "df2": r["df2"],
            "F": float(fval), "p_value": p, "gg_epsilon": eps if r["within"] else np.nan,
            "p_gg": p_gg, "correction": "Greenhouse-Geisser" if corrected else "none",
            "p_reported": p_gg if corrected else p,
            "partial_eta2": float(r["ss"] / (r["ss"] + r["ss_err"]))}


def recipe_timecourse(df, spec):
    """Draw a repeated-measures time course with a mixed ANOVA.

    Args:
        df: The tidy table, one row per subject and time.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a subject has two values at one
            time, if a subject is in two groups, or if too few subjects remain.
    """
    dat = spec["data"]
    x, y, idc, grp = dat.get("x"), dat.get("y"), dat.get("id"), dat.get("group")
    if not x or not y:
        raise ValueError("timecourse needs --x (the time) and --y (the measured value)")
    if not idc:
        raise ValueError("timecourse needs --id (the subject column)")
    has_group = bool(grp)
    for c in [idc, x, y] + ([grp] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({"id": df[idc].astype(str).where(df[idc].notna()), "time": df[x],
                      "y": pd.to_numeric(df[y], errors="coerce"),
                      "group": df[grp].astype(str).where(df[grp].notna()) if has_group else "all"})
    d = d.dropna().reset_index(drop=True)
    n_missing = n0 - len(d)
    if d.duplicated(["id", "time"]).any():
        raise ValueError("timecourse needs one value per subject and time; average any technical "
                         "replicates first")
    if (d.groupby("id")["group"].nunique() > 1).any():
        raise ValueError("each subject must belong to one group only")
    time_numeric = pd.api.types.is_numeric_dtype(d["time"])
    t_lv = sorted(d["time"].unique()) if time_numeric else list(pd.unique(d["time"].astype(str)))
    d["time_key"] = d["time"] if time_numeric else d["time"].astype(str)
    counts = d.groupby("id").size()
    complete = sorted(counts.index[counts == len(t_lv)])
    incomplete = int((counts != len(t_lv)).sum())
    d = d[d["id"].isin(complete)].reset_index(drop=True)
    g_lv = list(pd.unique(d["group"]))
    steps = []
    if n_missing:
        steps.append(f"Dropped {n_missing} row(s) with a missing subject, time, value or group.")
    if incomplete:
        steps.append(f"Dropped {incomplete} subject(s) without a value at every one of the "
                     f"{len(t_lv)} time points, because the repeated-measures ANOVA needs "
                     "complete subjects.")
    steps.append(f"Kept {len(complete)} subject(s) at {len(t_lv)} time points ({len(d)} values).")
    if len(t_lv) < 2:
        raise ValueError("timecourse needs at least 2 time points")
    if len(complete) - len(g_lv) < 1:
        raise ValueError("timecourse needs more complete subjects than groups")

    wide = d.pivot(index="id", columns="time_key", values="y").loc[complete, t_lv]
    ymat = wide.to_numpy(float)
    g = d.drop_duplicates("id").set_index("id").loc[complete, "group"].to_numpy()
    res = mixed_anova(ymat, g, g_lv)
    method = (spec["test"].get("method") or "auto").lower()
    if method not in ("auto", "gg", "uncorrected"):
        raise ValueError(f"unknown timecourse test '{method}' (use auto, gg or uncorrected)")
    use_gg = {"gg": True, "uncorrected": False,
              "auto": len(t_lv) > 2 and res["mauchly"]["p_value"] < 0.05}[method]
    labels = {"group": grp or "group", "time": x, "interaction": f"{grp or 'group'}:{x}"}
    stats_df = pd.DataFrame([effect_row(r, res["eps"], use_gg, labels) for r in res["rows"]])

    fig, w, h = _plot(ymat, g, g_lv, t_lv, time_numeric, stats_df, has_group, spec, x, y, grp)
    resid = ymat - res["cells"][[g_lv.index(v) for v in g]]
    name = "mixed ANOVA (group between, time within)" if has_group else "one-way repeated-measures ANOVA"
    test_meta = {
        "name": name, "n_subjects": len(complete), "n_times": len(t_lv),
        "n_dropped_subjects": incomplete,
        "sphericity": {"test": "Mauchly", **res["mauchly"]},
        "gg_epsilon": res["eps"], "correction_applied": use_gg,
        "effects": {r["effect"]: {k: r[k] for k in ("F", "df1", "df2", "p_value", "p_gg",
                                                     "p_reported", "partial_eta2")}
                    for _, r in stats_df.iterrows()},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, res, use_gg, method, len(complete), len(t_lv),
                            incomplete, has_group, spec, y),
        "df_used": d[["id", "group", "time", "y"]], "clean_steps": steps,
        "label": slugify(y), "width": w, "height": h,
        "qc": {"quantity": "residuals from the group-by-time means",
               "panels": {"residuals": resid.ravel()}},
        "build_script": functools.partial(_script, spec, has_group),
    }


def _plot(ymat, g, g_lv, t_lv, time_numeric, stats_df, has_group, spec, x, y, grp):
    """Draw the group means +/- SEM over time with thin subject lines.

    Args:
        ymat: The subjects-by-times matrix.
        g: The group of each row.
        g_lv: The group levels, in order.
        t_lv: The time levels, in order.
        time_numeric: Whether the time is numeric.
        stats_df: The ANOVA table.
        has_group: Whether there is a group column.
        spec: The figure specification.
        x: The time column.
        y: The value column.
        grp: The group column.

    Returns:
        tuple: The figure, its width and its height.
    """
    ap = spec["appearance"]
    pal = palette(len(g_lv), ap.get("palette"))
    pos = np.asarray(t_lv, dtype=float) if time_numeric else np.arange(1, len(t_lv) + 1, dtype=float)
    w, h = (6.0 if has_group else 5.0), 4.4
    fig, ax = new_fig(w, h)
    span = pos.max() - pos.min()
    for i, lv in enumerate(g_lv):
        yg = ymat[g == lv]
        if ap.get("show_points", True):
            for row in yg:
                ax.plot(pos, row, color=pal[i], alpha=0.18, lw=0.5)
        sem = yg.std(axis=0, ddof=1) / np.sqrt(len(yg)) if len(yg) > 1 else np.full(len(pos), np.nan)
        off = (i - (len(g_lv) - 1) / 2) * 0.012 * span if len(g_lv) > 1 else 0.0
        ax.errorbar(pos + off, yg.mean(axis=0), yerr=sem, color=pal[i], lw=1.6, marker="o",
                    ms=5, elinewidth=1.0, capsize=0, label=str(lv))
    sub = "; ".join(f"{r['effect']}: {p_label(r['p_reported'])}"
                    f"{'' if r['correction'] == 'none' else ' (GG)'}" for _, r in stats_df.iterrows())
    ax.set_title(sub, fontsize=9.5, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    if not time_numeric:
        ax.set_xticks(pos)
        ax.set_xticklabels([str(t) for t in t_lv])
    elif len(pos) <= 8:
        ax.set_xticks(pos)
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if has_group:
        ax.legend(title=grp, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()
    return fig, w, h


def _df_label(df, eps, correction):
    """Format a degree of freedom, scaled by epsilon when corrected.

    Args:
        df: The degrees of freedom.
        eps: The Greenhouse-Geisser epsilon.
        correction: The correction name, or "none".

    Returns:
        str: The label.
    """
    return f"{df:g}" if correction == "none" else f"{df * eps:.2f}"


def _methods(rows, res, use_gg, method, n_subj, n_time, n_drop, has_group, spec, y):
    """Write the methods paragraph.

    Args:
        rows: The ANOVA table.
        res: The ANOVA result with epsilon and Mauchly's test.
        use_gg: Whether the correction was applied.
        method: The requested test mode.
        n_subj: The count of complete subjects.
        n_time: The count of time points.
        n_drop: The count of dropped subjects.
        has_group: Whether there is a group column.
        spec: The figure specification.
        y: The value column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    design = ("a mixed ANOVA with the group as a between-subjects factor and time as a "
              "within-subjects factor") if has_group else "a one-way repeated-measures ANOVA over time"
    drop_txt = (f" {n_drop} subject(s) without a value at every time point were excluded."
                if n_drop else "")
    mt = res["mauchly"]
    if method == "auto":
        corr_txt = (f" Sphericity was assessed with Mauchly's test (W = {mt['W']:.3f}, "
                    f"P = {mt['p_value']:.2g}); " +
                    (f"because it was rejected, the within-subject effects carry the "
                     f"Greenhouse-Geisser correction (epsilon = {res['eps']:.3f})." if use_gg
                     else "it was not rejected, so no correction was applied."))
    elif use_gg:
        corr_txt = (f" The within-subject effects carry the Greenhouse-Geisser correction "
                    f"(epsilon = {res['eps']:.3f}).")
    else:
        corr_txt = " No sphericity correction was applied."
    eff = "; ".join(
        f"{r['effect']}, F({_df_label(r['df1'], r['gg_epsilon'], r['correction'])}, "
        f"{_df_label(r['df2'], r['gg_epsilon'], r['correction'])}) = {r['F']:.2f}, "
        f"{p_label(r['p_reported'])}, partial eta squared = {r['partial_eta2']:.3f}"
        for _, r in rows.iterrows())
    return (
        f"{cap_first(y_lab)} was analysed with {design} (n = {n_subj} subjects, {n_time} time "
        f"points).{drop_txt}{corr_txt} {eff}. The figure shows the mean +/- SEM at each time "
        "point, with one thin line per subject. Computed in Python "
        f"{platform.python_version()} with numpy {v('numpy')} and scipy {v('scipy')}, and drawn "
        f"with matplotlib {v('matplotlib')}."
    )


def _script(spec, has_group, in_name, fig_stub):
    """Write the standalone script that repeats the ANOVA with pingouin.

    Args:
        spec: The figure specification.
        has_group: Whether there is a group column.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    d = spec["data"]
    x, y, idc, grp = d["x"], d["y"], d["id"], d.get("group")
    anova = (f'pg.mixed_anova(data=raw, dv="{y}", within="{x}", subject="{idc}", between="{grp}")'
             if has_group else f'pg.rm_anova(data=raw, dv="{y}", within="{x}", subject="{idc}")')
    hue = f'"{grp}"' if has_group else "None"
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The ANOVA
# uses the subjects seen at every time point.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, pingouin as pg, seaborn as sns

raw = pd.read_csv("{in_name}").dropna(subset=["{idc}", "{x}", "{y}"])
n_times = raw["{x}"].nunique()
raw = raw[raw.groupby("{idc}")["{x}"].transform("count") == n_times]
print({anova})
fig, ax = plt.subplots(figsize=(6.0, 4.4))
sns.lineplot(data=raw, x="{x}", y="{y}", hue={hue}, errorbar="se", marker="o", ax=ax)
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
