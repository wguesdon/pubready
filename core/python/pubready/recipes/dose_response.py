"""Recipe: dose-response curves with a four-parameter logistic fit.

Mirrors the R recipe. Each curve is fitted by least squares to

    y = bottom + (top - bottom) / (1 + 10**((log10(EC50) - log10(dose)) * hill)),

the variable-slope model of GraphPad Prism. A positive Hill slope rises with the
dose and a negative one falls, so a falling curve reports its midpoint as an
IC50. The midpoint carries a 95% CI from the standard error of log10(EC50). With
two or more curves (--group), an extra-sum-of-squares F test asks whether one
shared EC50 fits the data as well as one EC50 per curve.

Both engines compute the same starting values and minimise the same residual sum
of squares, so they converge to the same estimates.
"""
import functools

import numpy as np
import pandas as pd
from scipy import optimize, stats

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette, plain_log_axis
from ..util import cap_first

_FIT_OPTS = {"maxfev": 20000, "ftol": 1e-12, "xtol": 1e-12, "gtol": 1e-12}


def fpl(lx, bottom, top, log_ec50, hill):
    """Evaluate the four-parameter logistic on log10 dose.

    Args:
        lx: The log10 doses.
        bottom: The lower plateau.
        top: The upper plateau.
        log_ec50: The log10 of the midpoint dose.
        hill: The Hill slope.

    Returns:
        The predicted responses.
    """
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((log_ec50 - lx) * hill))


def start_values(lx, y):
    """Compute the starting values, as the R engine does.

    The plateaus come from the lowest and highest mean response, the midpoint
    from the dose whose mean is nearest half way, and the Hill slope is +1 or -1
    from the sign of the correlation between log dose and response.

    Args:
        lx: The log10 doses.
        y: The responses.

    Returns:
        dict: bottom, top, log_ec50 and hill.
    """
    ux = np.unique(lx)
    means = np.array([y[lx == u].mean() for u in ux])
    bottom, top = float(means.min()), float(means.max())
    hill = -1.0 if np.corrcoef(lx, y)[0, 1] < 0 else 1.0
    log_ec50 = float(ux[int(np.argmin(np.abs(means - (bottom + top) / 2)))])
    return {"bottom": bottom, "top": top, "log_ec50": log_ec50, "hill": hill}


def fit_curve(lx, y):
    """Fit one curve, with Levenberg-Marquardt first and then a trust region.

    Args:
        lx: The log10 doses.
        y: The responses.

    Returns:
        tuple: The parameter vector, its covariance, the starting values and the
        routine that converged, or (None, None, start, None) when neither does.
    """
    st = start_values(lx, y)
    p0 = [st["bottom"], st["top"], st["log_ec50"], st["hill"]]
    for method in ("lm", "trf"):
        try:
            popt, pcov = optimize.curve_fit(fpl, lx, y, p0=p0, method=method, **_FIT_OPTS)
        except (RuntimeError, ValueError, optimize.OptimizeWarning):
            continue
        if np.all(np.isfinite(pcov)):
            return popt, pcov, st, method
    return None, None, st, None


def curve_row(g, lx, y, popt, pcov):
    """Summarise one fitted curve as one row of the stats table.

    Args:
        g: The curve name.
        lx: The log10 doses.
        y: The responses.
        popt: The fitted bottom, top, log10 EC50 and Hill slope.
        pcov: Their covariance.

    Returns:
        dict: The stats row.
    """
    bottom, top, log_ec50, hill = (float(v) for v in popt)
    se = np.sqrt(np.diag(pcov))
    # The plateaus swap when the slope changes sign; report bottom below top.
    if bottom > top:
        bottom, top, hill = top, bottom, -hill
    resid = y - fpl(lx, *popt)
    rss = float(np.sum(resid ** 2))
    dfr = len(y) - 4
    tq = stats.t.ppf(0.975, dfr)
    return {
        "group": g, "n": int(len(y)), "doses": int(len(np.unique(lx))),
        "midpoint": "IC50" if hill < 0 else "EC50",
        "bottom": bottom, "top": top, "hill": hill,
        "log_ec50": log_ec50, "se_log_ec50": float(se[2]),
        "ec50": 10 ** log_ec50,
        "ec50_ci_low": 10 ** (log_ec50 - tq * se[2]),
        "ec50_ci_high": 10 ** (log_ec50 + tq * se[2]),
        "r_squared": 1 - rss / float(np.sum((y - y.mean()) ** 2)),
        "rss": rss, "df_resid": dfr,
    }


def rss_fixed(u, sgn, lx, y, log_ec50):
    """Residual sum of squares of one curve at a fixed log10 EC50 and Hill slope.

    With both fixed, the 4PL is linear in the two plateaus, so they come from an
    exact least-squares solve (variable projection), as in the R engine.

    Args:
        u: The log of the absolute Hill slope.
        sgn: The sign of the Hill slope.
        lx: The log10 doses.
        y: The responses.
        log_ec50: The fixed log10 EC50.

    Returns:
        float: The residual sum of squares.
    """
    s = 1.0 / (1.0 + 10.0 ** ((log_ec50 - lx) * sgn * np.exp(u)))
    if s.max() - s.min() < 1e-10:
        return float(np.sum((y - y.mean()) ** 2))
    x = np.column_stack([1 - s, s])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    return float(np.sum((y - x @ coef) ** 2))


def rss_at(log_ec50, d, rows):
    """Residual sum of squares of the shared model at one fixed log10 EC50.

    Each curve keeps its own plateaus and slope. The slope comes from Brent's
    method over log(|hill|) in [log 0.05, log 20], with the sign of the separate fit.

    Args:
        log_ec50: The shared log10 EC50.
        d: The cleaned data with group, lx and y columns.
        rows: The stats table of the separate fits.

    Returns:
        float: The summed residual sum of squares.
    """
    total = 0.0
    for _, r in rows.iterrows():
        sel = (d["group"] == r["group"]).to_numpy()
        lx, y = d.loc[sel, "lx"].to_numpy(), d.loc[sel, "y"].to_numpy()
        sgn = -1.0 if r["hill"] < 0 else 1.0
        opt = optimize.minimize_scalar(rss_fixed, bounds=(np.log(0.05), np.log(20)),
                                       args=(sgn, lx, y, log_ec50), method="bounded",
                                       options={"xatol": 1e-10})
        total += float(opt.fun)
    return total


def shared_test(d, rows):
    """Extra-sum-of-squares F test of one shared log10 EC50, as the R engine does.

    The shared fit is a profile: the residual sum of squares over a grid of 31
    log10 EC50 values across the doses, then Brent's method in the interval around
    the best grid point. A joint fit from one start can stop in a local minimum
    when the curves sit far apart; the profile finds the global one.

    Args:
        d: The cleaned data with group, lx and y columns.
        rows: The stats table of the separate fits.

    Returns:
        dict or None: The F test, or None when no grid point converges.
    """
    grid = np.linspace(d["lx"].min(), d["lx"].max(), 31)
    prof = np.array([rss_at(g, d, rows) for g in grid])
    if not np.isfinite(prof).any():
        return None
    i = int(np.argmin(prof))
    bounds = (grid[max(0, i - 1)], grid[min(len(grid) - 1, i + 1)])
    opt = optimize.minimize_scalar(rss_at, bounds=bounds, args=(d, rows), method="bounded",
                                   options={"xatol": 1e-10})
    if opt.fun <= prof[i]:
        rss_sh, log_sh = float(opt.fun), float(opt.x)
    else:
        rss_sh, log_sh = float(prof[i]), float(grid[i])
    rss_sep = float(rows["rss"].sum())
    df_sep = int(rows["df_resid"].sum())
    df1 = len(rows) - 1
    fval = ((rss_sh - rss_sep) / df1) / (rss_sep / df_sep)
    return {"statistic": fval, "df1": int(df1), "df2": df_sep,
            "p_value": float(stats.f.sf(fval, df1, df_sep)),
            "shared_ec50": float(10 ** log_sh),
            "rss_shared": rss_sh, "rss_separate": rss_sep}


def dose_summary(d):
    """Mean and standard error of the response at each dose of each curve.

    Args:
        d: The cleaned data with group, lx and y columns.

    Returns:
        pandas.DataFrame: One row per curve and dose.
    """
    out = d.groupby(["group", "lx"], sort=False)["y"].agg(["mean", "count", "std"]).reset_index()
    out["sem"] = np.where(out["count"] > 1, out["std"] / np.sqrt(out["count"]), np.nan)
    return out.rename(columns={"count": "n"})


def recipe_dose_response(df, spec):
    """Fit and draw four-parameter logistic dose-response curves.

    Args:
        df: The tidy table, one row per measurement.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a curve has too few doses, or if a
            fit does not converge.
    """
    x, y = spec["data"].get("x"), spec["data"].get("y")
    grp = spec["data"].get("group")
    if not x or not y:
        raise ValueError("dose_response needs --x (the dose or concentration) and --y (the response)")
    has_group = bool(grp)
    for c in [x, y] + ([grp] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({
        "dose": pd.to_numeric(df[x], errors="coerce"),
        "y": pd.to_numeric(df[y], errors="coerce"),
        "group": df[grp].astype(str).where(df[grp].notna()) if has_group else "all",
    })
    ok = np.isfinite(d["dose"]) & np.isfinite(d["y"]) & d["group"].notna()
    n_missing = int((~ok).sum())
    d = d[ok]
    nonpos = d["dose"] <= 0
    d = d[~nonpos].reset_index(drop=True)
    steps = []
    if n_missing:
        steps.append(f"Dropped {n_missing} row(s) with a missing or non-numeric dose, response or group.")
    if int(nonpos.sum()):
        steps.append(f"Dropped {int(nonpos.sum())} row(s) with a dose of zero or below, "
                     "because the dose axis is logarithmic.")
    steps.append(f"Kept {len(d)} of {n0} row(s).")
    d["lx"] = np.log10(d["dose"])

    groups = list(pd.unique(d["group"]))
    rows, starts, methods = [], {}, {}
    for g in groups:
        sel = d["group"] == g
        lx, yy = d.loc[sel, "lx"].to_numpy(), d.loc[sel, "y"].to_numpy()
        if len(np.unique(lx)) < 4 or len(yy) < 5:
            raise ValueError(f"dose_response needs at least 4 distinct doses and 5 points per curve; "
                             f"curve '{g}' has {len(np.unique(lx))} dose(s) and {len(yy)} point(s)")
        popt, pcov, st, method = fit_curve(lx, yy)
        if popt is None:
            raise ValueError(f"the four-parameter logistic did not converge for curve '{g}'")
        rows.append(curve_row(g, lx, yy, popt, pcov))
        starts[g] = st
        methods[g] = method
    stats_df = pd.DataFrame(rows)
    shared = shared_test(d, stats_df) if len(groups) >= 2 else None

    ap = spec["appearance"]
    pal = palette(len(groups), ap.get("palette"))
    mid_name = "IC50" if (stats_df["midpoint"] == "IC50").all() else "EC50"
    w, h = (6.2 if has_group else 4.6), 4.2
    fig, ax = new_fig(w, h)
    summ = dose_summary(d)
    replicated = bool((summ["n"] > 1).any())
    for i, g in enumerate(groups):
        r = stats_df.iloc[i]
        lxg = d.loc[d["group"] == g, "lx"]
        grid = np.linspace(lxg.min(), lxg.max(), 200)
        cf = (r["bottom"], r["top"], r["log_ec50"], r["hill"])
        label = f"{g} ({r['midpoint']} = {r['ec50']:.3g})"
        ax.plot(10 ** grid, fpl(grid, *cf), color=pal[i], lw=1.6, label=label)
        s = summ[summ["group"] == g]
        if replicated:
            rep = s[s["n"] > 1]
            ax.errorbar(10 ** rep["lx"], rep["mean"], yerr=rep["sem"], fmt="none",
                        ecolor=pal[i], elinewidth=0.8, capsize=2)
        ax.scatter(10 ** s["lx"], s["mean"], color=pal[i], s=22, zorder=3)
    ax.set_xscale("log")
    plain_log_axis(ax.xaxis)
    if shared is not None:
        sub = (f"Shared {mid_name} F test: F({shared['df1']}, {shared['df2']}) = "
               f"{shared['statistic']:.2f}, P = {shared['p_value']:.2g}")
    else:
        r = stats_df.iloc[0]
        sub = (f"{r['midpoint']} = {r['ec50']:.3g} (95% CI {r['ec50_ci_low']:.3g} to "
               f"{r['ec50_ci_high']:.3g})")
    ax.set_title(sub, fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if has_group:
        ax.legend(title=grp, loc="center left", bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=9)
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    test_meta = {
        "name": "Four-parameter logistic fit",
        "n_curves": len(groups),
        "curves": {r["group"]: {"midpoint": r["midpoint"], "ec50": r["ec50"],
                                "ci_low": r["ec50_ci_low"], "ci_high": r["ec50_ci_high"],
                                "hill": r["hill"], "bottom": r["bottom"], "top": r["top"],
                                "r_squared": r["r_squared"]} for r in rows},
        "shared_ec50_test": shared,
    }
    df_used = pd.DataFrame({"group": d["group"], "dose": d["dose"], "response": d["y"]})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, shared, len(d), spec, x, y),
        "df_used": df_used, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, starts, methods, has_group),
    }


def _methods(rows, shared, n_points, spec, x, y):
    """Write the methods paragraph.

    Args:
        rows: The stats table.
        shared: The shared-EC50 F test, or None.
        n_points: The count of points fitted.
        spec: The figure specification.
        x: The dose column.
        y: The response column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    x_lab = spec["appearance"].get("x_label") or x
    multi = len(rows) > 1
    fits = "; ".join(
        f"{(str(r['group']) + ': ') if multi else ''}{r['midpoint']} = {r['ec50']:.3g} "
        f"(95% CI {r['ec50_ci_low']:.3g} to {r['ec50_ci_high']:.3g}, Hill slope {r['hill']:.2f}, "
        f"R2 = {r['r_squared']:.3f})" for _, r in rows.iterrows())
    shared_txt = ""
    if shared is not None:
        mid = rows["midpoint"].iloc[0]
        shared_txt = (f" An extra-sum-of-squares F test compared one shared {mid} with one {mid} "
                      "per curve, with the shared fit profiled over the midpoint and each Hill "
                      "slope held between 0.05 and 20 in magnitude: "
                      f"F({shared['df1']}, {shared['df2']}) = {shared['statistic']:.2f}, "
                      f"P = {shared['p_value']:.2g}.")
    elif len(rows) > 1:
        shared_txt = (" The shared-midpoint F test could not be computed, because the shared fit "
                      "did not converge.")
    return (
        f"{cap_first(y_lab)} was modelled as a function of {x_lab} with a four-parameter logistic "
        "curve (variable Hill slope) fitted by nonlinear least squares on log10 dose "
        f"(n = {n_points} points, {len(rows)} curve(s)). The 95% confidence interval of the "
        "midpoint was computed on the log10 scale from the standard error of the fit and "
        f"back-transformed. {fits}.{shared_txt} Fitted in Python {platform.python_version()} with "
        f"scipy {v('scipy')} (curve_fit) and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, starts, methods, has_group, in_name, fig_stub):
    """Write the standalone script that refits and redraws the curves.

    Args:
        spec: The figure specification.
        starts: The starting values of each curve.
        methods: The curve_fit routine that converged for each curve.
        has_group: Whether the data has a curve column.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    x, y, grp = spec["data"]["x"], spec["data"]["y"], spec["data"].get("group")
    ap = spec["appearance"]
    group_expr = f'raw["{grp}"].astype(str)' if has_group else '"all"'
    start_txt = ",\n    ".join(
        f'"{g}": [{s["bottom"]!r}, {s["top"]!r}, {s["log_ec50"]!r}, {s["hill"]!r}]'
        for g, s in starts.items())
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import curve_fit


def fpl(lx, bottom, top, log_ec50, hill):
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - lx) * hill))


raw = pd.read_csv("{in_name}")
d = pd.DataFrame({{"dose": pd.to_numeric(raw["{x}"], errors="coerce"),
                   "y": pd.to_numeric(raw["{y}"], errors="coerce"),
                   "group": {group_expr}}})
d = d[np.isfinite(d["dose"]) & np.isfinite(d["y"]) & (d["dose"] > 0)]
d["lx"] = np.log10(d["dose"])
methods = {methods!r}
starts = {{
    {start_txt}
}}
fig, ax = plt.subplots(figsize=({6.2 if has_group else 4.6}, 4.2))
for g, p0 in starts.items():
    s = d[d["group"] == g]
    popt, pcov = curve_fit(fpl, s["lx"], s["y"], p0=p0, method=methods[g], maxfev=20000,
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    print(g, "bottom, top, log10 EC50, hill:", popt)
    grid = np.linspace(s["lx"].min(), s["lx"].max(), 200)
    ax.plot(10 ** grid, fpl(grid, *popt), label=g)
    ax.scatter(s["dose"], s["y"], s=14, alpha=0.6)
ax.set_xscale("log")
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
