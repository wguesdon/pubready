"""Recipe: enzyme kinetics with a Michaelis-Menten fit. Mirrors the R recipe.

The equation ``v = Vmax * S / (Km + S)`` is fitted to each curve (``--group``,
for example one curve per inhibitor concentration) by nonlinear least squares.
Each parameter carries its asymptotic standard error and the 95% CI
estimate +/- t(0.975, n - 2) * SE, the interval GraphPad Prism reports. With two
or more curves, two extra-sum-of-squares F tests compare one shared Km (a Vmax
per curve) and one shared Vmax (a Km per curve) with the separate fits.

Both engines start from the same values and minimise the same residual sum of
squares, so they converge to the same estimates.
"""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats

from ..bundle import slugify
from ..theme import apply_pub_style, palette
from ..util import cap_first, sig_text
from .estimation_plot import nice_ticks, wrap

_FIT_OPTS = {"maxfev": 20000, "ftol": 1e-12, "xtol": 1e-12, "gtol": 1e-12}


def mm(s, vmax, km):
    """Evaluate the Michaelis-Menten equation.

    Args:
        s: The substrate concentrations.
        vmax: The maximum rate.
        km: The Michaelis constant.

    Returns:
        The rates.
    """
    return vmax * s / (km + s)


def mm_jac(s, vmax, km):
    """The Jacobian of the Michaelis-Menten equation in (Vmax, Km).

    Args:
        s: The substrate concentrations.
        vmax: The maximum rate.
        km: The Michaelis constant.

    Returns:
        numpy.ndarray: One row per concentration, one column per parameter.
    """
    s = np.asarray(s, dtype=float)
    return np.column_stack([s / (km + s), -vmax * s / (km + s) ** 2])


def start_values(s, v):
    """The starting values, as the R engine computes them.

    Vmax0 is the largest rate, and Km0 is the positive substrate concentration
    whose mean rate is nearest Vmax0 / 2.

    Args:
        s: The substrate concentrations.
        v: The rates.

    Returns:
        dict: vmax and km.
    """
    us = np.unique(s[s > 0])
    means = np.array([v[s == u].mean() for u in us])
    vmax = float(v.max())
    return {"vmax": vmax, "km": float(us[int(np.argmin(np.abs(means - vmax / 2)))])}


def fit_curve(s, v):
    """Fit one curve, with Levenberg-Marquardt first and then a trust region.

    Args:
        s: The substrate concentrations.
        v: The rates.

    Returns:
        tuple: The estimates, their covariance, the starting values and the
        routine that converged, or (None, None, start, None).
    """
    st = start_values(s, v)
    for method in ("lm", "trf"):
        try:
            popt, pcov = optimize.curve_fit(mm, s, v, p0=[st["vmax"], st["km"]], jac=mm_jac,
                                            method=method, **_FIT_OPTS)
        except (RuntimeError, ValueError, optimize.OptimizeWarning):
            continue
        if np.all(np.isfinite(pcov)):
            return popt, pcov, st, method
    return None, None, st, None


def curve_row(g, s, v, popt, pcov):
    """Summarise one fitted curve as one row of the stats table.

    Args:
        g: The curve name.
        s: The substrate concentrations.
        v: The rates.
        popt: The fitted Vmax and Km.
        pcov: Their covariance.

    Returns:
        dict: The stats row.
    """
    vm, km = (float(p) for p in popt)
    vse, kse = (float(e) for e in np.sqrt(np.diag(pcov)))
    rss = float(np.sum((v - mm(s, vm, km)) ** 2))
    dfr = len(v) - 2
    tq = stats.t.ppf(0.975, dfr)
    row = {
        "group": g, "n": int(len(v)), "concentrations": int(len(np.unique(s))),
        "vmax": vm, "vmax_se": vse, "vmax_ci_low": vm - tq * vse, "vmax_ci_high": vm + tq * vse,
        "km": km, "km_se": kse, "km_ci_low": km - tq * kse, "km_ci_high": km + tq * kse,
        "residual_sd": float(np.sqrt(rss / dfr)), "r_squared": 1 - rss / float(np.sum((v - v.mean()) ** 2)),
        "rss": rss, "df_resid": int(dfr), "fit_status": "fitted", "status_reason": None,
    }
    reason = status_reason(row, float(np.max(s)))
    if reason is None:
        return row
    for k in FIT_COLUMNS:
        row[k] = float("nan")
    row["fit_status"] = "not determined"
    row["status_reason"] = reason
    return row


#: The columns that a fit that is not determined leaves empty, as dose_response does.
FIT_COLUMNS = ("vmax", "vmax_se", "vmax_ci_low", "vmax_ci_high", "km", "km_se", "km_ci_low",
               "km_ci_high", "residual_sd", "r_squared", "rss", "df_resid")


def status_reason(r, smax):
    """Say why a Michaelis-Menten fit is "not determined", as .enz_status_reason() in R.

    The rule (F50): Vmax or Km is not positive, Km lies above the highest
    substrate concentration, so the data do not reach half of Vmax, or the 95% CI
    of Km or of Vmax holds 0.

    Args:
        r: The stats row of the fit.
        smax: The highest substrate concentration of the curve.

    Returns:
        str or None: The reason, or None when the fit stands.
    """
    if not np.isfinite(r["vmax"]) or r["vmax"] <= 0:
        return "Vmax is not positive"
    if not np.isfinite(r["km"]) or r["km"] <= 0:
        return "Km is not positive"
    if r["km"] > smax:
        return "Km lies above the highest substrate concentration"
    if not np.isfinite(r["km_ci_low"]) or r["km_ci_low"] <= 0:
        return "the 95% CI of Km holds 0"
    if not np.isfinite(r["vmax_ci_low"]) or r["vmax_ci_low"] <= 0:
        return "the 95% CI of Vmax holds 0"
    return None


def shared_km_model(x, *p):
    """One Km and a Vmax per curve.

    Args:
        x: A two-row array of the concentrations and the curve indices.
        *p: The Vmax of each curve, then the shared Km.

    Returns:
        numpy.ndarray: The rates.
    """
    s, gi = x[0], x[1].astype(int)
    return np.asarray(p[:-1])[gi] * s / (p[-1] + s)


def shared_vmax_model(x, *p):
    """One Vmax and a Km per curve.

    Args:
        x: A two-row array of the concentrations and the curve indices.
        *p: The shared Vmax, then the Km of each curve.

    Returns:
        numpy.ndarray: The rates.
    """
    s, gi = x[0], x[1].astype(int)
    return p[0] * s / (np.asarray(p[1:])[gi] + s)


def shared_test(d, rows, shared):
    """Extra-sum-of-squares F test of one shared parameter, as the R engine does.

    The shared model starts from the separate estimates and the median of the
    shared parameter across curves.

    Args:
        d: The cleaned data with s, v and group columns.
        rows: The stats table of the separate fits.
        shared: "km" or "vmax".

    Returns:
        dict or None: The F test, or None when the shared fit does not converge.
    """
    gi = pd.Index(rows["group"]).get_indexer(d["group"])
    xdat = np.vstack([d["s"].to_numpy(float), gi.astype(float)])
    v = d["v"].to_numpy(float)
    if shared == "km":
        model, p0 = shared_km_model, list(rows["vmax"]) + [float(np.median(rows["km"]))]
    else:
        model, p0 = shared_vmax_model, [float(np.median(rows["vmax"]))] + list(rows["km"])
    popt = None
    for method in ("lm", "trf"):
        try:
            popt, _ = optimize.curve_fit(model, xdat, v, p0=p0, method=method, **_FIT_OPTS)
            break
        except (RuntimeError, ValueError, optimize.OptimizeWarning):
            popt = None
    if popt is None:
        return None
    rss_sh = float(np.sum((v - model(xdat, *popt)) ** 2))
    rss_sep = float(rows["rss"].sum())
    df_sep = int(rows["df_resid"].sum())
    df1 = len(rows) - 1
    fval = ((rss_sh - rss_sep) / df1) / (rss_sep / df_sep)
    return {"statistic": fval, "df1": int(df1), "df2": df_sep,
            "p_value": float(stats.f.sf(fval, df1, df_sep)),
            "shared_value": float(popt[-1] if shared == "km" else popt[0]),
            "rss_shared": rss_sh, "rss_separate": rss_sep}


def conc_summary(d):
    """Mean and standard error of the rate at each concentration of each curve.

    Args:
        d: The cleaned data with s, v and group columns.

    Returns:
        pandas.DataFrame: One row per curve and concentration.
    """
    out = d.groupby(["group", "s"], sort=False)["v"].agg(["mean", "count", "std"]).reset_index()
    out["sem"] = np.where(out["count"] > 1, out["std"] / np.sqrt(out["count"]), np.nan)
    return out.rename(columns={"count": "n"})


def fmt_test(name, t):
    """Format one F test for the subtitle.

    Args:
        name: The name of the test.
        t: The test, or None.

    Returns:
        str: The line.
    """
    if t is None:
        return f"{name}: not computed"
    p = "P < 0.001" if t["p_value"] < 0.001 else f"P = {sig_text(t['p_value'], 2)}"
    return f"{name}: F({t['df1']}, {t['df2']}) = {t['statistic']:.2f}, {p}"


def recipe_enzyme_kinetics(df, spec):
    """Fit and draw Michaelis-Menten curves with the shared-parameter F tests.

    Args:
        df: The tidy table, one row per measured rate.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a curve has too few
            concentrations, or if a fit does not converge.
    """
    x, y = spec["data"].get("x"), spec["data"].get("y")
    grp = spec["data"].get("group")
    if not x or not y:
        raise ValueError("enzyme_kinetics needs --x (the substrate concentration) and --y (the rate)")
    has_group = bool(grp)
    for c in [x, y] + ([grp] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({
        "s": pd.to_numeric(df[x], errors="coerce"),
        "v": pd.to_numeric(df[y], errors="coerce"),
        "group": df[grp].astype(str).where(df[grp].notna()) if has_group else "all",
    })
    ok = np.isfinite(d["s"]) & np.isfinite(d["v"]) & d["group"].notna()
    n_missing = int((~ok).sum())
    d = d[ok]
    neg = d["s"] < 0
    d = d[~neg].reset_index(drop=True)
    steps = []
    if n_missing:
        steps.append(f"Dropped {n_missing} row(s) with a missing or non-numeric substrate, rate or group.")
    if int(neg.sum()):
        steps.append(f"Dropped {int(neg.sum())} row(s) with a negative substrate concentration.")
    steps.append(f"Kept {len(d)} of {n0} row(s).")

    groups = list(pd.unique(d["group"]))
    rows, starts, methods = [], {}, {}
    for g in groups:
        sel = (d["group"] == g).to_numpy()
        s, v = d.loc[sel, "s"].to_numpy(float), d.loc[sel, "v"].to_numpy(float)
        n_conc = len(np.unique(s[s > 0]))
        if n_conc < 3 or len(v) < 3:
            raise ValueError("enzyme_kinetics needs at least 3 distinct positive substrate "
                             f"concentrations and 3 points per curve; curve '{g}' has {n_conc} and {len(v)}")
        popt, pcov, st, method = fit_curve(s, v)
        if popt is None:
            raise ValueError(f"the Michaelis-Menten fit did not converge for curve '{g}'")
        row = curve_row(g, s, v, popt, pcov)
        rows.append(row)
        if row["fit_status"] == "fitted":
            starts[g] = st
            methods[g] = method
    stats_df = pd.DataFrame(rows)
    det = (stats_df["fit_status"] == "fitted").to_numpy()
    # A curve that is not determined stays out of the shared-parameter tests (F50).
    shared_km = shared_vmax = None
    if det.sum() >= 2:
        d_det = d[d["group"].isin(stats_df.loc[det, "group"])].reset_index(drop=True)
        rows_det = stats_df[det].reset_index(drop=True)
        shared_km = shared_test(d_det, rows_det, "km")
        shared_vmax = shared_test(d_det, rows_det, "vmax")

    fig, w, h = _plot(d, groups, stats_df, shared_km, shared_vmax, spec, x, y, grp, has_group)
    test_meta = {
        "name": "Michaelis-Menten fit", "n_curves": len(groups),
        "curves": {r["group"]: {k: r[k] for k in ("fit_status", "status_reason", "vmax",
                                                   "vmax_ci_low", "vmax_ci_high", "km",
                                                   "km_ci_low", "km_ci_high", "r_squared")}
                   for r in rows},
        "shared_km_test": shared_km, "shared_vmax_test": shared_vmax,
    }
    df_used = pd.DataFrame({"group": d["group"], "substrate": d["s"], "rate": d["v"]})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, shared_km, shared_vmax, len(d), spec, x, y, grp,
                            bool((conc_summary(d)["n"] > 1).any())),
        "df_used": df_used, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, starts, methods, has_group),
    }


def _plot(d, groups, rows, shared_km, shared_vmax, spec, x, y, grp, has_group):
    """Draw the curves, the points and the tests.

    Args:
        d: The cleaned data.
        groups: The curve names in order.
        rows: The stats table.
        shared_km: The shared-Km test, or None.
        shared_vmax: The shared-Vmax test, or None.
        spec: The figure specification.
        x: The substrate column.
        y: The rate column.
        grp: The curve column, or None.
        has_group: Whether there is a curve column.

    Returns:
        tuple: The figure, its width and its height.
    """
    ap = spec["appearance"]
    pal = palette(len(groups), ap.get("palette"))
    labels = [f"{wrap(g, 18)}\nVmax {sig_text(r['vmax'])}, Km {sig_text(r['km'])}"
              if r["fit_status"] == "fitted" else f"{wrap(g, 18)}\nnot determined"
              for g, (_, r) in zip(groups, rows.iterrows())]
    key_lines = max(lab.count("\n") for lab in labels) + 1
    w, h = (6.6 if has_group else 4.8), max(4.4, 1.4 + 0.165 * len(groups) * key_lines)
    fig, ax = plt.subplots(figsize=(w, h))
    summ = conc_summary(d)
    replicated = bool((summ["n"] > 1).any())
    grid = np.linspace(0, d["s"].max(), 200)
    for i, g in enumerate(groups):
        r = rows.iloc[i]
        raw = d[d["group"] == g]
        s = summ[summ["group"] == g]
        if replicated:
            ax.scatter(raw["s"], raw["v"], s=6, color=pal[i], alpha=0.35, lw=0, zorder=2)
            rep = s[s["n"] > 1]
            ax.errorbar(rep["s"], rep["mean"], yerr=rep["sem"], fmt="none", ecolor=pal[i],
                        elinewidth=1.0, capsize=0, zorder=3)
        if r["fit_status"] == "fitted":
            ax.plot(grid, mm(grid, r["vmax"], r["km"]), color=pal[i], lw=1.6, label=labels[i],
                    zorder=4)
            ax.scatter(s["s"], s["mean"], s=24, color=pal[i], zorder=5)
        else:
            ax.scatter(s["s"], s["mean"], s=24, color=pal[i], zorder=5, label=labels[i])
    ax.set_title(subtitle(rows, shared_km, shared_vmax), fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    top = max(float(d["v"].max()), float(np.nanmax(summ["mean"] + summ["sem"].fillna(0))))
    ax.set_xlim(0, d["s"].max() * 1.04)
    ax.set_ylim(min(0.0, float(d["v"].min())), top * 1.06)
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if has_group:
        ax.legend(title=grp, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
                  fontsize=9.5, title_fontsize=11, labelspacing=1.0, alignment="left")
    apply_pub_style(ax, ap.get("theme"))
    nice_ticks(ax.xaxis)
    nice_ticks(ax.yaxis)
    fig.tight_layout()
    return fig, w, h


def subtitle(rows, shared_km, shared_vmax):
    """Write the subtitle, as .enz_subtitle() of the R engine does.

    The subtitle gives the two shared-parameter tests, or the estimates of the
    one curve, and a line that names the curves that are not determined (F50).

    Args:
        rows: The stats table.
        shared_km: The shared-Km test, or None.
        shared_vmax: The shared-Vmax test, or None.

    Returns:
        str: The subtitle.
    """
    det = (rows["fit_status"] == "fitted").to_numpy()
    if len(rows) == 1:
        r = rows.iloc[0]
        if not det[0]:
            return f"Fit not determined: {r['status_reason']}"
        return (f"Vmax = {sig_text(r['vmax'])} (95% CI {sig_text(r['vmax_ci_low'])} to "
                f"{sig_text(r['vmax_ci_high'])})\nKm = {sig_text(r['km'])} (95% CI "
                f"{sig_text(r['km_ci_low'])} to {sig_text(r['km_ci_high'])})")
    out = fmt_test("Shared Km", shared_km) + "\n" + fmt_test("Shared Vmax", shared_vmax)
    nd = [str(g) for g in rows.loc[~det, "group"]]
    if not nd:
        return out
    if len(nd) <= 2:
        return f"{out}\nNot determined, left out of the tests: {', '.join(nd)}"
    return f"{out}\n{len(nd)} curves not determined, left out of the tests"


def _test_text(name, t):
    """Write one F test for the methods paragraph.

    Args:
        name: The name of the test.
        t: The test, or None.

    Returns:
        str: The text.
    """
    if t is None:
        return f"{name} not computed, because the shared fit did not converge"
    return f"{name} F({t['df1']}, {t['df2']}) = {t['statistic']:.2f}, P = {sig_text(t['p_value'], 2)}"


def _points_text(replicated):
    """Say what each point of the figure shows.

    The figure draws the replicates and the error bars only when a
    concentration has replicates.

    Args:
        replicated: True when a concentration of a curve has more than one rate.

    Returns:
        str: The sentence, with a trailing space.
    """
    if replicated:
        return ("Large points show the mean +/- SEM of the replicates at each concentration, "
                "and small points show each replicate. ")
    return "Each point is the single rate at its concentration. "


def _methods(rows, shared_km, shared_vmax, n_points, spec, x, y, grp, replicated):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        rows: The stats table.
        shared_km: The shared-Km test, or None.
        shared_vmax: The shared-Vmax test, or None.
        n_points: The count of points fitted.
        spec: The figure specification.
        x: The substrate column.
        y: The rate column.
        grp: The curve column, or None.
        replicated: True when a concentration of a curve has more than one rate.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    x_lab = spec["appearance"].get("x_label") or x
    multi = len(rows) > 1
    det = rows["fit_status"] == "fitted"
    fits = "; ".join(
        f"{(str(r['group']) + ': ') if multi else ''}Vmax = {sig_text(r['vmax'])} (95% CI "
        f"{sig_text(r['vmax_ci_low'])} to {sig_text(r['vmax_ci_high'])}), Km = {sig_text(r['km'])} (95% CI "
        f"{sig_text(r['km_ci_low'])} to {sig_text(r['km_ci_high'])}), R2 = {r['r_squared']:.3f}"
        for _, r in rows[det].iterrows())
    fits = f"{fits}." if fits else ""
    nd = rows[~det]
    nd_txt = ""
    if len(nd):
        nd_txt = (" " if fits else "") + " ".join(
            f"The fit {('of ' + str(r['group']) + ' ') if multi else ''}was not determined, "
            f"because {r['status_reason']}." for _, r in nd.iterrows())
        nd_txt += (" A curve that is not determined is drawn as its points without a fit"
                   + (" and is left out of the F tests." if multi else "."))
    shared_txt = ""
    if multi and int(det.sum()) >= 2:
        shared_txt = (" Extra-sum-of-squares F tests compared the separate fits with one shared Km "
                      "(a Vmax for each curve) and with one shared Vmax (a Km for each curve): "
                      f"{_test_text('shared Km', shared_km)}; {_test_text('shared Vmax', shared_vmax)}.")
    elif multi:
        shared_txt = (" The F tests of a shared Km and a shared Vmax were not run, because fewer "
                      "than two curves were determined.")
    rule_txt = ("A fit was reported as not determined when Vmax or Km was not positive, when Km "
                "lay above the highest substrate concentration, or when the 95% CI of Km or of "
                "Vmax held 0. ")
    target = f"each level of {grp} ({len(rows)} curves)" if multi else "one curve"
    return (
        f"{cap_first(y_lab)} was modelled as a function of {x_lab} with the Michaelis-Menten "
        f"equation v = Vmax S / (Km + S), fitted by nonlinear least squares to {target} "
        f"(n = {n_points} points). The starting values were the largest rate for Vmax and, for "
        "Km, the substrate concentration whose mean rate is nearest half of it. Each 95% "
        "confidence interval is the estimate plus or minus t(0.975, n - 2) times its asymptotic "
        f"standard error. {_points_text(replicated)}{rule_txt}{fits}{nd_txt}{shared_txt} Fitted in Python "
        f"{platform.python_version()} with "
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
    start_txt = ",\n    ".join(f'"{g}": [{s["vmax"]!r}, {s["km"]!r}]' for g, s in starts.items())
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy.optimize import curve_fit


def mm(s, vmax, km):
    return vmax * s / (km + s)


raw = pd.read_csv("{in_name}")
d = pd.DataFrame({{"s": pd.to_numeric(raw["{x}"], errors="coerce"),
                   "v": pd.to_numeric(raw["{y}"], errors="coerce"),
                   "group": {group_expr}}})
d = d[np.isfinite(d["s"]) & np.isfinite(d["v"]) & (d["s"] >= 0)]
methods = {methods!r}
starts = {{
    {start_txt}
}}
fig, ax = plt.subplots(figsize=({6.6 if has_group else 4.8}, 4.4))
grid = np.linspace(0, d["s"].max(), 200)
for g, p0 in starts.items():
    s = d[d["group"] == g]
    popt, pcov = curve_fit(mm, s["s"], s["v"], p0=p0, method=methods[g], maxfev=20000,
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    print(g, "Vmax, Km:", popt, "SE:", np.sqrt(np.diag(pcov)))
    ax.plot(grid, mm(grid, *popt), label=g)
    ax.scatter(s["s"], s["v"], s=14, alpha=0.6)
# A curve that is not determined is drawn as points only, as in the figure.
for g in d["group"].unique():
    if g not in starts:
        s = d[d["group"] == g]
        ax.scatter(s["s"], s["v"], s=14, alpha=0.6, label=f"{{g}} (not determined)")
ax.set_xlim(left=0); ax.set_ylim(bottom=0)
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
