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

A curve whose mean responses vary by less than 10% of the response range, or
whose responses do not differ between the doses in a one-way ANOVA, is flat and
gets no fit (flat_check). A curve whose fit fails the rule of status_reason is
"not determined" too. Both are drawn as their mean +/- SEM without a curve.
"""
import functools

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import optimize, stats

from ..bundle import slugify
from ..theme import R_POINT_PT, apply_pub_style, new_fig, palette, plain_log_axis
from ..util import cap_first, sig_text

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


MIN_SPAN = 0.10
"""The share of the response range below which a fitted curve counts as flat."""

FLAT_REASON = "the mean responses vary by less than 10% of the data range"
"""The reason of a curve that is_flat finds flat, as .drc_flat_reason() in R."""


def is_flat(lx, y, span_ref):
    """Check, before the fit, whether a curve is flat, as .drc_is_flat() in R.

    A curve whose mean responses vary by less than 10% of the range of the mean
    responses of all curves has no midpoint to estimate. The check uses no fitted
    value, so both engines decide the same way whatever their optimisers do on
    flat data. A single curve is never flat by this rule, because its own range
    is the reference.

    Args:
        lx: The log10 doses of the curve.
        y: The responses of the curve.
        span_ref: The range of the mean responses of all curves.

    Returns:
        bool: Whether the curve is flat.
    """
    means = np.array([y[lx == u].mean() for u in np.unique(lx)])
    return bool(np.isfinite(span_ref) and span_ref > 0
                and means.max() - means.min() < MIN_SPAN * span_ref)


ANOVA_REASON = "the responses do not differ between the doses (one-way ANOVA, P >= 0.05)"
"""The reason of a curve that the dose ANOVA finds flat, as .drc_anova_reason() in R."""


def dose_anova_p(lx, y):
    """The P value of a one-way ANOVA of one curve across its doses.

    The sums of squares are written out as in .drc_dose_anova_p() of the R
    engine, so both engines compute the same value.

    Args:
        lx: The log10 doses of the curve.
        y: The responses of the curve.

    Returns:
        float: The P value, or NaN when no dose has a replicate, as the ANOVA
        then has no residual degree of freedom.
    """
    u = np.unique(lx)
    k, n = len(u), len(y)
    if k < 2 or n - k < 1:
        return float("nan")
    means = {v: y[lx == v].mean() for v in u}
    fitted = np.array([means[v] for v in lx])
    ssw = float(np.sum((y - fitted) ** 2))
    ssb = float(np.sum((fitted - y.mean()) ** 2))
    if ssw == 0:
        return 0.0 if ssb > 0 else float("nan")
    return float(stats.f.sf((ssb / (k - 1)) / (ssw / (n - k)), k - 1, n - k))


def flat_check(lx, y, span_ref):
    """Run the checks before the fit, as .drc_flat_check() in R.

    The checks, in order, are the 10% range rule of is_flat and the one-way
    ANOVA across the doses, which also covers a single curve.

    Args:
        lx: The log10 doses of the curve.
        y: The responses of the curve.
        span_ref: The range of the mean responses of all curves.

    Returns:
        str or None: The reason of a flat curve, or None when the curve goes to
        the fit.
    """
    if is_flat(lx, y, span_ref):
        return FLAT_REASON
    p = dose_anova_p(lx, y)
    if np.isfinite(p) and p >= 0.05:
        return ANOVA_REASON
    return None


def status_reason(converged, bottom, top, log_ec50, lx, ci_low, ci_high, span_ref):
    """Say why a curve is "not determined", as .drc_status_reason() of the R engine.

    The rule: the fit did not converge, the fitted span (top minus bottom) is
    below 10% of the range of the mean responses of all curves, the EC50 lies
    outside the tested doses, or the 95% CI of the EC50 is not finite.

    Args:
        converged: Whether the fit converged.
        bottom: The fitted lower plateau.
        top: The fitted upper plateau.
        log_ec50: The fitted log10 EC50.
        lx: The log10 doses of the curve.
        ci_low: The lower bound of the 95% CI of the EC50.
        ci_high: The upper bound.
        span_ref: The range of the mean responses of all curves.

    Returns:
        str or None: The reason, or None when the fit stands.
    """
    if not converged:
        return "the fit did not converge"
    if np.isfinite(span_ref) and span_ref > 0 and abs(top - bottom) < MIN_SPAN * span_ref:
        return "the response span is below 10% of the data range"
    if not np.isfinite(log_ec50) or log_ec50 < np.min(lx) or log_ec50 > np.max(lx):
        return "the EC50 lies outside the tested doses"
    if not np.isfinite(ci_low) or not np.isfinite(ci_high) or ci_low <= 0:
        return "the 95% CI of the EC50 is not finite"
    return None


def curve_row(g, lx, y, popt, pcov, start, span_ref, flat_reason=None):
    """Summarise one curve as one row of the stats table.

    A flat curve, or a curve whose fit fails the rule of status_reason, keeps its
    sizes and its direction but no fitted value, as in the R engine.

    Args:
        g: The curve name.
        lx: The log10 doses.
        y: The responses.
        popt: The fitted bottom, top, log10 EC50 and Hill slope, or None.
        pcov: Their covariance, or None.
        start: The starting values, which give the direction of a failed fit.
        span_ref: The range of the mean responses of all curves.
        flat_reason: The reason of flat_check when it found the curve flat, so that
            it was not fitted, or None.

    Returns:
        dict: The stats row.
    """
    nan = float("nan")
    base = {
        "group": g, "n": int(len(y)), "doses": int(len(np.unique(lx))),
        "midpoint_type": "IC50" if start["hill"] < 0 else "EC50",
        "bottom": nan, "top": nan, "hill": nan, "log_midpoint": nan, "se_log_midpoint": nan,
        "midpoint": nan, "midpoint_ci_low": nan, "midpoint_ci_high": nan, "r_squared": nan,
        "rss": nan, "df_resid": nan, "fit_status": "not determined",
        "status_reason": flat_reason or "the fit did not converge",
    }
    if popt is None:
        return base
    bottom, top, log_ec50, hill = (float(v) for v in popt)
    with np.errstate(invalid="ignore"):
        se = np.sqrt(np.diag(pcov))
    # The plateaus swap when the slope changes sign; report bottom below top.
    if bottom > top:
        bottom, top, hill = top, bottom, -hill
    resid = y - fpl(lx, *popt)
    rss = float(np.sum(resid ** 2))
    dfr = len(y) - 4
    tq = stats.t.ppf(0.975, dfr)
    with np.errstate(over="ignore", invalid="ignore"):
        ci_low = float(10 ** (log_ec50 - tq * se[2]))
        ci_high = float(10 ** (log_ec50 + tq * se[2]))
    reason = status_reason(True, bottom, top, log_ec50, lx, ci_low, ci_high, span_ref)
    if reason is not None:
        base["status_reason"] = reason
        return base
    return {
        "group": g, "n": int(len(y)), "doses": int(len(np.unique(lx))),
        "midpoint_type": "IC50" if hill < 0 else "EC50",
        "bottom": bottom, "top": top, "hill": hill,
        "log_midpoint": log_ec50, "se_log_midpoint": float(se[2]),
        "midpoint": 10 ** log_ec50, "midpoint_ci_low": ci_low, "midpoint_ci_high": ci_high,
        "r_squared": 1 - rss / float(np.sum((y - y.mean()) ** 2)),
        "rss": rss, "df_resid": dfr, "fit_status": "fitted", "status_reason": None,
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


def response_scale(y):
    """The power of ten that the response is divided by before the fit (F80).

    The scale is 10**e for a largest absolute response of the order 10**e when
    e is 3 or more or -3 or less, else 1, as .drc_response_scale() in R.

    Args:
        y: The responses.

    Returns:
        float: The scale.
    """
    m = float(np.max(np.abs(y)))
    if not np.isfinite(m) or m == 0:
        return 1.0
    e = int(np.floor(np.log10(m)))
    return float(10.0 ** e) if abs(e) >= 3 else 1.0


def unscale_row(r, rs):
    """Put the plateaus and the residual sum of squares of one row back on scale.

    Args:
        r: The stats row, changed in place.
        rs: The scale of the fit.
    """
    r["bottom"] = r["bottom"] * rs
    r["top"] = r["top"] * rs
    r["rss"] = r["rss"] * rs ** 2


def unscale_shared(shared, rs):
    """Put the sums of squares of the shared test back on the scale of the response.

    Args:
        shared: The shared-EC50 test, or None.
        rs: The scale of the fit.

    Returns:
        dict or None: The test.
    """
    if shared is None:
        return None
    shared["rss_shared"] = shared["rss_shared"] * rs ** 2
    shared["rss_separate"] = shared["rss_separate"] * rs ** 2
    return shared


def _scale_text(rs):
    """Say that the response was scaled for the fit, as .drc_scale_text() in R.

    Args:
        rs: The scale of the fit.

    Returns:
        str: The sentence with a trailing space, or "" for a scale of 1.
    """
    if rs == 1:
        return ""
    return (f"The responses were divided by {rs:g} for the fit, so that the optimiser works near "
            "1, and the plateaus are reported on the original scale. ")


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

    A curve whose fit fails, or whose midpoint is flat, outside the tested doses
    or without a finite CI, is reported as not determined (see status_reason):
    it is drawn as its mean +/- SEM without a curve and left out of the
    shared-midpoint F test.

    Args:
        df: The tidy table, one row per measurement.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, or if a curve has too few doses.
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
    # A response on a scale far from 1, such as 1e12 RLU, can keep a fit from
    # converging. The fit runs on the response divided by a power of ten, and the
    # plateaus and the sums of squares go back to the original scale (F80).
    rs = response_scale(d["y"].to_numpy())
    y_orig = d["y"].copy()
    d["y"] = d["y"] / rs

    groups = list(pd.unique(d["group"]))
    # The span rule compares each fitted span with the range of the mean
    # responses of all curves.
    summ = dose_summary(d)
    span_ref = float(summ["mean"].max() - summ["mean"].min())
    rows, starts, methods = [], {}, {}
    for g in groups:
        sel = d["group"] == g
        lx, yy = d.loc[sel, "lx"].to_numpy(), d.loc[sel, "y"].to_numpy()
        if len(np.unique(lx)) < 4 or len(yy) < 5:
            raise ValueError(f"dose_response needs at least 4 distinct doses and 5 points per curve; "
                             f"curve '{g}' has {len(np.unique(lx))} dose(s) and {len(yy)} point(s)")
        flat_reason = flat_check(lx, yy, span_ref)
        if flat_reason:
            popt, pcov, st, method = None, None, start_values(lx, yy), None
        else:
            popt, pcov, st, method = fit_curve(lx, yy)
        row = curve_row(g, lx, yy, popt, pcov, st, span_ref, flat_reason=flat_reason)
        rows.append(row)
        if row["fit_status"] == "fitted":
            starts[g] = st
            methods[g] = method
    stats_df = pd.DataFrame(rows)
    det = (stats_df["fit_status"] == "fitted").to_numpy()
    det_groups = [g for g, ok in zip(groups, det) if ok]
    nd_groups = [g for g, ok in zip(groups, det) if not ok]
    # A curve that is not determined stays out of the shared-midpoint test.
    shared = None
    if len(det_groups) >= 2:
        shared = shared_test(d[d["group"].isin(det_groups)].reset_index(drop=True),
                             stats_df[det].reset_index(drop=True))
    for r in rows:
        unscale_row(r, rs)
    stats_df = pd.DataFrame(rows)
    shared = unscale_shared(shared, rs)
    d["y"] = y_orig
    summ["mean"] = summ["mean"] * rs
    summ["sem"] = summ["sem"] * rs

    ap = spec["appearance"]
    pal = palette(len(groups), ap.get("palette"))
    mid_name = _mid_name(stats_df)
    w, h = (6.2 if has_group else 4.6), 4.2
    fig, ax = new_fig(w, h)
    replicated = bool((summ["n"] > 1).any())
    handles = []
    for i, g in enumerate(groups):
        r = stats_df.iloc[i]
        s = summ[summ["group"] == g]
        if replicated:
            rep = s[s["n"] > 1]
            ax.errorbar(10 ** rep["lx"], rep["mean"], yerr=rep["sem"], fmt="none",
                        ecolor=pal[i], elinewidth=0.8, capsize=2)
        # R draws the means with geom_point(size = 2).
        pts = ax.scatter(10 ** s["lx"], s["mean"], color=pal[i], s=(2 * R_POINT_PT) ** 2,
                         zorder=3, label=f"{g} (not determined)")
        if det[i]:
            lxg = d.loc[d["group"] == g, "lx"]
            grid = np.linspace(lxg.min(), lxg.max(), 200)
            cf = (r["bottom"], r["top"], r["log_midpoint"], r["hill"])
            ax.plot(10 ** grid, fpl(grid, *cf), color=pal[i], lw=1.6)
            # The key of a fitted curve shows the point on the line, as in R.
            handles.append(Line2D([], [], color=pal[i], lw=1.6, marker="o",
                                  markersize=2 * R_POINT_PT,
                                  label=f"{g} ({r['midpoint_type']} = {sig_text(r['midpoint'])})"))
        else:
            handles.append(pts)
    ax.set_xscale("log")
    plain_log_axis(ax.xaxis)
    ax.set_title(_subtitle(stats_df, shared, mid_name, has_group), fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if has_group:
        ax.legend(handles=handles, labels=[hd.get_label() for hd in handles], title=grp,
                  loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    test_meta = {
        "name": "Four-parameter logistic fit",
        "n_curves": len(groups),
        "n_determined": len(det_groups),
        "not_determined": nd_groups,
        "curves": {r["group"]: {"midpoint_type": r["midpoint_type"], "fit_status": r["fit_status"],
                                "status_reason": r["status_reason"], "midpoint": r["midpoint"],
                                "ci_low": r["midpoint_ci_low"], "ci_high": r["midpoint_ci_high"],
                                "hill": r["hill"], "bottom": r["bottom"], "top": r["top"],
                                "r_squared": r["r_squared"]} for r in rows},
        "shared_ec50_test": shared,
    }
    df_used = pd.DataFrame({"group": d["group"], "dose": d["dose"], "response": d["y"]})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, shared, len(d), spec, x, y, replicated,
                            int(nonpos.sum()), rs),
        "df_used": df_used, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, starts, methods, has_group, rs),
    }


def _mid_name(rows):
    """IC50 when every determined curve falls, EC50 otherwise.

    Args:
        rows: The stats table.

    Returns:
        str: "IC50" or "EC50".
    """
    mids = rows.loc[rows["fit_status"] == "fitted", "midpoint_type"]
    if len(mids) == 0:
        mids = rows["midpoint_type"]
    return "IC50" if (mids == "IC50").all() else "EC50"


def _subtitle(rows, shared, mid_name, has_group):
    """Write the subtitle, as .drc_subtitle() of the R engine does.

    The first line gives the shared-midpoint test, or the midpoint of the one
    determined curve; a second line names the curves left out as not determined.

    Args:
        rows: The stats table.
        shared: The shared-EC50 F test, or None.
        mid_name: "IC50" or "EC50".
        has_group: Whether the data has a curve column.

    Returns:
        str: The subtitle.
    """
    det = (rows["fit_status"] == "fitted").to_numpy()
    if shared is not None:
        first = (f"Shared {mid_name} F test: F({shared['df1']}, {shared['df2']}) = "
                 f"{shared['statistic']:.2f}, P = {sig_text(shared['p_value'], 2)}")
    elif det.any():
        r = rows[det].iloc[0]
        first = (f"{(str(r['group']) + ' ') if has_group else ''}{r['midpoint_type']} = {sig_text(r['midpoint'])} "
                 f"(95% CI {sig_text(r['midpoint_ci_low'])} to {sig_text(r['midpoint_ci_high'])})")
    elif has_group:
        first = f"{mid_name} not determined for any curve"
    else:
        first = f"{mid_name} not determined: {rows['status_reason'].iloc[0]}"
    nd = list(rows.loc[~det, "group"])
    if not has_group or not nd or not det.any():
        return first
    if len(nd) <= 2:
        second = f"Not determined, left out of the test: {', '.join(str(g) for g in nd)}"
    else:
        second = f"{len(nd)} curves not determined, left out of the test"
    return f"{first}\n{second}"


def _nonpos_text(n):
    """Say how many rows with a dose of zero or below were left out (F76).

    Args:
        n: The count of such rows.

    Returns:
        str: The sentence with a trailing space, or "" when no row was left out.
    """
    if n > 0:
        return (f"{n} row(s) with a dose of zero or below were left out of the fit and the "
                "figure, because the dose axis is logarithmic. ")
    return ""


def _points_text(replicated):
    """Say what each point of the figure shows.

    The figure draws error bars only when a dose has replicates.

    Args:
        replicated: True when a dose of a curve has more than one response.

    Returns:
        str: The sentence, with a trailing space.
    """
    if replicated:
        return "Points show the mean +/- SEM of the replicates at each dose. "
    return "Each point is the single response at its dose. "


def _methods(rows, shared, n_points, spec, x, y, replicated, n_nonpos=0, rs=1.0):
    """Write the methods paragraph.

    Args:
        rows: The stats table.
        shared: The shared-EC50 F test, or None.
        n_points: The count of points fitted.
        spec: The figure specification.
        x: The dose column.
        y: The response column.
        replicated: True when a dose of a curve has more than one response.
        n_nonpos: The count of rows dropped for a dose of zero or below.
        rs: The power of ten that the response was divided by for the fit.

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
        f"{(str(r['group']) + ': ') if multi else ''}{r['midpoint_type']} = {sig_text(r['midpoint'])} "
        f"(95% CI {sig_text(r['midpoint_ci_low'])} to {sig_text(r['midpoint_ci_high'])}, Hill slope {r['hill']:.2f}, "
        f"R2 = {r['r_squared']:.3f})" for _, r in rows[det].iterrows())
    fits = f"{fits}." if fits else ""
    nd = rows[~det]
    nd_txt = ""
    if len(nd):
        nd_txt = " " + " ".join(
            f"The {r['midpoint_type']}{(' of ' + str(r['group'])) if multi else ''} was not determined, "
            f"because {r['status_reason']}." for _, r in nd.iterrows())
        nd_txt += (" A curve that is not determined is drawn as its mean +/- SEM without a fit"
                   + (" and is left out of the shared-midpoint F test." if multi else "."))
    shared_txt = ""
    if shared is not None:
        mid = _mid_name(rows)
        shared_txt = (f" An extra-sum-of-squares F test compared one shared {mid} with one {mid} "
                      f"per curve over the {int(det.sum())} determined curves, with the shared fit "
                      "profiled over the midpoint and each Hill slope held between 0.05 and 20 in "
                      f"magnitude: F({shared['df1']}, {shared['df2']}) = {shared['statistic']:.2f}, "
                      f"P = {sig_text(shared['p_value'], 2)}.")
    elif multi and int(det.sum()) < 2:
        shared_txt = (" The shared-midpoint F test was not run, because fewer than two curves "
                      "were determined.")
    elif multi:
        shared_txt = (" The shared-midpoint F test could not be computed, because the shared fit "
                      "did not converge.")
    return (
        f"{cap_first(y_lab)} was modelled as a function of {x_lab} with a four-parameter logistic "
        "curve (variable Hill slope) fitted by nonlinear least squares on log10 dose "
        f"(n = {n_points} points, {len(rows)} curve(s)). The 95% confidence interval of the "
        "midpoint was computed on the log10 scale from the standard error of the fit and "
        f"back-transformed. {_points_text(replicated)}{_nonpos_text(n_nonpos)}{_scale_text(rs)}A curve whose mean "
        "responses varied by "
        "less than 10% of the range "
        "of the mean responses of all curves, or whose responses did not differ between the "
        "doses in a one-way ANOVA (P >= 0.05), was not fitted, and its midpoint was reported as "
        "not determined. A midpoint was also reported as not determined when the fit did not "
        "converge, when the fitted span (top minus bottom) was below 10% of the range of the mean "
        "responses, when the midpoint lay outside the tested doses, or when its 95% CI was not "
        f"finite. {fits}{nd_txt}{shared_txt} Fitted in Python {platform.python_version()} with "
        f"scipy {v('scipy')} (curve_fit) and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, starts, methods, has_group, rs, in_name, fig_stub):
    """Write the standalone script that refits and redraws the curves.

    Args:
        spec: The figure specification.
        starts: The starting values of each curve.
        methods: The curve_fit routine that converged for each curve.
        has_group: Whether the data has a curve column.
        rs: The power of ten that the response is divided by for the fit.
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
# The fit runs on the response divided by a power of ten, as in the figure.
rs = {rs!r}
methods = {methods!r}
starts = {{
    {start_txt}
}}
fig, ax = plt.subplots(figsize=({6.2 if has_group else 4.6}, 4.2))
for g, p0 in starts.items():
    s = d[d["group"] == g]
    popt, pcov = curve_fit(fpl, s["lx"], s["y"] / rs, p0=p0, method=methods[g], maxfev=20000,
                           ftol=1e-12, xtol=1e-12, gtol=1e-12)
    print(g, "bottom, top, log10 EC50, hill:", popt)
    grid = np.linspace(s["lx"].min(), s["lx"].max(), 200)
    ax.plot(10 ** grid, fpl(grid, *popt) * rs, label=g)
    ax.scatter(s["dose"], s["y"], s=14, alpha=0.6)
# A curve whose midpoint is not determined is drawn as points only, as in the figure.
for g in d["group"].unique():
    if g not in starts:
        s = d[d["group"] == g]
        ax.scatter(s["dose"], s["y"], s=14, alpha=0.6, label=f"{{g}} (not determined)")
ax.set_xscale("log")
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
