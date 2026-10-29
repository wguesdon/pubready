"""Recipe: standard curve interpolation. Mirrors the R recipe.

Rows with a concentration in --x are standards. Rows with an empty concentration
are unknowns, named by --id (or by row number), whose concentration is read off
the fitted curve. ``--test auto|linear|4pl`` picks the curve: auto fits the
four-parameter logistic of dose_response (same parameterisation, same starting
values, same fitting routine) when the standards have at least 5 distinct
positive concentrations and the fit converges, and a straight line otherwise.
The 4PL is fitted on log10 concentration, so a zero standard is left out of it.

An unknown whose signal lies outside the signal range of the curve between the
lowest and the highest standard is flagged "below range" or "above range" and
gets no concentration: nothing is extrapolated. Each standard level is read back
through the curve, and its recovery (back-calculated over nominal, in %) is the
QC of the fit.
"""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, NullLocator

from ..bundle import slugify
from ..theme import apply_pub_style, palette, plain_log_axis
from ..util import cap_first
from .dose_response import fit_curve, fpl
from .estimation_plot import nice_ticks

STD_LABEL = "Standard (mean ± SD)"
UNK_LABEL = "Unknown (interpolated)"


def predict(model, conc):
    """The fitted signal at each concentration, for either model.

    Args:
        model: The fitted model from ``fit_4pl`` or ``fit_linear``.
        conc: The concentrations.

    Returns:
        numpy.ndarray: The signals.
    """
    p = model["parameters"]
    conc = np.asarray(conc, dtype=float)
    if model["type"] == "4pl":
        return fpl(np.log10(conc), p["bottom"], p["top"], p["log_ec50"], p["hill"])
    return p["intercept"] + p["slope"] * conc


def invert(model, sig):
    """Read each signal off the curve, with its range flag, as the R engine does.

    Args:
        model: The fitted model.
        sig: The signals.

    Returns:
        tuple: The concentrations (NaN when out of range) and the flags.
    """
    sig = np.asarray(sig, dtype=float)
    ends = predict(model, model["conc_range"])
    lo, hi = ends.min(), ends.max()
    inside = (sig >= lo) & (sig <= hi)
    # A signal past the end of the lowest standard means a concentration below it.
    low_side = (sig - ends[0]) * (ends[1] - ends[0]) < 0
    flag = np.where(inside, "in range", np.where(low_side, "below range", "above range"))
    p = model["parameters"]
    conc = np.full(len(sig), np.nan)
    s = sig[inside]
    if model["type"] == "4pl":
        conc[inside] = 10 ** (p["log_ec50"]
                              - np.log10((p["top"] - p["bottom"]) / (s - p["bottom"]) - 1) / p["hill"])
    else:
        conc[inside] = (s - p["intercept"]) / p["slope"]
    return conc, flag


def fit_4pl(conc, sig):
    """Fit the 4PL to the positive standards with the dose_response routine.

    Args:
        conc: The standard concentrations.
        sig: The standard signals.

    Returns:
        dict or None: The model, or None when the fit does not converge.
    """
    keep = conc > 0
    lx, y = np.log10(conc[keep]), sig[keep]
    popt, _, st, method = fit_curve(lx, y)
    if popt is None:
        return None
    bottom, top, log_ec50, hill = (float(v) for v in popt)
    if bottom > top:
        bottom, top, hill = top, bottom, -hill
    fitted = fpl(lx, bottom, top, log_ec50, hill)
    return {
        "type": "4pl", "label": "four-parameter logistic",
        "parameters": {"bottom": bottom, "top": top, "ec50": 10 ** log_ec50,
                       "log_ec50": log_ec50, "hill": hill},
        "r_squared": float(1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)),
        "n_fit": int(len(y)), "conc_range": [float(conc[keep].min()), float(conc[keep].max())],
        "start": st, "algorithm": method,
    }


def fit_linear(conc, sig):
    """Fit a straight line to the standards by least squares.

    Args:
        conc: The standard concentrations.
        sig: The standard signals.

    Returns:
        dict: The model.
    """
    mx, my = conc.mean(), sig.mean()
    slope = float(np.sum((conc - mx) * (sig - my)) / np.sum((conc - mx) ** 2))
    intercept = float(my - slope * mx)
    resid = sig - (intercept + slope * conc)
    return {
        "type": "linear", "label": "linear",
        "parameters": {"intercept": intercept, "slope": slope},
        "r_squared": float(1 - np.sum(resid ** 2) / np.sum((sig - my) ** 2)),
        "n_fit": int(len(sig)), "conc_range": [float(conc.min()), float(conc.max())],
    }


def recovery(model, conc, sig):
    """The recovery of each positive standard level, read back through the curve.

    Args:
        model: The fitted model.
        conc: The standard concentrations.
        sig: The standard signals.

    Returns:
        dict: The levels and the summary, with the keys of the R engine.
    """
    lv = np.unique(conc[conc > 0])
    m = np.array([sig[conc == u].mean() for u in lv])
    back, _ = invert(model, m)
    rec = 100 * back / lv
    ok = np.isfinite(rec)
    levels = [{"concentration": float(c), "mean_signal": float(s),
               "back_calculated": float(b) if np.isfinite(b) else None,
               "recovery_pct": float(r) if np.isfinite(r) else None}
              for c, s, b, r in zip(lv, m, back, rec)]
    return {
        "levels": levels, "n_levels": int(len(lv)), "n_read_back": int(ok.sum()),
        "min_pct": float(rec[ok].min()) if ok.any() else None,
        "max_pct": float(rec[ok].max()) if ok.any() else None,
        "n_within_80_120": int(np.sum(ok & (rec >= 80) & (rec <= 120))),
    }


def log_breaks(lo, hi):
    """Ticks of the 1-2-5 sequence on a log axis, as the R engine computes them.

    The 1-10 sequence replaces it when the range spans more than 2.5 decades.

    Args:
        lo: The smallest concentration.
        hi: The largest concentration.

    Returns:
        list: The tick positions.
    """
    mult = [1] if np.log10(hi / lo) > 2.5 else [1, 2, 5]
    seq = sorted(m * 10.0 ** k for k in range(-6, 10) for m in mult)
    inside = [v for v in seq if lo * (1 - 1e-9) <= v <= hi * (1 + 1e-9)]
    if len(inside) < 2:
        below = [v for v in seq if v < lo]
        above = [v for v in seq if v > hi]
        inside = below[-1:] + inside + above[:1]
    return inside


def recipe_standard_curve(df, spec):
    """Fit a standard curve and interpolate the unknowns.

    Args:
        df: The tidy table: standards with a concentration, unknowns without.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if the standards are too few for
            the requested model, or if the requested fit does not converge.
    """
    x, y, id_col = spec["data"].get("x"), spec["data"].get("y"), spec["data"].get("id")
    if not x or not y:
        raise ValueError("standard_curve needs --x (the concentration, empty for an unknown) "
                         "and --y (the signal)")
    for c in [x, y] + ([id_col] if id_col else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")
    method = (spec["test"].get("method") or "auto").lower()
    if method not in ("auto", "linear", "4pl"):
        raise ValueError(f"standard_curve --test must be auto, linear or 4pl; got '{method}'")

    n0 = len(df)
    raw = df[x]
    blank = (raw.isna() | (raw.astype(str).str.strip() == "")).to_numpy()
    conc = pd.to_numeric(raw, errors="coerce").to_numpy(float)
    sig = pd.to_numeric(df[y], errors="coerce").to_numpy(float)
    ids = (df[id_col].astype(str).to_numpy() if id_col
           else np.array([f"row {i}" for i in range(1, n0 + 1)]))
    bad_sig = ~np.isfinite(sig)
    with np.errstate(invalid="ignore"):
        bad_conc = ~blank & (~np.isfinite(conc) | (conc < 0))
    keep = ~bad_sig & ~bad_conc
    steps = []
    if bad_sig.sum():
        steps.append(f"Dropped {int(bad_sig.sum())} row(s) with a missing or non-numeric signal in '{y}'.")
    if (bad_conc & ~bad_sig).sum():
        steps.append(f"Dropped {int((bad_conc & ~bad_sig).sum())} row(s) with a concentration that is "
                     "not a number of zero or more.")
    steps.append(f"Kept {int(keep.sum())} of {n0} row(s): {int((keep & ~blank).sum())} standard(s) and "
                 f"{int((keep & blank).sum())} unknown(s).")
    is_std, is_unk = keep & ~blank, keep & blank
    s_conc, s_sig = conc[is_std], sig[is_std]
    n_pos_levels = len(np.unique(s_conc[s_conc > 0]))

    model = None
    if method in ("auto", "4pl"):
        enough = (n_pos_levels >= 5 if method == "auto"
                  else n_pos_levels >= 4 and int((s_conc > 0).sum()) >= 5)
        if method == "4pl" and not enough:
            raise ValueError("standard_curve --test 4pl needs at least 4 distinct positive standard "
                             f"concentrations and 5 positive standards; found {n_pos_levels} and "
                             f"{int((s_conc > 0).sum())}")
        if enough:
            model = fit_4pl(s_conc, s_sig)
        if method == "4pl" and model is None:
            raise ValueError("the four-parameter logistic did not converge on the standards")
    if model is None:
        if len(np.unique(s_conc)) < 2 or len(s_conc) < 3:
            raise ValueError(f"standard_curve needs at least 3 standards at 2 or more concentrations; "
                             f"found {len(s_conc)} at {len(np.unique(s_conc))}")
        model = fit_linear(s_conc, s_sig)
        if model["parameters"]["slope"] == 0:
            raise ValueError("the linear standard curve is flat, so no unknown can be read off it")
    if model["type"] == "4pl" and (s_conc == 0).any():
        steps.append(f"Left {int((s_conc == 0).sum())} standard(s) at a concentration of zero out of the "
                     "four-parameter logistic fit, because it is fitted on log10 concentration.")

    conc_u, flag = invert(model, sig[is_unk])
    stats_df = pd.DataFrame({"id": ids[is_unk], "signal": sig[is_unk],
                             "concentration": conc_u, "flag": flag})
    rec = recovery(model, s_conc, s_sig)
    counts = {"n_unknowns": int(len(stats_df)),
              "n_in_range": int((stats_df["flag"] == "in range").sum()),
              "n_below": int((stats_df["flag"] == "below range").sum()),
              "n_above": int((stats_df["flag"] == "above range").sum())}
    ends = predict(model, model["conc_range"])
    test_meta = {
        "name": "Standard curve interpolation", "model": model["label"],
        "model_choice": "auto" if method == "auto" else "requested",
        "parameters": model["parameters"], "r_squared": model["r_squared"],
        "n_standards": int(len(s_conc)), "n_fitted": model["n_fit"],
        "n_levels": int(len(np.unique(s_conc))),
        "concentration_range": model["conc_range"], "signal_range": sorted(float(e) for e in ends),
        **counts, "recovery": rec,
    }
    fig = _plot(model, s_conc, s_sig, stats_df, counts, spec, x, y)
    df_used = pd.DataFrame({"id": ids[keep], "role": np.where(blank[keep], "unknown", "standard"),
                            "concentration": conc[keep], "signal": sig[keep]})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(model, test_meta, method, spec, x, y),
        "df_used": df_used, "clean_steps": steps, "label": slugify(y),
        "width": 5.2, "height": 4.8,
        "build_script": functools.partial(_script, spec, model),
    }


def _plot(model, s_conc, s_sig, unk, counts, spec, x, y):
    """Draw the standards, the curve and the interpolated unknowns.

    Args:
        model: The fitted model.
        s_conc: The standard concentrations.
        s_sig: The standard signals.
        unk: The stats table of the unknowns.
        counts: The in-range and out-of-range counts.
        spec: The figure specification.
        x: The concentration column.
        y: The signal column.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    ap = spec["appearance"]
    pal = palette(2, ap.get("palette"))
    is_4pl = model["type"] == "4pl"
    use = s_conc > 0 if is_4pl else np.ones(len(s_conc), dtype=bool)
    lv = np.unique(s_conc[use])
    means = np.array([s_sig[s_conc == u].mean() for u in lv])
    n = np.array([int((s_conc == u).sum()) for u in lv])
    sds = np.array([s_sig[s_conc == u].std(ddof=1) if k > 1 else np.nan for u, k in zip(lv, n)])
    lo, hi = model["conc_range"]
    grid = np.logspace(np.log10(lo), np.log10(hi), 200) if is_4pl else np.linspace(lo, hi, 200)
    inr = unk[unk["flag"] == "in range"]

    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    if is_4pl:
        ax.set_xscale("log")
    ax.plot(grid, predict(model, grid), color=pal[0], lw=1.6, zorder=2)
    ax.scatter(s_conc[use], s_sig[use], s=8, color=pal[0], alpha=0.4, lw=0, zorder=3)
    rep = n > 1
    if rep.any():
        ax.errorbar(lv[rep], means[rep], yerr=sds[rep], fmt="none", ecolor=pal[0], elinewidth=1.0,
                    capsize=0, zorder=3)
    ax.scatter(lv, means, s=36, color=pal[0], zorder=6)
    ys = np.concatenate([s_sig[use], means + np.nan_to_num(sds), means - np.nan_to_num(sds),
                         inr["signal"].to_numpy(float), predict(model, grid)])
    span = ys.max() - ys.min()
    ylo, yhi = ys.min() - 0.05 * span, ys.max() + 0.05 * span
    if len(inr):
        # Many unknowns would hide the standards and turn the drop lines into a
        # band, so a large set gets smaller diamonds and fainter lines.
        many = len(inr) > 50
        ax.vlines(inr["concentration"], ylo, inr["signal"], colors="0.6", lw=0.4 if many else 0.6,
                  linestyles="--", alpha=0.3 if many else 0.8, zorder=1)
        ax.scatter(inr["concentration"], inr["signal"], s=12 if many else 30, marker="D",
                   facecolors="white", edgecolors=pal[1], linewidths=0.7 if many else 1.0, zorder=5)
    ax.set_ylim(ylo, yhi)
    if is_4pl:
        # ggplot2 adds 5% of the log10 range on each side.
        pad = (hi / lo) ** 0.05
        ax.set_xticks(log_breaks(lo, hi))
        ax.set_xlim(lo / pad, hi * pad)
        plain_log_axis(ax.xaxis)
        ax.xaxis.set_minor_locator(NullLocator())
    else:
        nice_ticks(ax.xaxis)
        ax.set_xlim(lo - 0.02 * (hi - lo), hi + 0.04 * (hi - lo))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    ax.set_title(f"{cap_first(model['label'])} fit, R² = {model['r_squared']:.4f}\n"
                 f"{counts['n_in_range']} of {counts['n_unknowns']} unknowns interpolated; "
                 f"{counts['n_below']} below range, {counts['n_above']} above range",
                 fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    handles = [Line2D([], [], marker="o", ls="none", color=pal[0], markersize=7, label=STD_LABEL),
               Line2D([], [], marker="D", ls="none", markerfacecolor="white", markeredgecolor=pal[1],
                      markersize=6, label=UNK_LABEL)]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2,
              frameon=False, fontsize=10)
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()
    return fig


def _fmt(v):
    """Format a number with three significant digits, as R's %.3g does.

    Args:
        v: The number.

    Returns:
        str: The text.
    """
    return f"{v:.3g}"


def _methods(model, meta, method, spec, x, y):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        model: The fitted model.
        meta: The test metadata.
        method: The requested model (auto, linear or 4pl).
        spec: The figure specification.
        x: The concentration column.
        y: The signal column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    x_lab = spec["appearance"].get("x_label") or x
    p = model["parameters"]
    rec = meta["recovery"]
    if method == "auto":
        choice = ("The model was chosen automatically: the four-parameter logistic when the standards "
                  "have at least five distinct positive concentrations and the fit converges, and a "
                  "straight line otherwise. ")
    else:
        choice = "The model was set by the analyst. "
    if model["type"] == "4pl":
        fit_txt = ("The four-parameter logistic, signal = bottom + (top - bottom) / (1 + 10^((log10 "
                   "EC50 - log10 concentration) x Hill slope)), was fitted by nonlinear least squares "
                   f"on log10 concentration to {model['n_fit']} standards: bottom = {_fmt(p['bottom'])}, "
                   f"top = {_fmt(p['top'])}, EC50 = {_fmt(p['ec50'])}, Hill slope = {p['hill']:.3f}, "
                   f"R2 = {model['r_squared']:.4f}. ")
    else:
        fit_txt = ("A straight line, signal = intercept + slope x concentration, was fitted by least "
                   f"squares to {model['n_fit']} standards: intercept = {_fmt(p['intercept'])}, "
                   f"slope = {_fmt(p['slope'])}, R2 = {model['r_squared']:.4f}. ")
    if rec["n_read_back"] > 0:
        rec_txt = (f"Read back through the curve, the standard levels had a recovery of "
                   f"{rec['min_pct']:.1f}% to {rec['max_pct']:.1f}% of the nominal concentration, and "
                   f"{rec['n_within_80_120']} of {rec['n_levels']} levels fell within 80% to 120%. ")
    else:
        rec_txt = "No standard level could be read back through the curve. "
    lib = f"scipy {v('scipy')} (curve_fit)" if model["type"] == "4pl" else f"numpy {v('numpy')}"
    return (
        f"Concentrations ({x_lab}) were interpolated from the signal ({y_lab}) on a standard curve of "
        f"{meta['n_standards']} standards at {meta['n_levels']} concentrations. {choice}{fit_txt}"
        "Each unknown was interpolated by inverting the fitted curve. An unknown whose signal lies "
        "outside the signal range of the curve between the lowest and the highest standard "
        f"({_fmt(meta['signal_range'][0])} to {_fmt(meta['signal_range'][1])}) was flagged below range "
        "or above range and given no concentration, so no value was extrapolated. "
        f"Of {meta['n_unknowns']} unknowns, {meta['n_in_range']} were in range, {meta['n_below']} "
        f"below range and {meta['n_above']} above range. {rec_txt}"
        f"Fitted in Python {platform.python_version()} with {lib} "
        f"and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, model, in_name, fig_stub):
    """Write the standalone script that refits the curve and interpolates.

    Args:
        spec: The figure specification.
        model: The fitted model, with its starting values for a 4PL.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    x, y = spec["data"]["x"], spec["data"]["y"]
    ap = spec["appearance"]
    if model["type"] == "4pl":
        st = model["start"]
        fit = f'''pos = std[std["conc"] > 0]
p0 = [{st["bottom"]!r}, {st["top"]!r}, {st["log_ec50"]!r}, {st["hill"]!r}]
popt, _ = curve_fit(fpl, np.log10(pos["conc"]), pos["sig"], p0=p0, method="{model["algorithm"]}",
                    maxfev=20000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
b, t, le, h = popt
if b > t:
    b, t, h = t, b, -h


def f(c):
    return fpl(np.log10(c), b, t, le, h)


def finv(s):
    return 10 ** (le - np.log10((t - b) / (s - b) - 1) / h)


rng = (pos["conc"].min(), pos["conc"].max())
print("bottom, top, EC50, hill:", b, t, 10 ** le, h)'''
    else:
        fit = '''slope, intercept = np.polyfit(std["conc"], std["sig"], 1)


def f(c):
    return intercept + slope * c


def finv(s):
    return (s - intercept) / slope


rng = (std["conc"].min(), std["conc"].max())
print("intercept, slope:", intercept, slope)'''
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
d = pd.DataFrame({{"conc": pd.to_numeric(raw["{x}"], errors="coerce"),
                   "sig": pd.to_numeric(raw["{y}"], errors="coerce")}})
std = d[np.isfinite(d["conc"]) & np.isfinite(d["sig"])]
unk = d[d["conc"].isna() & np.isfinite(d["sig"])].copy()
{fit}
ends = [f(rng[0]), f(rng[1])]
inside = (unk["sig"] >= min(ends)) & (unk["sig"] <= max(ends))
unk["concentration"] = np.where(inside, finv(unk["sig"].where(inside)), np.nan)
print(unk)
fig, ax = plt.subplots(figsize=(5.2, 4.8))
grid = np.{"logspace(np.log10(rng[0]), np.log10(rng[1]), 200)" if model["type"] == "4pl" else "linspace(rng[0], rng[1], 200)"}
ax.plot(grid, f(grid)); ax.scatter(std["conc"], std["sig"])
ax.scatter(unk["concentration"], unk["sig"], marker="D", facecolors="white", edgecolors="red")
{"ax.set_xscale('log')" if model["type"] == "4pl" else ""}
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
