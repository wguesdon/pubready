"""Recipe: the shape of one numeric variable (--y), overall or per group (--x).

Mirrors core/r/R/distribution.R. The figure is a histogram with a Gaussian
kernel density (solid) and a normal curve with the sample mean and SD (dashed),
both scaled to counts, and a rug when a panel holds 500 values or fewer. Groups
are panels stacked on a shared x axis. The stats table gives the moments, the
quantiles and the Shapiro-Wilk test of each group; no comparison test is run.

The bins are the same as in R: the Freedman-Diaconis width on the pooled data
with type-7 quantiles (Sturges when the IQR is zero), anchored at
floor(min / width) * width and closed on the left. The density uses the
bandwidth of R's stats::bw.nrd0, restated here, and is evaluated directly on
512 points; scipy's gaussian_kde bandwidth is not used.
"""
import functools
import math
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy import stats

from ..bundle import slugify
from ..qc import qc_moment_stats
from ..theme import apply_pub_style, palette
from .correlation_heatmap import check_finite_variance

MAX_BINS = 200
GRID_N = 512
RUG_MAX = 500


def q7(v, p):
    """Type-7 quantile with the arithmetic of R's stats::quantile.

    numpy interpolates from the upper value when the fraction is 0.5 or more,
    which can differ from R in the last bit; this keeps the bins identical.

    Args:
        v: The values.
        p: The probability, between 0 and 1.

    Returns:
        float: The quantile.
    """
    s = np.sort(np.asarray(v, dtype=float))
    n = len(s)
    index = 1 + max(n - 1, 0) * p
    lo = math.floor(index)
    hi = math.ceil(index)
    q = s[lo - 1]
    if index > lo and s[hi - 1] != q:
        h = index - lo
        q = (1 - h) * q + h * s[hi - 1]
    return float(q)


def bw_nrd0(v):
    """Silverman's rule-of-thumb bandwidth, as R's stats::bw.nrd0 computes it.

    Args:
        v: The values, at least two.

    Returns:
        float: The bandwidth, or NaN for fewer than two values.
    """
    v = np.asarray(v, dtype=float)
    if len(v) < 2:
        return float("nan")
    hi = float(np.std(v, ddof=1))
    lo = min(hi, (q7(v, 0.75) - q7(v, 0.25)) / 1.34)
    if not lo:
        lo = hi or abs(float(v[0])) or 1.0
    return 0.9 * lo * len(v) ** (-0.2)


def fd_bins(v):
    """Freedman-Diaconis bins on the pooled values, as in the R engine.

    Args:
        v: The pooled values.

    Returns:
        dict: The width, the start, the bin count, the rule and whether the
        width was widened to keep 200 bins.
    """
    v = np.asarray(v, dtype=float)
    n = len(v)
    lo, hi = float(v.min()), float(v.max())
    iqr = q7(v, 0.75) - q7(v, 0.25)
    if iqr > 0:
        width, rule = 2 * iqr * n ** (-1.0 / 3.0), "Freedman-Diaconis"
    elif hi > lo:
        width, rule = (hi - lo) / math.ceil(math.log2(n) + 1), "Sturges"
    else:
        width, rule = 1.0, "unit"
    capped = hi > lo and (hi - lo) / width > MAX_BINS
    if capped:
        width = (hi - lo) / MAX_BINS
    start = math.floor(lo / width) * width
    return {"width": width, "start": start, "n_bins": int(math.floor((hi - start) / width)) + 1,
            "rule": rule, "capped": bool(capped)}


def bin_counts(v, bins):
    """Left-closed bin counts on the shared bins.

    Args:
        v: The values of one group.
        bins: The bins from ``fd_bins``.

    Returns:
        numpy.ndarray: The count of each bin.
    """
    idx = np.floor((np.asarray(v, dtype=float) - bins["start"]) / bins["width"]).astype(int)
    idx = np.clip(idx, 0, bins["n_bins"] - 1)
    return np.bincount(idx, minlength=bins["n_bins"])


def kde(grid, v, bw):
    """The Gaussian kernel density of ``v`` on ``grid``, summed directly.

    Args:
        grid: The evaluation points.
        v: The values.
        bw: The bandwidth.

    Returns:
        numpy.ndarray: The density at each grid point.
    """
    v = np.asarray(v, dtype=float)
    out = np.zeros(len(grid))
    for k in range(0, len(v), 20000):
        z = (grid[:, None] - v[None, k:k + 20000]) / bw
        out += np.exp(-0.5 * z * z).sum(axis=1)
    return out / (len(v) * bw * math.sqrt(2 * math.pi))


def shapiro(v):
    """Shapiro-Wilk W and P for 3 to 5000 values that are not all equal.

    Args:
        v: The values.

    Returns:
        tuple: W and P, or two NaN values outside that range.
    """
    v = np.asarray(v, dtype=float)
    if len(v) < 3 or len(v) > 5000 or np.std(v) == 0:
        return float("nan"), float("nan")
    try:
        res = stats.shapiro(v)
    except ValueError:
        return float("nan"), float("nan")
    return float(res.statistic), float(res.pvalue)


def _num_label(v):
    """Format one number as R's as.character does, with 15 significant digits.

    Args:
        v: The number.

    Returns:
        str: The label, "2" for 2.0 and "0.5" for 0.5.
    """
    return f"{v:.15g}"


def as_labels(s):
    """Turn a column into text labels that match the R engine.

    pandas reads a numeric column with a missing value as float, so ``astype(str)``
    gives "2.0" where R gives "2". Numbers are therefore formatted as R does.

    Args:
        s: The column, a pandas Series without missing values.

    Returns:
        pandas.Series: The labels.
    """
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        return s.map(_num_label)
    return s.astype(str)


def p_label(p):
    """Format a P value as the R engine does.

    Args:
        p: The P value.

    Returns:
        str: The label.
    """
    if p is None or not np.isfinite(p):
        return "P = NA"
    return "P < 0.001" if p < 0.001 else f"P = {p:.2g}"


def group_row(gl, values, bins):
    """Summarise one group for the stats table.

    Args:
        gl: The group label.
        values: The values of the group.
        bins: The shared bins.

    Returns:
        dict: One row of the stats table.
    """
    v = np.asarray(values, dtype=float)
    n = len(v)
    mom = qc_moment_stats(v)
    w, p = shapiro(v)
    return {"group": gl, "n": n, "mean": float(v.mean()),
            "sd": float(np.std(v, ddof=1)) if n >= 2 else float("nan"),
            "median": q7(v, 0.5), "q1": q7(v, 0.25), "q3": q7(v, 0.75),
            "min": float(v.min()), "max": float(v.max()),
            "skewness": mom["skewness"], "excess_kurtosis": mom["kurtosis"],
            "shapiro_w": w, "shapiro_p": p, "bin_width": bins["width"],
            "bandwidth": bw_nrd0(v) if n >= 2 and np.std(v) > 0 else float("nan")}


def strip_label(row, width):
    """The panel title of one group: its name, n and the Shapiro-Wilk P.

    Args:
        row: The stats row of the group.
        width: The character count that fits on one line.

    Returns:
        str: The title, on two lines when it is long.
    """
    name = textwrap.fill(str(row["group"]), 40, break_on_hyphens=False, break_long_words=False)
    info = (f"n = {row['n']}" if not np.isfinite(row["shapiro_p"])
            else f"n = {row['n']}, Shapiro-Wilk {p_label(row['shapiro_p'])}")
    return f"{name}\n{info}" if len(name) + len(info) + 3 > width else f"{name}: {info}"


def layout(k):
    """The panel grid: one column up to six groups, then two, then three.

    Args:
        k: The number of groups.

    Returns:
        dict: The columns, the rows, the figure size and the title width.
    """
    ncol = 1 if k <= 6 else 2 if k <= 12 else 3
    nrow = math.ceil(k / ncol)
    panel_h = 1.6 if nrow <= 3 else 1.25
    return {"ncol": ncol, "nrow": nrow, "width": 5.6 if ncol == 1 else 1.2 + 3.6 * ncol,
            "height": 1.5 + nrow * panel_h, "strip_chars": 52 if ncol == 1 else 34}


def recipe_distribution(df, spec):
    """Describe the distribution of one numeric variable, overall or per group.

    Args:
        df: The tidy data.
        spec: The figure specification.

    Returns:
        dict: The recipe result for the bundle writer.
    """
    y, x = spec["data"].get("y"), spec["data"].get("x")
    has_group = bool(x)
    if not y:
        raise ValueError("distribution needs --y (the numeric variable)")
    for c in [y] + ([x] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    steps = []
    n0 = len(df)
    yv = pd.to_numeric(df[y], errors="coerce")
    coerced = int((yv.isna() & df[y].notna()).sum())
    if coerced:
        steps.append(f"Coerced column '{y}' to numeric; {coerced} value(s) became NA.")
    grp = df[x] if has_group else pd.Series(["all"] * n0, index=df.index)
    keep = np.isfinite(yv.to_numpy(dtype=float)) & grp.notna().to_numpy()
    if (~keep).sum():
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing value"
                     f"{' or group' if has_group else ''} ({n0} -> {int(keep.sum())} rows).")
    d = pd.DataFrame({"group": as_labels(grp[keep]).to_numpy(),
                      "value": yv[keep].to_numpy(dtype=float)})
    if len(d) < 2:
        raise ValueError("distribution needs at least 2 finite values")
    if not steps:
        steps = ["No cleaning applied; input used as-is."]
    lv = sorted(d["group"].unique()) if has_group else ["all"]
    # Values near 1e300 or 1e-300 give a variance that is not finite at double
    # precision, and with it a silent NaN or a wrong curve. The run stops (F83).
    check_finite_variance(d["value"], y, "distribution")
    for g in lv:
        check_finite_variance(d.loc[d["group"] == g, "value"], y, "distribution")

    bins = fd_bins(d["value"].to_numpy())
    vals = {g: d.loc[d["group"] == g, "value"].to_numpy() for g in lv}
    stats_df = pd.DataFrame([group_row(g, vals[g], bins) for g in lv])
    lay = layout(len(lv))
    strips = [strip_label(r, lay["strip_chars"]) for _, r in stats_df.iterrows()]
    fig, w, h = _plot(vals, stats_df, strips, lv, lay, bins, has_group, spec, y)

    test_meta = {
        "name": "Shapiro-Wilk normality test per group (descriptive, no comparison)",
        "n": int(len(d)), "n_groups": len(lv),
        "bins": {"rule": bins["rule"], "capped_at_200": bins["capped"], "width": bins["width"],
                 "start": bins["start"], "n_bins": bins["n_bins"]},
        "groups": {r["group"]: {k: r[k] for k in stats_df.columns if k != "group"}
                   for _, r in stats_df.iterrows()},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, bins, has_group, spec, x, y),
        "df_used": d, "clean_steps": steps, "label": slugify(y), "width": w, "height": h,
        "qc": {"quantity": "values of each group" if has_group else "values",
               "panels": {g: vals[g] for g in lv}},
        "build_script": functools.partial(_script, spec, has_group),
    }


def _curves(v, bins):
    """The density and the normal curve of one group, scaled to counts.

    Args:
        v: The values of the group.
        bins: The shared bins.

    Returns:
        tuple: The grid, the density and the normal curve, or None for fewer
        than two values or no spread.
    """
    # A group with no spread has no density: bw_nrd0 falls back to a width from
    # the value itself, which draws a wide curve that stretches the shared axis
    # (F125).
    if len(v) < 2 or np.std(v) == 0:
        return None
    bw = bw_nrd0(v)
    grid = np.linspace(v.min() - 3 * bw, v.max() + 3 * bw, GRID_N)
    to_count = len(v) * bins["width"]
    s = float(np.std(v, ddof=1))
    normal = stats.norm.pdf(grid, v.mean(), s) * to_count if s > 0 else None
    return grid, kde(grid, v, bw) * to_count, normal


def header(fig, title, sub, h):
    """Draw the bold title and the grey subtitle at the top of a figure.

    Args:
        fig: The figure.
        title: The title, or None.
        sub: The subtitle.
        h: The figure height in inches.

    Returns:
        float: The top of the plot area as a fraction of the figure height.
    """
    top_in = 0.06
    if title:
        fig.text(0.5, 1 - (top_in + 0.12) / h, title, ha="center", va="center",
                 fontsize=13, fontweight="bold")
        top_in += 0.3
    n_lines = sub.count("\n") + 1
    fig.text(0.5, 1 - top_in / h, sub, ha="center", va="top", fontsize=9.5, color="#4d4d4d",
             linespacing=1.3)
    top_in += 0.17 * n_lines + 0.02
    return 1 - top_in / h


def _plot(vals, stats_df, strips, lv, lay, bins, has_group, spec, y):
    """Draw the histogram panels.

    Args:
        vals: The values of each group.
        stats_df: The stats table.
        strips: The panel titles.
        lv: The group levels.
        lay: The panel layout.
        bins: The shared bins.
        has_group: Whether there is a group column.
        spec: The figure specification.
        y: The value column.

    Returns:
        tuple: The figure, its width and its height.
    """
    ap = spec["appearance"]
    pal = palette(len(lv), ap.get("palette"))
    w, h = (lay["width"], lay["height"]) if has_group else (5.4, 4.2)
    nrow, ncol = (lay["nrow"], lay["ncol"]) if has_group else (1, 1)
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, h), sharex=True, squeeze=False)
    flat = [ax for row in axes for ax in row]
    edges = bins["start"] + bins["width"] * np.arange(bins["n_bins"])
    for i, g in enumerate(lv):
        ax = flat[i]
        v = vals[g]
        counts = bin_counts(v, bins)
        ax.bar(edges, counts, width=bins["width"], align="edge", color=pal[i], alpha=0.65,
               edgecolor="white", linewidth=0.4)
        top = float(counts.max()) if len(counts) else 1.0
        cur = _curves(v, bins)
        if cur is not None:
            grid, dens, normal = cur
            ax.plot(grid, dens, color="#262626", lw=1.1)
            top = max(top, float(dens.max()))
            if normal is not None:
                ax.plot(grid, normal, color="#262626", lw=1.1, ls=(0, (2.5, 2.5)))
                top = max(top, float(normal.max()))
        top = top if top > 0 else 1.0
        # The limits of R's expansion(mult = c(0.07, 0.06)); the rug fills the
        # bottom 4% of the panel, below the zero line, as geom_rug does.
        lo = -0.07 * top
        if len(v) <= RUG_MAX:
            ax.vlines(v, lo, lo + 0.04 * 1.13 * top, color=pal[i], alpha=0.6, lw=0.5)
        ax.set_ylim(lo, top * 1.06)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10], integer=True))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
        if has_group:
            ax.set_title(strips[i], fontsize=10, loc="left", linespacing=1.15)
        apply_pub_style(ax, ap.get("theme"))
        if not has_group:
            ax.set_ylabel("Count")
    xlab = ap.get("x_label") or ap.get("y_label") or y
    for j, ax in enumerate(flat):
        if j >= len(lv):
            ax.set_visible(False)
            flat[j - ncol].tick_params(axis="x", labelbottom=True)
    if not has_group:
        flat[0].set_xlabel(xlab)
    if has_group:
        sub = f"Shapiro-Wilk normality test per group; bin width {bins['width']:.3g}"
    else:
        r = stats_df.iloc[0]
        sub = f"n = {r['n']}"
        if np.isfinite(r["shapiro_w"]):
            sub += f", Shapiro-Wilk W = {r['shapiro_w']:.3f}, {p_label(r['shapiro_p'])}"
        if np.isfinite(r["skewness"]):
            sub += f"\nskewness = {r['skewness']:.2f}, excess kurtosis = {r['excess_kurtosis']:.2f}"
    top = header(fig, ap.get("title"), sub, h)
    handles = [Line2D([0], [0], color="#262626", lw=1.1, label="Kernel density"),
               Line2D([0], [0], color="#262626", lw=1.1, ls=(0, (2.5, 2.5)), label="Normal fit")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=10,
               handlelength=2.6, bbox_to_anchor=(0.5, 0.0))
    left, bottom = 0.0, 0.35 / h
    if has_group:
        # One "Count" title and one x title for the panel grid, as ggplot2
        # draws them, centred on the panels.
        left, bottom = 0.3 / w, 0.66 / h
        fig.text(0.012, (bottom + top) / 2, "Count", rotation=90, ha="left", va="center",
                 fontsize=13)
    fig.tight_layout(rect=(left, bottom, 1, top), h_pad=0.8)
    if has_group:
        x0 = min(ax.get_position().x0 for ax in flat[:len(lv)])
        x1 = max(ax.get_position().x1 for ax in flat[:len(lv)])
        fig.text((x0 + x1) / 2, 0.34 / h, xlab, ha="center", va="bottom", fontsize=13)
    return fig, w, h


def _methods(stats_df, bins, has_group, spec, x, y):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        stats_df: The stats table.
        bins: The shared bins.
        has_group: Whether there is a group column.
        spec: The figure specification.
        x: The group column.
        y: The value column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    if has_group:
        scope = (f"in each of {len(stats_df)} groups of {x} (n = "
                 f"{', '.join(str(int(n)) for n in stats_df['n'])})")
    else:
        scope = f"(n = {int(stats_df['n'].iloc[0])})"
    bin_txt = {
        "Freedman-Diaconis": "the Freedman-Diaconis rule (2 x IQR x n^(-1/3)) on the pooled values",
        "Sturges": "the Sturges rule, because the interquartile range of the pooled values is zero",
    }.get(bins["rule"], "a unit width, because every value is equal")
    cap_txt = (" The rule asked for more than 200 bins, so the width was widened to give 200 bins."
               if bins["capped"] else "")
    sw = stats_df[np.isfinite(stats_df["shapiro_p"].astype(float))]
    if len(sw):
        parts = "; ".join(f"{(r['group'] + ': ') if has_group else ''}W = {r['shapiro_w']:.3f}, "
                          f"{p_label(r['shapiro_p'])}" for _, r in sw.iterrows())
        sw_txt = f" Normality was assessed with the Shapiro-Wilk test ({parts})."
    else:
        sw_txt = (" The Shapiro-Wilk test was not run, because it needs 3 to 5000 values that are "
                  "not all equal.")
    cmp_txt = (" No test compared the groups; use two_group_compare or multi_group_compare to test "
               "a difference in location." if has_group else "")
    flat = stats_df.loc[stats_df["sd"].astype(float) == 0, "group"].tolist()
    if flat and has_group:
        flat_txt = (f" The group(s) {', '.join(str(g) for g in flat)} have no spread, so no density "
                    "curve or normal fit is drawn for them.")
    elif flat:
        flat_txt = " The values have no spread, so no density curve or normal fit is drawn."
    else:
        flat_txt = ""
    return (
        f"The distribution of {y_lab} was described {scope}."
        f" The histogram uses a bin width of {bins['width']:.4g} from {bin_txt}, with bins anchored "
        "at a multiple of the width." + cap_txt +
        " The solid curve is a Gaussian kernel density with the bandwidth of Silverman's rule of "
        "thumb (stats::bw.nrd0), and the dashed curve is a normal distribution with the sample mean "
        "and SD; both are scaled to counts." + flat_txt + sw_txt +
        " Skewness and excess kurtosis are the moment estimates on standardised values (sum(z^3) / n "
        "and sum(z^4) / n - 3, with z scaled by the SD)." + cmp_txt +
        f" Computed in Python {platform.python_version()} with scipy {v('scipy')} and drawn with "
        f"matplotlib {v('matplotlib')}."
    )


def _script(spec, has_group, in_name, fig_stub):
    """Write the standalone reproduction script.

    Args:
        spec: The figure specification.
        has_group: Whether there is a group column.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script.
    """
    y, x = spec["data"]["y"], spec["data"].get("x")
    ap = spec["appearance"]
    grp = f'raw["{x}"].astype(str)' if has_group else '"all"'
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("{in_name}")
d = pd.DataFrame({{"group": {grp}, "value": pd.to_numeric(raw["{y}"], errors="coerce")}})
d = d[np.isfinite(d["value"]) & d["group"].notna()]
v = d["value"].to_numpy(); n = len(v)
iqr = np.quantile(v, 0.75) - np.quantile(v, 0.25)
w = 2 * iqr * n ** (-1 / 3) if iqr > 0 else ((v.max() - v.min()) / math.ceil(math.log2(n) + 1) if v.max() > v.min() else 1.0)
if v.max() > v.min() and (v.max() - v.min()) / w > 200:
    w = (v.max() - v.min()) / 200
start = math.floor(v.min() / w) * w; nb = int((v.max() - start) // w) + 1
groups = sorted(d["group"].unique())
fig, axes = plt.subplots(len(groups), 1, figsize=(5.6, 1.5 + 1.6 * len(groups)), sharex=True, squeeze=False)
for ax, g in zip(axes[:, 0], groups):
    vg = d.loc[d["group"] == g, "value"].to_numpy()
    idx = np.clip(np.floor((vg - start) / w).astype(int), 0, nb - 1)
    ax.bar(start + w * np.arange(nb), np.bincount(idx, minlength=nb), width=w, align="edge", alpha=0.65, edgecolor="white")
    if len(vg) >= 2 and vg.std() > 0:
        sd = vg.std(ddof=1); iq = np.quantile(vg, 0.75) - np.quantile(vg, 0.25)
        bw = 0.9 * (min(sd, iq / 1.34) or sd or 1.0) * len(vg) ** -0.2
        grid = np.linspace(vg.min() - 3 * bw, vg.max() + 3 * bw, 512)
        dens = np.exp(-0.5 * ((grid[:, None] - vg[None, :]) / bw) ** 2).sum(axis=1) / (len(vg) * bw * math.sqrt(2 * math.pi))
        ax.plot(grid, dens * len(vg) * w, color="#262626")
        ax.plot(grid, stats.norm.pdf(grid, vg.mean(), sd) * len(vg) * w, color="#262626", ls="--")
    if 3 <= len(vg) <= 5000:
        print(g, stats.shapiro(vg))
    ax.set_title(str(g) if {has_group} else "", loc="left", fontsize=10); ax.set_ylabel("Count")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1, 0].set_xlabel("{ap.get('x_label') or ap.get('y_label') or y}")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
