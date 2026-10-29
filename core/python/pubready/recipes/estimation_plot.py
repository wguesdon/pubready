"""Recipe: estimation plot in the style of DABEST. Mirrors the R recipe.

Each group is compared with one reference group, and the figure shows the mean
difference with its bootstrap confidence interval next to the raw data. Two
groups give the Gardner-Altman layout (a difference axis on the right whose zero
sits at the reference mean); three or more give the Cumming layout (raw data on
top, the differences against the reference below). ``--paired`` draws a
slopegraph of each subject and resamples whole pairs.

The bootstrap gives identical resamples in the R and the Python engine, because
neither engine uses its own random number generator. Both run the Park-Miller
minimal standard generator, which is exact in 64-bit integers and in double
precision::

    state = (48271 * state) % (2**31 - 1);  u = state / (2**31 - 1)
    index = floor(u * n)          # 0-based here, +1 in R
    seed 20261029, restarted for each comparison

Each resample draws the n_ref reference indices first and then the n_group
indices of the compared group (paired: n pair indices). The 95% interval is the
bias-corrected and accelerated (BCa) interval: z0 is the normal quantile of the
proportion of bootstrap differences strictly below the estimate, the
acceleration comes from the jackknife of the mean difference over every
observation of both groups (paired: every pair), and the ends are type-7
quantiles (numpy's default ``linear`` method) of the 5000 bootstrap differences.
"""
import functools

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import special, stats

from ..bundle import slugify
from ..clean import clean_xy
from ..theme import apply_pub_style, palette
from ..util import cap_first

SEED = 20261029
MODULUS = 2147483647
MULTIPLIER = 48271
N_BOOT = 5000
BLOCK = 4096


def block_multipliers(b):
    """The multipliers 48271**j mod m for j = 1..b, as in the R engine.

    Args:
        b: The block length.

    Returns:
        numpy.ndarray: The multipliers as int64.
    """
    out = np.empty(b, dtype=np.int64)
    v = MULTIPLIER
    for j in range(b):
        out[j] = v
        v = (v * MULTIPLIER) % MODULUS
    return out


def park_miller_uniforms(n, state, mult):
    """Draw n uniforms of the Park-Miller generator that continue from state.

    A block of the generator is ``state * mult % m``: every product is below
    2**62, so int64 arithmetic is exact and gives the states of the R engine.

    Args:
        n: The number of uniforms.
        state: The current state of the generator.
        mult: The block multipliers from ``block_multipliers``.

    Returns:
        tuple: The uniforms as float64 and the last state.
    """
    out = np.empty(n, dtype=np.int64)
    pos, b = 0, len(mult)
    state = np.int64(state)
    while pos < n:
        blk = (state * mult) % MODULUS
        k = min(b, n - pos)
        out[pos:pos + k] = blk[:k]
        state = blk[k - 1]
        pos += k
    return out / float(MODULUS), int(state)


def boot_diffs(r, g, paired, n_boot=N_BOOT):
    """The bootstrap distribution of the mean difference (group minus reference).

    Args:
        r: The reference values (paired: matched by position with ``g``).
        g: The values of the compared group.
        paired: Whether the values are pairs.
        n_boot: The number of resamples.

    Returns:
        numpy.ndarray: The bootstrap differences in resample order.
    """
    mult = block_multipliers(BLOCK)
    state = SEED
    n_r, n_g = len(r), len(g)
    per = n_r if paired else n_r + n_g
    chunk = max(1, int(np.floor(2e6 / per)))
    out = np.empty(n_boot)
    done = 0
    while done < n_boot:
        m = min(chunk, n_boot - done)
        u, state = park_miller_uniforms(m * per, state, mult)
        u = u.reshape(m, per)
        if paired:
            idx = np.floor(u * n_r).astype(np.int64)
            diffs = (g[idx] - r[idx]).sum(axis=1) / n_r
        else:
            ir = np.floor(u[:, :n_r] * n_r).astype(np.int64)
            ig = np.floor(u[:, n_r:] * n_g).astype(np.int64)
            diffs = g[ig].sum(axis=1) / n_g - r[ir].sum(axis=1) / n_r
        out[done:done + m] = diffs
        done += m
    return out


def jackknife(r, g, paired):
    """Jackknife values of the mean difference.

    Args:
        r: The reference values.
        g: The values of the compared group.
        paired: Whether to leave one pair out rather than one observation.

    Returns:
        numpy.ndarray: The leave-one-out mean differences.
    """
    if paired:
        d = g - r
        return (d.sum() - d) / (len(d) - 1)
    mr, mg = r.sum() / len(r), g.sum() / len(g)
    return np.concatenate([mg - (r.sum() - r) / (len(r) - 1),
                           (g.sum() - g) / (len(g) - 1) - mr])


def bca_interval(boot, est, jack, scale):
    """The BCa 95% interval, as the R engine computes it.

    A bootstrap difference counts as below the estimate only when it is below by
    more than 1e-11 of the data scale, so that rounding in the last bit of a sum
    cannot decide a tie that is exact in real arithmetic.

    Args:
        boot: The bootstrap differences.
        est: The observed mean difference.
        jack: The jackknife values.
        scale: The largest absolute data value.

    Returns:
        tuple: The lower and upper ends.
    """
    n_boot = len(boot)
    if np.std(boot, ddof=1) == 0:
        return est, est
    below = np.sum(boot < est - 1e-11 * scale) / n_boot
    below = min(max(below, 0.5 / n_boot), 1 - 0.5 / n_boot)
    z0 = stats.norm.ppf(below)
    dj = jack.mean() - jack
    den = 6 * np.sum(dj ** 2) ** 1.5
    acc = np.sum(dj ** 3) / den if den > 0 else 0.0
    zq = stats.norm.ppf([0.025, 0.975])
    alpha = stats.norm.cdf(z0 + (z0 + zq) / (1 - acc * (z0 + zq)))
    lo, hi = np.quantile(boot, alpha)
    return float(lo), float(hi)


def hedges_g(r, g, paired):
    """Hedges' g as DABEST computes it.

    Cohen's d over the pooled SD (paired: over the average SD of the two
    conditions), times the exact correction J for n1 + n2 - 2 degrees of freedom.

    Args:
        r: The reference values.
        g: The values of the compared group.
        paired: Whether the values are pairs.

    Returns:
        float: The effect size, positive when the group is above the reference.
    """
    n1, n2 = len(r), len(g)
    v1, v2 = np.var(r, ddof=1), np.var(g, ddof=1)
    s = np.sqrt((v1 + v2) / 2) if paired else np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    dfj = n1 + n2 - 2
    j = np.exp(special.gammaln(dfj / 2) - special.gammaln((dfj - 1) / 2)) / np.sqrt(dfj / 2)
    return float((g.mean() - r.mean()) / s * j)


def t_test(r, g, paired):
    """The Welch two-sample t-test, or the paired t-test, of g against r.

    Args:
        r: The reference values.
        g: The values of the compared group.
        paired: Whether the values are pairs.

    Returns:
        tuple: The t statistic, the degrees of freedom and the two-sided P.
    """
    if paired:
        res = stats.ttest_rel(g, r)
        return float(res.statistic), float(len(r) - 1), float(res.pvalue)
    res = stats.ttest_ind(g, r, equal_var=False)
    a, b = np.var(g, ddof=1) / len(g), np.var(r, ddof=1) / len(r)
    df = (a + b) ** 2 / (a ** 2 / (len(g) - 1) + b ** 2 / (len(r) - 1))
    return float(res.statistic), float(df), float(res.pvalue)


def compare(r, g, paired, ref_name, grp_name):
    """Compare one group with the reference.

    Args:
        r: The reference values.
        g: The values of the compared group.
        paired: Whether the values are pairs.
        ref_name: The reference level.
        grp_name: The compared level.

    Returns:
        tuple: The stats row and the bootstrap differences.
    """
    est = g.sum() / len(g) - r.sum() / len(r)
    boot = boot_diffs(r, g, paired)
    lo, hi = bca_interval(boot, est, jackknife(r, g, paired), float(np.max(np.abs(np.concatenate([r, g])))))
    t, df, p = t_test(r, g, paired)
    row = {
        "reference": ref_name, "group": grp_name, "n_ref": int(len(r)), "n_group": int(len(g)),
        "mean_ref": float(r.mean()), "mean_group": float(g.mean()), "mean_diff": float(est),
        "ci_low": lo, "ci_high": hi, "hedges_g": hedges_g(r, g, paired),
        "test": "paired t-test" if paired else "Welch two-sample t-test",
        "statistic": t, "df": df, "p_value": p,
    }
    return row, boot


def swarm_offsets(v, vmin, bin_h, step=0.07, half=0.3):
    """Swarm offsets, as the R engine computes them.

    The values are cut into bins of height bin_h, and the points of one bin
    spread out from the centre (0, +1, -1, +2, ...) in value order.

    Args:
        v: The values of one group.
        vmin: The smallest value of all groups.
        bin_h: The bin height.
        step: The largest offset between neighbours.
        half: The largest offset from the centre.

    Returns:
        numpy.ndarray: The x offsets.
    """
    v = np.asarray(v, dtype=float)
    off = np.zeros(len(v))
    if len(v) == 0:
        return off
    bins = np.floor((v - vmin) / bin_h)
    for b in pd.unique(bins):
        sel = np.flatnonzero(bins == b)
        sel = sel[np.argsort(v[sel], kind="stable")]
        k = len(sel)
        if k < 2:
            continue
        i = np.arange(1, k + 1)
        slots = np.ceil((i - 1) / 2) * np.where(i % 2 == 0, 1.0, -1.0)
        slots = slots - slots.mean()
        off[sel] = slots * min(step, 2 * half / (k - 1))
    return off


def density(boot):
    """Gaussian kernel density of the bootstrap differences, as in the R engine.

    Args:
        boot: The bootstrap differences.

    Returns:
        tuple: The 128 grid values and the density scaled to a maximum of 1.
    """
    grid = np.linspace(boot.min(), boot.max(), 128)
    q1, q3 = np.quantile(boot, [0.25, 0.75])
    s = min(np.std(boot, ddof=1), (q3 - q1) / 1.34)
    if not np.isfinite(s) or s <= 0:
        s = np.std(boot, ddof=1)
    bw = 0.9 * s * len(boot) ** (-0.2)
    dens = stats.norm.pdf(np.subtract.outer(boot, grid) / bw).mean(axis=0)
    return grid, dens / dens.max()


def nice_ticks(axis):
    """Put about five ticks on 1-2-5 steps, as ggplot2 does, on an axis.

    Args:
        axis: The axis, for example ``ax.yaxis``.
    """
    axis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))


def shift(v, offset):
    """Add a constant; the transform of the secondary difference axis.

    Args:
        v: The values.
        offset: The constant.

    Returns:
        The shifted values.
    """
    return np.asarray(v) + offset


def ga_width(labels):
    """The width of the Gardner-Altman figure, grown for long column labels.

    Each of the three columns gets 1.1 in, or 0.1 in per character of the longest
    label line when that is wider, plus 1.6 in for the two axes. The R engine
    uses the same rule.

    Args:
        labels: The three wrapped column labels.

    Returns:
        float: The width in inches, at least 4.9.
    """
    chars = max(len(line) for lab in labels for line in str(lab).split("\n"))
    return max(4.9, 3 * max(1.1, 0.1 * chars) + 1.6)

def wrap(label, width=14):
    """Wrap a category label into lines shorter than width, as the R engine does.

    The rule is that of R's scales::label_wrap (each line shorter than the
    width, words kept whole), and a line may also break after an underscore, so
    that a name such as ``stimulated_plus_inhibitor`` wraps too.

    Args:
        label: The label.
        width: The width; every line is shorter than it unless one word is longer.

    Returns:
        str: The wrapped label.
    """
    out, line = [], ""
    for w in str(label).split(" "):
        if not w:
            continue
        parts = [p for p in w.replace("_", "_\x01").split("\x01") if p]
        for j, p in enumerate(parts):
            sep = " " if j == 0 and line else ""
            if line and len(line) + len(sep) + len(p) >= width:
                out.append(line)
                line = p
            else:
                line += sep + p
    return "\n".join(out + [line])


def point_area(d):
    """The marker area of the raw data, smaller for a large group.

    It matches the R sizes 2, 1.6 and 1.2 mm.

    Args:
        d: The cleaned data with a group column.

    Returns:
        float: The marker area in points squared.
    """
    n = int(d["group"].value_counts().max())
    return 8.0 if n > 100 else (14.0 if n > 40 else 22.0)


def fmt_p(p):
    """Format a P value for the figure.

    Args:
        p: The P value.

    Returns:
        str: The label.
    """
    return "P < 0.001" if p < 0.001 else f"P = {p:.2g}"


def id_text(s):
    """Subject ids as text, with integer-valued numbers written as R writes them.

    Args:
        s: The id column.

    Returns:
        pandas.Series: The ids as strings, missing values kept as NaN.
    """
    if pd.api.types.is_numeric_dtype(s):
        out = s.map(_number_text)
        return out.where(s.notna())
    return s.astype(str).where(s.notna())


def _number_text(v):
    """Write one number as R's as.character does for the usual ids.

    Args:
        v: The number.

    Returns:
        str: The text.
    """
    if pd.isna(v):
        return ""
    return str(int(v)) if float(v).is_integer() else repr(float(v))


def recipe_estimation_plot(df, spec):
    """Draw an estimation plot: each group against the reference with a BCa CI.

    Args:
        df: The tidy table, one row per value.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a group is too small, or if the
            reference level is not in the data.
    """
    x, y = spec["data"].get("x"), spec["data"].get("y")
    if not x or not y:
        raise ValueError("estimation_plot needs --x (the group column) and --y (the value column)")
    paired = bool(spec["test"].get("paired"))
    id_col = spec["data"].get("id")
    if paired and not id_col:
        raise ValueError("estimation_plot --paired needs --id (the subject column)")
    method = (spec["test"].get("method") or "auto").lower()
    if method not in ("auto", "welch_t", "t"):
        raise ValueError("estimation_plot reports a Welch t-test (a paired t-test with --paired); "
                         "use --test auto")
    if paired and id_col not in df.columns:
        raise ValueError(f"column '{id_col}' not found in data")

    d0, steps = clean_xy(df, x, y)
    steps = [s for s in steps if s != "No cleaning applied; input used as-is."]
    d = pd.DataFrame({"group": d0[x].astype(str), "y": d0[y].astype(float),
                      "id": id_text(d0[id_col]) if paired else np.nan})
    if paired:
        miss = d["id"].isna() | (d["id"] == "")
        if miss.any():
            steps.append(f"Dropped {int(miss.sum())} row(s) with a missing '{id_col}'.")
            d = d[~miss].reset_index(drop=True)
        if d.duplicated(["id", "group"]).any():
            raise ValueError(f"estimation_plot --paired needs one value per subject and group; "
                             f"'{id_col}' repeats a subject within a group")

    lv_sorted = sorted(d["group"].unique())
    ref = spec["data"].get("ref") or lv_sorted[0]
    if ref not in lv_sorted:
        raise ValueError(f"--ref '{ref}' is not a level of '{x}' (levels: {', '.join(lv_sorted)})")
    lv = [ref] + [g for g in lv_sorted if g != ref]
    if len(lv) < 2:
        raise ValueError(f"estimation_plot needs at least 2 groups in '{x}'")

    rows, boots = [], []
    for g in lv[1:]:
        if paired:
            a = d.loc[d["group"] == ref, ["id", "y"]]
            b = d.loc[d["group"] == g, ["id", "y"]]
            m = a.merge(b, on="id", suffixes=("_r", "_g")).sort_values("id", kind="stable")
            lost = len(set(a["id"]) | set(b["id"])) - len(m)
            if lost > 0:
                steps.append(f"Comparison {g} against {ref}: left out {lost} subject(s) without a "
                             "value in both groups.")
            if len(m) < 3:
                raise ValueError(f"estimation_plot --paired needs at least 3 subjects with a value in "
                                 f"both '{ref}' and '{g}'")
            r_v, g_v = m["y_r"].to_numpy(float), m["y_g"].to_numpy(float)
        else:
            r_v = d.loc[d["group"] == ref, "y"].to_numpy(float)
            g_v = d.loc[d["group"] == g, "y"].to_numpy(float)
            if len(r_v) < 2 or len(g_v) < 2:
                raise ValueError(f"estimation_plot needs at least 2 values in each group; '{ref}' has "
                                 f"{len(r_v)} and '{g}' has {len(g_v)}")
        row, boot = compare(r_v, g_v, paired, ref, g)
        rows.append(row)
        boots.append(boot)
    if not steps:
        steps = ["No cleaning applied; input used as-is."]
    stats_df = pd.DataFrame(rows)

    ap = spec["appearance"]
    pal = dict(zip(lv, palette(len(lv), ap.get("palette"))))
    y_lab = ap.get("y_label") or y
    if len(lv) == 2:
        fig, w, h = _gardner_altman(d, lv, boots, stats_df, pal, paired, ap, y_lab)
    else:
        fig, w, h = _cumming(d, lv, boots, stats_df, pal, paired, ap, y_lab)

    test_name = ("paired mean difference with a BCa bootstrap CI" if paired
                 else "mean difference with a BCa bootstrap CI")
    test_meta = {
        "name": test_name, "reference": ref, "paired": paired,
        "bootstrap": {"resamples": N_BOOT,
                      "generator": "Park-Miller minimal standard (48271, 2^31 - 1)",
                      "seed": SEED, "interval": "bias-corrected and accelerated (BCa), 95%"},
        "test": stats_df["test"].iloc[0],
        "comparisons": {r["group"]: {k: r[k] for k in ("mean_diff", "ci_low", "ci_high", "hedges_g",
                                                        "statistic", "df", "p_value")}
                        for r in rows},
    }
    used = d[["group", "y"] + (["id"] if paired else [])].copy()
    used.columns = [x, y] + ([id_col] if paired else [])
    if paired:
        panels = {f"{g} - {ref}": _paired_diffs(d, ref, g) for g in lv[1:]}
        qc = {"quantity": "paired differences", "panels": panels}
    else:
        qc = {"quantity": "group values",
              "panels": {g: d.loc[d["group"] == g, "y"].to_numpy() for g in lv}}
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(stats_df, ref, paired, spec, y),
        "df_used": used, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h, "qc": qc,
        "build_script": functools.partial(_script, spec, lv, paired),
    }


def _paired_diffs(d, ref, g):
    """Paired differences of one group against the reference.

    Args:
        d: The cleaned data with group, y and id columns.
        ref: The reference level.
        g: The compared level.

    Returns:
        numpy.ndarray: The differences, group minus reference.
    """
    m = d.loc[d["group"] == ref, ["id", "y"]].merge(d.loc[d["group"] == g, ["id", "y"]], on="id")
    return (m["y_y"] - m["y_x"]).to_numpy()


def _axis_labels(lv, d):
    """Wrapped group labels with the group size.

    Args:
        lv: The levels in display order.
        d: The cleaned data.

    Returns:
        list: The labels.
    """
    return [f"{wrap(g)}\n(n = {int((d['group'] == g).sum())})" for g in lv]


def _draw_raw(ax, d, lv, pos, pal, paired, size):
    """Draw the raw values: a swarm, or the subject lines of a paired design.

    Args:
        ax: The axes.
        d: The cleaned data.
        lv: The levels in display order.
        pos: The x position of each level.
        pal: The colour of each level.
        paired: Whether to draw the subject lines.
        size: The marker area.

    Returns:
        list: The x of each group's mean +/- SD line, just right of its swarm.
    """
    vmin = d["y"].min()
    bin_h = (d["y"].max() - vmin) / 40
    if not np.isfinite(bin_h) or bin_h <= 0:
        bin_h = 1.0
    gap_x = []
    if paired:
        for i in range(len(lv) - 1):
            a = d.loc[d["group"] == lv[i], ["id", "y"]]
            b = d.loc[d["group"] == lv[i + 1], ["id", "y"]]
            m = a.merge(b, on="id")
            for ya, yb in zip(m["y_x"], m["y_y"]):
                ax.plot([pos[i], pos[i + 1]], [ya, yb], color="0.55", lw=0.6, alpha=0.7, zorder=1)
    for i, g in enumerate(lv):
        v = d.loc[d["group"] == g, "y"].to_numpy()
        off = np.zeros(len(v)) if paired else swarm_offsets(v, vmin, bin_h)
        ax.scatter(pos[i] + off, v, s=size, color=pal[g], edgecolors="black", linewidths=0.25,
                   alpha=0.85, zorder=3)
        gap_x.append(pos[i] + (np.abs(off).max() if len(off) else 0.0) + 0.15)
    return gap_x


def _draw_gapped(ax, d, lv, xpos, gap):
    """Draw each group's mean +/- SD as a line with a gap at the mean.

    Args:
        ax: The axes.
        d: The cleaned data.
        lv: The levels in display order.
        xpos: The x position of each line.
        gap: The half gap around the mean.
    """
    for g, xp in zip(lv, xpos):
        v = d.loc[d["group"] == g, "y"].to_numpy()
        m = v.mean()
        s = np.std(v, ddof=1) if len(v) > 1 else 0.0
        if m - s < m - gap:
            ax.plot([xp, xp], [m - s, m - gap], color="black", lw=1.6, solid_capstyle="butt", zorder=4)
        if m + gap < m + s:
            ax.plot([xp, xp], [m + gap, m + s], color="black", lw=1.6, solid_capstyle="butt", zorder=4)


def _draw_difference(ax, boot, pos, offset, est, lo, hi, colour):
    """Draw one bootstrap half violin with the estimate and its CI.

    Args:
        ax: The axes.
        boot: The bootstrap differences.
        pos: The x position.
        offset: The y shift (the reference mean on a Gardner-Altman plot).
        est: The mean difference.
        lo: The lower end of the CI.
        hi: The upper end of the CI.
        colour: The fill colour.
    """
    grid, w = density(boot)
    xs = np.concatenate([[pos], pos + 0.38 * w, [pos]])
    ys = offset + np.concatenate([[grid[0]], grid, [grid[-1]]])
    ax.fill(xs, ys, color=colour, alpha=0.45, lw=0, zorder=2)
    ax.plot([pos, pos], [offset + lo, offset + hi], color="black", lw=1.8, zorder=5)
    ax.scatter([pos], [offset + est], s=40, color="black", zorder=6)


def _gardner_altman(d, lv, boots, rows, pal, paired, ap, y_lab):
    """Two groups: the Gardner-Altman layout in one axes.

    Args:
        d: The cleaned data.
        lv: The two levels, reference first.
        boots: The bootstrap differences of the comparison.
        rows: The stats table.
        pal: The colour of each level.
        paired: Whether the design is paired.
        ap: The appearance spec.
        y_lab: The y-axis title.

    Returns:
        tuple: The figure, its width and its height.
    """
    labels = _axis_labels(lv, d) + [wrap(f"{lv[1]} minus {lv[0]}")]
    w, h = ga_width(labels), 4.4
    fig, ax = plt.subplots(figsize=(w, h))
    r = rows.iloc[0]
    m_ref, m_grp = r["mean_ref"], r["mean_group"]
    yr = np.ptp(np.concatenate([d["y"].to_numpy(), m_ref + boots[0]]))
    for yv in (m_ref, m_grp):
        ax.plot([2.55, 3.5], [yv, yv], color="0.6", lw=0.6, ls="--", zorder=1)
    gx = _draw_raw(ax, d, lv, [1, 2], pal, paired, point_area(d))
    if paired:
        gx = [1 - 0.16, 2 + 0.16]
    _draw_gapped(ax, d, lv, gx, yr * 0.012)
    _draw_difference(ax, boots[0], 3, m_ref, r["mean_diff"], r["ci_low"], r["ci_high"], pal[lv[1]])
    ax.set_xlim(0.55, 3.55)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(labels)
    ax.set_ylabel(y_lab)
    if ap.get("x_label"):
        ax.set_xlabel(ap["x_label"])
    ax.set_title(f"Mean difference {r['mean_diff']:.3g} (95% CI {r['ci_low']:.3g} to "
                 f"{r['ci_high']:.3g})\n{r['test']}, {fmt_p(r['p_value'])}",
                 fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    nice_ticks(ax.yaxis)
    sec = ax.secondary_yaxis("right", functions=(functools.partial(shift, offset=-m_ref),
                                                 functools.partial(shift, offset=m_ref)))
    sec.set_ylabel("Mean difference", fontsize=13, color="black", rotation=270, labelpad=16)
    sec.tick_params(colors="black", labelsize=11)
    sec.spines["right"].set_linewidth(0.8)
    fig.tight_layout()
    return fig, w, h


def _cumming(d, lv, boots, rows, pal, paired, ap, y_lab):
    """Three or more groups: the Cumming layout, raw data above the differences.

    Args:
        d: The cleaned data.
        lv: The levels, reference first.
        boots: The bootstrap differences of each comparison.
        rows: The stats table.
        pal: The colour of each level.
        paired: Whether the design is paired.
        ap: The appearance spec.
        y_lab: The y-axis title of the raw data.

    Returns:
        tuple: The figure, its width and its height.
    """
    k = len(lv)
    w, h = max(5.2, 1.15 * k + 1.6), 6.4
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(w, h),
                                      gridspec_kw={"height_ratios": [1.25, 1]})
    pos = np.arange(1, k + 1, dtype=float)
    xlim = (0.5, k + 0.6)
    gap = np.ptp(d["y"].to_numpy()) * 0.012
    gx = _draw_raw(top, d, lv, pos, pal, paired, point_area(d))
    if not paired:
        # In a paired design the subject lines leave each group on both sides, so
        # the mean +/- SD lines would cross them; the slopegraph shows the spread.
        _draw_gapped(top, d, lv, gx, gap)
    top.set_xlim(*xlim)
    top.set_xticks(pos)
    top.set_xticklabels(_axis_labels(lv, d))
    top.set_ylabel(y_lab)
    top.set_title(f"Mean difference against {lv[0]} with its 95% BCa bootstrap CI",
                  fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    apply_pub_style(top, ap.get("theme"))
    nice_ticks(top.yaxis)

    bottom.axhline(0, color="0.5", lw=0.8, ls="--", zorder=1)
    span = [0.0]
    for i, (boot, (_, r)) in enumerate(zip(boots, rows.iterrows())):
        _draw_difference(bottom, boot, pos[i + 1], 0.0, r["mean_diff"], r["ci_low"], r["ci_high"],
                         pal[lv[i + 1]])
        span += [boot.min(), boot.max(), r["ci_low"], r["ci_high"]]
    lo, hi = min(span), max(span)
    for i, (_, r) in enumerate(rows.iterrows()):
        bottom.text(pos[i + 1], hi + (hi - lo) * 0.1, fmt_p(r["p_value"]), ha="center",
                    va="center", fontsize=9, color="0.3")
    bottom.set_ylim(lo - (hi - lo) * 0.06, hi + (hi - lo) * 0.2)
    bottom.set_xlim(*xlim)
    bottom.set_xticks(pos)
    bottom.set_xticklabels([])
    bottom.tick_params(axis="x", length=0)
    bottom.set_ylabel("Mean difference\n" + wrap(f"vs {lv[0]}", 20))
    if ap.get("x_label"):
        bottom.set_xlabel(ap["x_label"])
    apply_pub_style(bottom, ap.get("theme"))
    nice_ticks(bottom.yaxis)
    fig.tight_layout()
    fig.align_ylabels([top, bottom])
    return fig, w, h


def _p_text(test, k, ref):
    """Name the test of the P values and their adjustment.

    Each group is tested against the reference on its own, so with two or more
    groups the sentence says that the P values carry no adjustment.

    Args:
        test: The name of the t-test.
        k: The count of groups compared with the reference.
        ref: The reference group.

    Returns:
        str: The sentence, with a trailing space.
    """
    if k > 1:
        return (f"A {test} gives each P value, with no adjustment for the {k} comparisons "
                f"against {ref}. ")
    return f"A {test} gives the P value. "


def _methods(rows, ref, paired, spec, y):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        rows: The stats table.
        ref: The reference level.
        paired: Whether the design is paired.
        spec: The figure specification.
        y: The value column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or y
    each = "; ".join(
        f"{r['group']}: mean difference {r['mean_diff']:.3g} (95% CI {r['ci_low']:.3g} to "
        f"{r['ci_high']:.3g}), Hedges' g = {r['hedges_g']:.2f}, t({r['df']:.3g}) = "
        f"{r['statistic']:.2f}, P = {r['p_value']:.2g}" for _, r in rows.iterrows())
    draw = ("each resample drawing the subject indices of the pairs as floor(u n). " if paired else
            "each resample drawing the reference indices first and then the group indices, as "
            "floor(u n). ")
    return (
        f"{cap_first(y_lab)} was analysed with estimation statistics, with {ref} as the reference "
        f"group (n = {int(rows['n_ref'].iloc[0])}). For each group the "
        f"{'paired ' if paired else ''}mean difference against {ref} (group minus reference) is "
        "reported with a 95% bias-corrected and accelerated (BCa) bootstrap confidence interval "
        f"from {N_BOOT} resamples. The resamples were drawn with the Park-Miller minimal standard "
        "generator (multiplier 48271, modulus 2^31 - 1, seed 20261029, restarted for each "
        f"comparison), {draw}The bias correction is the normal quantile of the proportion of "
        "bootstrap differences strictly below the estimate, the acceleration comes from the "
        "jackknife of the mean difference over every "
        f"{'pair' if paired else 'observation of both groups'}, and the interval ends are type-7 "
        "quantiles. Hedges' g is Cohen's d over the "
        f"{'average standard deviation of the two conditions' if paired else 'pooled standard deviation'}"
        " with the exact small-sample correction. "
        f"{_p_text('paired t-test' if paired else 'Welch two-sample t-test', len(rows), ref)}"
        f"{each}. "
        f"Computed in Python {platform.python_version()} with numpy {v('numpy')} and scipy "
        f"{v('scipy')} and drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, lv, paired, in_name, fig_stub):
    """Write the standalone script that recomputes the intervals and redraws.

    Args:
        spec: The figure specification.
        lv: The levels, reference first.
        paired: Whether the design is paired.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    x, y, id_col = spec["data"]["x"], spec["data"]["y"], spec["data"].get("id")
    y_lab = spec["appearance"].get("y_label") or y
    if paired:
        pick = (f'    m = d[d["{x}"] == lv[0]][["{id_col}", "{y}"]].merge(\n'
                f'        d[d["{x}"] == g][["{id_col}", "{y}"]], on="{id_col}")\n'
                f'    m["{id_col}"] = m["{id_col}"].astype(str)\n'
                f'    m = m.sort_values("{id_col}", kind="stable")\n'
                f'    r, v = m["{y}_x"].to_numpy(float), m["{y}_y"].to_numpy(float)\n')
    else:
        pick = (f'    r = d.loc[d["{x}"] == lv[0], "{y}"].to_numpy(float)\n'
                f'    v = d.loc[d["{x}"] == g, "{y}"].to_numpy(float)\n')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. It
# recomputes the BCa bootstrap intervals with the same Park-Miller resamples and
# draws a simple version of the estimation plot.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

M, A, SEED, B = 2147483647, 48271, 20261029, 5000


def uniforms(n):
    out, s = np.empty(n), SEED
    for i in range(n):
        s = (s * A) % M
        out[i] = s
    return out / M


def bca(r, g, paired):
    per = len(r) if paired else len(r) + len(g)
    u = uniforms(B * per).reshape(B, per)
    if paired:
        i = np.floor(u * len(r)).astype(int)
        boot = (g[i] - r[i]).sum(axis=1) / len(r)
        dd = g - r
        jack = (dd.sum() - dd) / (len(dd) - 1)
    else:
        ir = np.floor(u[:, :len(r)] * len(r)).astype(int)
        ig = np.floor(u[:, len(r):] * len(g)).astype(int)
        boot = g[ig].sum(axis=1) / len(g) - r[ir].sum(axis=1) / len(r)
        jack = np.concatenate([g.mean() - (r.sum() - r) / (len(r) - 1),
                               (g.sum() - g) / (len(g) - 1) - r.mean()])
    est = g.sum() / len(g) - r.sum() / len(r)
    p0 = np.sum(boot < est - 1e-11 * np.abs(np.concatenate([r, g])).max()) / B
    z0 = stats.norm.ppf(min(max(p0, 0.5 / B), 1 - 0.5 / B))
    dj = jack.mean() - jack
    acc = np.sum(dj ** 3) / (6 * np.sum(dj ** 2) ** 1.5)
    zq = stats.norm.ppf([0.025, 0.975])
    return est, *np.quantile(boot, stats.norm.cdf(z0 + (z0 + zq) / (1 - acc * (z0 + zq))))


d = pd.read_csv("{in_name}")
d = d[d["{x}"].notna() & np.isfinite(pd.to_numeric(d["{y}"], errors="coerce"))]
d["{x}"] = d["{x}"].astype(str)
lv = {lv!r}
res = []
for g in lv[1:]:
{pick}    est, lo, hi = bca(r, v, {paired})
    p = (stats.ttest_rel(v, r) if {paired} else stats.ttest_ind(v, r, equal_var=False)).pvalue
    res.append((g, est, lo, hi, p))
    print(f"{{g}}: mean difference {{est:.6g}} (95% CI {{lo:.6g}} to {{hi:.6g}}), P = {{p:.4g}}")
fig, (a1, a2) = plt.subplots(2, 1, figsize=(5, 6))
for i, g in enumerate(lv):
    vals = d.loc[d["{x}"] == g, "{y}"]
    a1.scatter(np.full(len(vals), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(vals)), vals)
a1.set_xticks(range(len(lv))); a1.set_xticklabels(lv); a1.set_ylabel("{y_lab}")
for i, (g, est, lo, hi, p) in enumerate(res, start=1):
    a2.errorbar([i], [est], yerr=[[est - lo], [hi - est]], fmt="o", color="black")
a2.axhline(0, ls="--", color="grey"); a2.set_xlim(a1.get_xlim())
a2.set_xticks(range(len(lv))); a2.set_xticklabels([""] + lv[1:])
a2.set_ylabel("Mean difference vs {lv[0]}")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
