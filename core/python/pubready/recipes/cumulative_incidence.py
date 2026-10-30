"""Recipe: cumulative incidence of competing risks.

Mirrors the R recipe, which calls ``cmprsk::cuminc``. The two routines of cmprsk
are restated here so that both engines report the same numbers:

* ``cif_curve`` follows the Fortran routine ``cinc``: the Aalen-Johansen
  cumulative incidence function (CIF) of one cause, as the step corners that
  cmprsk returns. Tied times are processed together, as one time point.
* ``gray_test`` follows the Fortran routines ``crstm`` and ``crst``: the score
  and the variance of Gray's k-sample test (rho = 0) of one cause across groups,
  with one stratum.

lifelines' AalenJohansenFitter is not used, because it jitters tied times.
"""
import functools
import inspect
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style, palette
from ..util import level_labels, level_order
from .roc_curve import level_strings

#: A group label longer than this is shown in the legend only, and the risk
#: table marks its row with a short line in the group colour.
KEY_CHARS = 14


def nice_step(span, target):
    """A step of 1, 2, 2.5 or 5 times a power of ten that cuts a span into intervals.

    Args:
        span: The length to cut.
        target: The approximate count of intervals.

    Returns:
        float: The step, matched to the R engine.
    """
    if not np.isfinite(span) or span <= 0:
        return 1.0
    raw = span / target
    mag = 10.0 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * mag >= raw * (1 - 1e-9):
            return m * mag
    return 10 * mag


def axis_breaks(tmax):
    """The x breaks: multiples of the step from 0 up to the largest follow-up time.

    Args:
        tmax: The largest follow-up time.

    Returns:
        tuple: The step and the array of breaks.
    """
    step = nice_step(tmax, 6)
    return step, np.arange(math.floor(tmax / step + 1e-9) + 1) * step


def time_points(tmax, step):
    """The three tabulated time points.

    Each is 25%, 50% or 75% of the largest follow-up time, rounded to a multiple
    of half the x step.

    Args:
        tmax: The largest follow-up time.
        step: The x step.

    Returns:
        list: The three time points.
    """
    half = step / 2
    return [math.floor(q * tmax / half + 0.5) * half for q in (0.25, 0.5, 0.75)]


def fmt_time(t):
    """Format a time point for a column name, as R's ``sprintf("%g")`` does.

    Args:
        t: The time point.

    Returns:
        str: The text.
    """
    return f"{t:g}"


def tick_labels(values, step):
    """Tick labels with the decimals of the step, as the R engine writes them.

    Args:
        values: The tick positions.
        step: The step between ticks.

    Returns:
        list: The labels, for example 0.0, 0.2, 0.4 or 0, 20, 40.
    """
    txt = f"{step:g}"
    dec = len(txt.split(".")[1]) if "." in txt else 0
    return [f"{v:.{dec}f}" for v in values]


def join_words(items):
    """Join items as "a, b and c".

    Args:
        items: The texts.

    Returns:
        str: The joined text.
    """
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def cif_curve(time, any_event, cause_event):
    """The Aalen-Johansen CIF of one cause in one group, as cmprsk's ``cinc``.

    Args:
        time: The follow-up times, sorted in increasing order.
        any_event: 1 when the subject failed from any cause, else 0.
        cause_event: 1 when the subject failed from the cause, else 0.

    Returns:
        tuple: The times and the estimates at the step corners. The curve starts
        at (0, 0), each jump time appears twice (before and after the jump), and
        the last corner is the largest follow-up time.
    """
    n = len(time)
    fk = 1.0
    rs = float(n)
    xs, fs = [0.0], [0.0]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and time[j + 1] == time[i]:
            j += 1
        nd1 = int(cause_event[i:j + 1].sum())
        nd = int(any_event[i:j + 1].sum())
        if nd > 0:
            fkn = fk * (rs - nd) / rs
            if nd1 > 0:
                xs += [float(time[i]), float(time[i])]
                fs += [fs[-1], fs[-1] + fk * nd1 / rs]
            fk = fkn
        rs = float(n - (j + 1))
        i = j + 1
    xs.append(float(time[n - 1]))
    fs.append(fs[-1])
    return np.array(xs), np.array(fs)


def cif_value_at(t, xs, fs):
    """The CIF at time t, as ``cmprsk::timepoints`` reads it.

    Args:
        t: The time point.
        xs: The corner times of the curve.
        fs: The corner estimates of the curve.

    Returns:
        float: The value of the right-continuous step function, or NaN past the
        last follow-up time of the group.
    """
    if t > xs[-1]:
        return np.nan
    return float(fs[np.nonzero(xs <= t)[0][-1]])


def gray_test(time, code, group, ng, rho=0.0):
    """Gray's k-sample test of one cause, restated from cmprsk's ``crst``.

    Args:
        time: The follow-up times.
        code: 0 = censored, 1 = failed from the cause, 2 = failed from another cause.
        group: The group of each subject, coded 0 to ng - 1.
        ng: The count of groups.
        rho: The power of the weight function.

    Returns:
        tuple: The scores of the first ng - 1 groups and their variance matrix.
    """
    order = np.argsort(time, kind="stable")
    y = np.asarray(time, dtype=float)[order]
    m = np.asarray(code, dtype=int)[order]
    ig = np.asarray(group, dtype=int)[order]
    n = len(y)
    ng1 = ng - 1
    rs = np.bincount(ig, minlength=ng).astype(float)
    s = np.zeros(ng1)
    v = np.zeros((ng1, ng1))
    f1m, f1 = np.zeros(ng), np.zeros(ng)
    skmm, skm = np.ones(ng), np.ones(ng)
    v3 = np.zeros(ng)
    v2 = np.zeros((ng1, ng))
    c = np.zeros((ng, ng))
    a = np.zeros((ng, ng))
    fm = 0.0
    f = 0.0
    ll = 0
    while True:
        lu = ll
        while lu + 1 < n and not y[lu + 1] > y[ll]:
            lu += 1
        d = np.zeros((3, ng), dtype=int)
        for i in range(ll, lu + 1):
            d[m[i], ig[i]] += 1
        nd1 = int(d[1].sum())
        nd2 = int(d[2].sum())
        if nd1 > 0 or nd2 > 0:
            f = _gray_step(d, nd1, nd2, rs, f1m, f1, skmm, skm, fm, rho,
                           s, v, v3, v2, c, a, ng, ng1)
        if lu >= n - 1:
            break
        for i in range(ll, lu + 1):
            rs[ig[i]] -= 1
        fm = f
        f1m[:] = f1
        skmm[:] = skm
        ll = lu + 1
    for i in range(ng1):
        for j in range(i + 1):
            for k in range(ng):
                v[i, j] += c[i, k] * c[j, k] * v3[k]
                v[i, j] += c[i, k] * v2[j, k]
                v[i, j] += c[j, k] * v2[i, k]
    for i in range(ng1):
        for j in range(i):
            v[j, i] = v[i, j]
    return s, v


def _gray_step(d, nd1, nd2, rs, f1m, f1, skmm, skm, fm, rho, s, v, v3, v2, c, a, ng, ng1):
    """Process one failure time of Gray's test, updating the arrays in place.

    Args:
        d: The counts censored, failed from the cause and failed from another
            cause at this time, by group.
        nd1: The count failed from the cause at this time.
        nd2: The count failed from another cause at this time.
        rs: The risk set size of each group.
        f1m: The left-continuous CIF of each group.
        f1: The right-continuous CIF of each group.
        skmm: The left-continuous Kaplan-Meier estimate of each group.
        skm: The right-continuous Kaplan-Meier estimate of each group.
        fm: The left-continuous pooled CIF.
        rho: The power of the weight function.
        s: The scores.
        v: The variance matrix, lower triangle.
        v3: A variance work array.
        v2: A variance work array.
        c: A work matrix.
        a: A work matrix.
        ng: The count of groups.
        ng1: The count of groups less one.

    Returns:
        float: The right-continuous pooled CIF.
    """
    tr = 0.0
    tq = 0.0
    for i in range(ng):
        if rs[i] <= 0:
            continue
        td = float(d[1, i] + d[2, i])
        skm[i] = skmm[i] * (rs[i] - td) / rs[i]
        f1[i] = f1m[i] + (skmm[i] * d[1, i]) / rs[i]
        tr += rs[i] / skmm[i]
        tq += rs[i] * (1 - f1m[i]) / skmm[i]
    f = fm + nd1 / tr
    fb = (1 - fm) ** rho
    for i in range(ng):
        a[i, i:] = 0.0
        if rs[i] <= 0:
            continue
        t1 = rs[i] / skmm[i]
        a[i, i] = fb * t1 * (1 - t1 / tr)
        if a[i, i] != 0:
            c[i, i] += a[i, i] * nd1 / (tr * (1 - fm))
        for j in range(i + 1, ng):
            if rs[j] <= 0:
                continue
            a[i, j] = -fb * t1 * rs[j] / (skmm[j] * tr)
            if a[i, j] != 0:
                c[i, j] += a[i, j] * nd1 / (tr * (1 - fm))
    for i in range(1, ng):
        for j in range(i):
            a[i, j] = a[j, i]
            c[i, j] = c[j, i]
    for i in range(ng1):
        if rs[i] <= 0:
            continue
        s[i] += fb * (d[1, i] - nd1 * rs[i] * (1 - f1m[i]) / (skmm[i] * tq))
    if nd1 > 0:
        for k in range(ng):
            if rs[k] <= 0:
                continue
            t4 = 1.0
            if skm[k] > 0:
                t4 = 1 - (1 - f) / skm[k]
            t5 = 1.0
            if nd1 > 1:
                t5 = 1 - (nd1 - 1) / (tr * skmm[k] - 1)
            t3 = t5 * skmm[k] * nd1 / (tr * rs[k])
            v3[k] += t4 * t4 * t3
            for i in range(ng1):
                t1 = a[i, k] - t4 * c[i, k]
                v2[i, k] += t1 * t4 * t3
                for j in range(i + 1):
                    t2 = a[j, k] - t4 * c[j, k]
                    v[i, j] += t1 * t2 * t3
    if nd2 > 0:
        for k in range(ng):
            if skm[k] <= 0 or d[2, k] <= 0:
                continue
            t4 = (1 - f) / skm[k]
            t5 = 1.0
            if d[2, k] > 1:
                t5 = 1 - (d[2, k] - 1.0) / (rs[k] - 1.0)
            t6 = rs[k]
            t3 = t5 * ((skmm[k] ** 2) * d[2, k]) / (t6 ** 2)
            v3[k] += t4 * t4 * t3
            for i in range(ng1):
                t1 = t4 * c[i, k]
                v2[i, k] -= t1 * t4 * t3
                for j in range(i + 1):
                    t2 = t4 * c[j, k]
                    v[i, j] += t1 * t2 * t3
    return f


def gray_statistic(s, v):
    """The chi-square statistic s' inv(V) s, or NaN when V is singular.

    cmprsk checks the rank of V with a QR decomposition; the rank here comes from
    the singular values with the same relative tolerance of 1e-7.

    Args:
        s: The scores.
        v: The variance matrix.

    Returns:
        float: The statistic.
    """
    sv = np.linalg.svd(v, compute_uv=False)
    if len(sv) == 0 or sv[-1] <= 1e-7 * sv[0]:
        return np.nan
    return float(s @ np.linalg.solve(v, s))


def fit_cause(cause, d, lv):
    """The CIF of one cause in each group and Gray's test across the groups.

    Args:
        cause: The event code.
        d: The cleaned table (time, status, group), sorted by time.
        lv: The group levels.

    Returns:
        dict: The cause, the curves of each group and the test.
    """
    is_event = (d["status"] != "0").to_numpy()
    is_cause = (d["status"] == cause).to_numpy()
    curves = {}
    for g in lv:
        sel = (d["group"] == g).to_numpy()
        curves[g] = cif_curve(d["time"].to_numpy()[sel], is_event[sel].astype(int),
                              is_cause[sel].astype(int))
    test = {"statistic": np.nan, "df": np.nan, "p_value": np.nan}
    if len(lv) > 1:
        code = np.where(~is_event, 0, np.where(is_cause, 1, 2))
        gidx = d["group"].map({g: i for i, g in enumerate(lv)}).to_numpy()
        s, v = gray_test(d["time"].to_numpy(float), code, gidx, len(lv))
        stat = gray_statistic(s, v)
        if np.isfinite(stat) and stat >= 0:
            df = float(len(lv) - 1)
            test = {"statistic": stat, "df": df, "p_value": float(stats.chi2.sf(stat, df))}
    return {"cause": cause, "curves": curves, "test": test}


def p_label(test):
    """The Gray's test line of a panel.

    Args:
        test: The test of one cause.

    Returns:
        str: The label.
    """
    if not np.isfinite(test["p_value"]):
        return "One group; no test"
    if test["p_value"] < 0.001:
        return "Gray's test P < 0.001"
    return f"Gray's test P = {test['p_value']:.2g}"


def legend_ncol(lv, width):
    """The legend columns: as many as fit across the figure width.

    Args:
        lv: The group levels.
        width: The figure width in inches.

    Returns:
        int: The count of columns, matched to the R engine.
    """
    per = 0.085 * max(len(g) for g in lv) + 0.7
    ncol = min(len(lv), max(1, math.floor((width - 0.5) / per)))
    # Both engines fill the legend column by column, so the column count is the
    # one that the row count implies.
    return int(math.ceil(len(lv) / math.ceil(len(lv) / ncol)))


def figure_width(n_causes, lv):
    """The figure width: one panel of 4.2 in per cause, wider for full row labels.

    Args:
        n_causes: The count of causes.
        lv: The group levels.

    Returns:
        float: The width in inches, matched to the R engine.
    """
    chars = max(len(g) for g in lv)
    extra = 0.0 if chars > KEY_CHARS else max(0.0, 0.075 * chars - 0.45)
    return n_causes * (4.2 + extra) + 0.2


def figure_dims(width, lv, legend_rows):
    """The figure size: taller with more groups and more legend rows.

    Args:
        width: The figure width in inches.
        lv: The group levels.
        legend_rows: The rows of the legend, or 0 without a legend.

    Returns:
        dict: The width, the height and the legend height, in inches.
    """
    legend_h = 0.45 + 0.27 * legend_rows if legend_rows > 0 else 0.0
    return {"width": width, "height": legend_h + 4.4 + 0.24 * len(lv), "legend_h": legend_h}


def recipe_cumulative_incidence(df, spec):
    """Draw the cumulative incidence of each cause with Gray's test.

    Args:
        df: The tidy table, one row per subject.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if no event is coded, or if --cause
            is not an event code.
    """
    dat = spec["data"]
    tcol, ecol, gcol = dat.get("time"), dat.get("event"), dat.get("x")
    if not tcol or not ecol:
        raise ValueError("cumulative_incidence needs --time and --event "
                         "(0 = censored, other values = the cause)")
    has_group = bool(gcol)
    for cn in [tcol, ecol] + ([gcol] if has_group else []):
        if cn not in df.columns:
            raise ValueError(f"column '{cn}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({
        "time": pd.to_numeric(df[tcol], errors="coerce"),
        "status": level_strings(df[ecol]).str.strip(),
        "group": level_labels(df[gcol]).str.strip() if has_group else "all",
    })
    keep = np.isfinite(d["time"]) & (d["time"] >= 0) & d["status"].notna() & d["group"].notna()
    d = d[keep].reset_index(drop=True)
    if len(d) < 2:
        raise ValueError("cumulative_incidence needs at least 2 subjects with a time and a status")
    causes = sorted(d.loc[d["status"] != "0", "status"].unique())
    if not causes:
        raise ValueError(f"column '{ecol}' holds no event (every value is 0 = censored)")

    note = []
    given = dat.get("cause")
    if given is not None and str(given).strip():
        cause = str(given).strip()
        if cause not in causes:
            raise ValueError(f"--cause '{cause}' is not an event code of '{ecol}' "
                             f"({', '.join(causes)})")
    elif "1" in causes:
        cause = "1"
    else:
        cause = causes[0]
        note.append("No --cause was given and no event has code 1, so the cause of interest "
                    f"is the first code in sorted order, {cause}.")
    causes = [cause] + [c for c in causes if c != cause]
    # The groups sort as numbers when every label is a number, else in code-point
    # order, by the level rule that both engines share (level_order).
    lv = level_order(d["group"])
    d = d.iloc[np.argsort(d["time"].to_numpy(), kind="stable")].reset_index(drop=True)
    n_drop = int((~keep).sum())
    counts = ", ".join([f"{int((d['status'] == '0').sum())} censored"] +
                       [f"{int((d['status'] == c).sum())} with {ecol} = {c}" for c in causes])
    clean_steps = ([f"Dropped {n_drop} row(s) with a missing or negative time, a missing "
                    "status or a missing group."] if n_drop else []) + \
        [f"Kept {len(d)} of {n0} row(s): {counts}."] + note

    fits = [fit_cause(c, d, lv) for c in causes]
    tmax = float(d["time"].max())
    step, brks = axis_breaks(tmax)
    tps = time_points(tmax, step)
    rows = []
    for f in fits:
        for g in lv:
            xs, fs = f["curves"][g]
            row = {"cause": f["cause"], "group": g, "n": int((d["group"] == g).sum()),
                   "n_events": int(((d["group"] == g) & (d["status"] == f["cause"])).sum())}
            for t in tps:
                row[f"cif_{fmt_time(t)}"] = cif_value_at(t, xs, fs)
            row.update({"gray_statistic": f["test"]["statistic"], "gray_df": f["test"]["df"],
                        "gray_p": f["test"]["p_value"]})
            rows.append(row)
    stats_df = pd.DataFrame(rows)

    ap = spec["appearance"]
    pal = palette(len(lv), ap.get("palette"))
    width = figure_width(len(causes), lv)
    ncol_leg = legend_ncol(lv, width)
    legend_rows = math.ceil(len(lv) / ncol_leg) if has_group else 0
    dims = figure_dims(width, lv, legend_rows)
    fig = draw_figure(fits, d, lv, pal, brks, tmax, ap, ecol, tcol,
                      gcol if has_group else None, ncol_leg, dims)

    name = ("Cumulative incidence (Aalen-Johansen) with Gray's test" if len(lv) > 1
            else "Cumulative incidence (Aalen-Johansen)")
    test_meta = {
        "name": name, "cause_of_interest": cause, "causes": causes, "groups": lv,
        "time_points": tps, "gray": {f["cause"]: f["test"] for f in fits},
    }
    curves = pd.DataFrame([{"cause": f["cause"], "group": g, "time": t, "est": e}
                           for f in fits for g in lv
                           for t, e in zip(*f["curves"][g])])
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(fits, d, lv, causes, tps, tcol, ecol, gcol, has_group),
        "df_used": d, "clean_steps": clean_steps, "label": slugify(ecol),
        "width": dims["width"], "height": dims["height"], "curves": curves,
        "build_script": functools.partial(_script, spec, has_group, tps, dims, lv),
    }


def y_axis(fits):
    """The y step and the top of the y axis, matched to the R engine.

    Args:
        fits: The fitted causes.

    Returns:
        tuple: The step and the top.
    """
    ymax = max(float(np.max(fs)) for f in fits for _, fs in f["curves"].values())
    ymax = max(ymax, 0.05)
    ystep = nice_step(ymax, 4)
    ytop = min(1.0, math.ceil(ymax / ystep - 1e-9) * ystep)
    return ystep, ytop


def draw_figure(fits, d, lv, pal, brks, tmax, ap, ecol, tcol, gcol, ncol_leg, dims):
    """Draw one panel per cause with the number at risk below it.

    Args:
        fits: The fitted causes.
        d: The cleaned table.
        lv: The group levels.
        pal: The colour of each group.
        brks: The x breaks.
        tmax: The largest follow-up time.
        ap: The appearance options.
        ecol: The status column.
        tcol: The time column.
        gcol: The grouping column, or None for one group.
        ncol_leg: The legend columns.
        dims: The figure size.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    ystep, ytop = y_axis(fits)
    yticks = np.arange(int(math.floor(ytop / ystep + 1e-9)) + 1) * ystep
    ng = len(lv)
    fig = plt.figure(figsize=(dims["width"], dims["height"]), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.06, h_pad=0.06, wspace=0.06, hspace=0.04)
    gs = fig.add_gridspec(2, len(fits), height_ratios=[3.2, 0.55 + 0.24 * ng])
    xlo, xhi = -0.02 * tmax, tmax * 1.04
    at_risk = {g: [int((d.loc[d["group"] == g, "time"] >= b).sum()) for b in brks] for g in lv}
    for k, f in enumerate(fits):
        ax = fig.add_subplot(gs[0, k])
        for i, g in enumerate(lv):
            xs, fs = f["curves"][g]
            ax.plot(xs, fs, color=pal[i], lw=1.6, solid_joinstyle="miter")
        ax.set_xlim(xlo, xhi)
        ax.set_xticks(brks)
        ax.set_xticklabels(tick_labels(brks, brks[1] - brks[0] if len(brks) > 1 else 1.0))
        ax.set_ylim(-0.01 * ytop, ytop * 1.03)
        ax.set_yticks(yticks)
        ax.set_yticklabels(tick_labels(yticks, ystep))
        ax.set_xlabel(ap.get("x_label") or tcol)
        ax.set_ylabel(ap.get("y_label") or "Cumulative incidence")
        ax.set_title(f"{ecol} = {f['cause']}", fontsize=12, fontweight="bold", pad=22)
        ax.annotate(p_label(f["test"]), xy=(0.5, 1.0), xycoords="axes fraction",
                    xytext=(0, 7), textcoords="offset points", ha="center", va="bottom",
                    fontsize=10, color="#4d4d4d")
        apply_pub_style(ax, ap.get("theme"))
        axt = fig.add_subplot(gs[1, k])
        for i, g in enumerate(lv):
            for b, nr in zip(brks, at_risk[g]):
                axt.text(b, ng - 1 - i, str(nr), ha="center", va="center", fontsize=9,
                         color=pal[i])
        axt.set_xlim(xlo, xhi)
        axt.set_ylim(-0.6, ng - 0.4)
        axt.set_yticks(range(ng))
        if max(len(g) for g in lv) > KEY_CHARS:
            axt.set_yticklabels([])
            for i in range(ng):
                axt.plot([-0.11 * tmax, -0.05 * tmax], [ng - 1 - i] * 2, color=pal[i], lw=2.2,
                         solid_capstyle="butt", clip_on=False)
        else:
            axt.set_yticklabels(list(reversed(lv)), fontsize=9)
        axt.set_xticks([])
        axt.tick_params(length=0, pad=14)
        for side in ("top", "right", "bottom", "left"):
            axt.spines[side].set_visible(False)
        axt.set_title("Number at risk", fontsize=10, loc="left", pad=4)
    if gcol:
        handles = [Line2D([0], [0], color=pal[i], lw=1.6) for i in range(ng)]
        fig.legend(handles, lv, title=gcol, loc="outside upper center", ncol=ncol_leg,
                   frameon=False, fontsize=11, title_fontsize=12)
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    return fig


def _methods(fits, d, lv, causes, tps, tcol, ecol, gcol, has_group):
    """Write the methods paragraph.

    Args:
        fits: The fitted causes.
        d: The cleaned table.
        lv: The group levels.
        causes: The event codes, the cause of interest first.
        tps: The tabulated time points.
        tcol: The time column.
        ecol: The status column.
        gcol: The grouping column.
        has_group: Whether a grouping column was given.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    grp = ""
    if has_group:
        sizes = "; ".join(f"{g}, n = {int((d['group'] == g).sum())}" for g in lv)
        grp = f" in each group of {gcol} ({sizes})"
    if has_group and len(lv) > 1:
        tests = (" The groups were compared for each event type with Gray's k-sample test "
                 "(rho = 0): " + "; ".join(_test_text(f, ecol) for f in fits) + ".")
    else:
        tests = " With one group, no test was applied."
    return (
        f"The cumulative incidence of each event type ({ecol} = {', '.join(causes)}; 0 = "
        f"censored) was estimated with the Aalen-Johansen estimator{grp}, which treats the "
        "other event types as competing events rather than as censoring."
        f" The event type of interest is {ecol} = {causes[0]}.{tests}"
        f" The cumulative incidence at {tcol} = {join_words(fmt_time(t) for t in tps)} is "
        "given in the stats table, and the number at risk (the subjects with a follow-up "
        "time of t or more) is shown below each panel."
        f" Computed in Python {platform.python_version()} with numpy {v('numpy')}, restating "
        f"the cinc and crstm routines of the R package cmprsk, and drawn with matplotlib "
        f"{v('matplotlib')}."
    )


def _test_text(f, ecol):
    """Describe the test of one cause.

    Args:
        f: The fitted cause.
        ecol: The status column.

    Returns:
        str: The text.
    """
    t = f["test"]
    if not np.isfinite(t["p_value"]):
        return f"{ecol} = {f['cause']}, not estimable"
    p = "< 0.001" if t["p_value"] < 0.001 else f"= {t['p_value']:.2g}"
    return (f"{ecol} = {f['cause']}, chi-square({int(t['df'])}) = {t['statistic']:.2f}, "
            f"P {p}")


def _script(spec, has_group, tps, dims, lv, in_name, fig_stub):
    """Write the standalone script that recomputes the curves and Gray's test.

    Args:
        spec: The figure specification.
        has_group: Whether a grouping column was given.
        tps: The tabulated time points.
        dims: The figure size.
        lv: The group levels, in order.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    dat = spec["data"]
    tcol, ecol, gcol = dat["time"], dat["event"], dat.get("x")
    helpers = "\n\n".join(inspect.getsource(fn) for fn in
                          (cif_curve, gray_test, _gray_step, gray_statistic))
    grp_line = (f'raw["{gcol}"].astype(str).str.replace(r"\\.0$", "", regex=True).str.strip()'
                if has_group
                else '"all"')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The
# functions below restate the cinc and crstm routines of the R package cmprsk.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

{helpers}

raw = pd.read_csv("{in_name}")
st = raw["{ecol}"].astype(str).str.replace(r"\\.0$", "", regex=True).str.strip()
d = pd.DataFrame({{"time": pd.to_numeric(raw["{tcol}"], errors="coerce"), "status": st,
                   "group": {grp_line}}}).dropna()
d = d[d["time"] >= 0].sort_values("time", kind="stable").reset_index(drop=True)
lv = {lv!r}
causes = sorted(d.loc[d["status"] != "0", "status"].unique())
fig, axes = plt.subplots(1, len(causes), figsize=({dims["width"]:.2f}, {dims["height"]:.2f}), squeeze=False)
for ax, cause in zip(axes[0], causes):
    ev, cs = (d["status"] != "0").to_numpy(), (d["status"] == cause).to_numpy()
    for g in lv:
        sel = (d["group"] == g).to_numpy()
        xs, fs = cif_curve(d["time"].to_numpy()[sel], ev[sel].astype(int), cs[sel].astype(int))
        ax.plot(xs, fs, label=g)
        print(cause, g, [float(fs[np.nonzero(xs <= t)[0][-1]]) if t <= xs[-1] else None
                         for t in {list(tps)}])
    if len(lv) > 1:
        code = np.where(~ev, 0, np.where(cs, 1, 2))
        s, v = gray_test(d["time"].to_numpy(float), code, d["group"].map(lv.index).to_numpy(), len(lv))
        stat = gray_statistic(s, v)
        print("Gray's test", cause, stat, stats.chi2.sf(stat, len(lv) - 1))
    ax.set_title("{ecol} = " + cause); ax.set_xlabel("{tcol}"); ax.set_ylabel("Cumulative incidence")
    ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
