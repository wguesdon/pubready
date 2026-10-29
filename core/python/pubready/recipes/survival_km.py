"""Recipe: Kaplan-Meier survival curves compared between groups (lifelines).
Step curves with censoring ticks, the log-rank p-value on the plot, and a
number-at-risk table below. Mirrors the R survival_km."""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from matplotlib.lines import Line2D

from ..theme import apply_pub_style, longest_line, palette, wrap_labels
from ..util import p_phrase

LEGEND_WRAP = 24
"""The legend labels wrap at this many characters, as in the R engine."""


def km_break_step(tmax):
    """The step between the x breaks, as km_break_step() of the R engine.

    Args:
        tmax: The largest time.

    Returns:
        float: The smallest value of the 1-2-5 sequence that cuts the follow-up
        into at most 6 intervals.
    """
    if not np.isfinite(tmax) or tmax <= 0:
        return 1.0
    cand = sorted(m * 10.0 ** k for k in range(-3, 7) for m in (1, 2, 5))
    return next(c for c in cand if tmax / c <= 6 + 1e-9)


def km_breaks(tmax, step):
    """The breaks from 0 to the last time, by ``step``.

    Args:
        tmax: The largest time.
        step: The step between breaks.

    Returns:
        list: The breaks.
    """
    return [round(step * k, 10) for k in range(int(np.floor(tmax / step + 1e-9)) + 1)]


def km_at_risk(time, group, levels, breaks):
    """Subjects at risk at each break: those whose time is the break or later.

    This is the definition of survival::summary.survfit(times = t) in R. The
    lifelines helper add_at_risk_counts counts the subjects with a time after the
    break instead, which gives smaller numbers whenever an event falls on it.

    Args:
        time: The follow-up times.
        group: The group of each subject.
        levels: The groups, in the order of the table.
        breaks: The time points.

    Returns:
        numpy.ndarray: One row per group and one column per break.
    """
    time = np.asarray(time, dtype=float)
    group = np.asarray(group).astype(str)
    return np.array([[int(np.sum((group == g) & (time >= b))) for b in breaks] for g in levels])


def km_width(leg_labs, title):
    """The figure width, as km_width() of the R engine.

    Args:
        leg_labs: The wrapped legend labels.
        title: The legend title.

    Returns:
        float: 4.6 in for the panel and the axis, plus the legend.
    """
    chars = max([longest_line(s) for s in leg_labs] + [len(str(title))])
    return 4.6 + 0.55 + 0.075 * chars


def km_table_share(k):
    """The share of the figure height that the risk table takes, as in R.

    Args:
        k: The count of groups, one row each.

    Returns:
        float: The share.
    """
    return 0.08 + 0.045 * k


def _risk_table(axr, risk, lv, pal, breaks, tmax):
    """Draw the number-at-risk table on the x scale of the curves.

    A short line in the colour of each group starts its row; the legend names
    the groups, as in the R engine.

    Args:
        axr: The axes of the table, on the x limits of the curves.
        risk: The at-risk counts, one row per group.
        lv: The groups.
        pal: The colours of the groups.
        breaks: The time points.
        tmax: The largest time.
    """
    k = len(lv)
    for i in range(k):
        row = k - i
        for j, b in enumerate(breaks):
            axr.text(b, row, str(int(risk[i, j])), ha="center", va="center", fontsize=9)
        axr.plot([-0.13 * tmax, -0.08 * tmax], [row, row], color=pal[i], lw=3,
                 solid_capstyle="butt", clip_on=False)
    axr.set_xlim(-0.05 * tmax, 1.05 * tmax)
    axr.set_ylim(0.5, k + 0.5)
    axr.set_title("Number at risk", fontsize=10, fontweight="bold", loc="left")
    for side in ("top", "right", "left", "bottom"):
        axr.spines[side].set_visible(False)
    axr.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def recipe_survival_km(df, spec):
    """Draw Kaplan-Meier curves with the log-rank test and a number-at-risk table.

    Args:
        df: The tidy table, one row per subject.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing or the event column is not coded 1/0.
    """
    t, e, g = spec["data"]["time"], spec["data"]["event"], spec["data"]["x"]
    if not t or not e:
        raise ValueError("survival_km needs --time and --event columns")
    if not g:
        raise ValueError("survival_km needs a grouping column via --x")
    for c in (t, e, g):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found")

    df = df.copy()
    df[t] = pd.to_numeric(df[t], errors="coerce")
    ev = pd.to_numeric(df[e], errors="coerce")
    if not ev.dropna().isin([0, 1]).all():
        raise ValueError(f"event column '{e}' must be coded 1 = event, 0 = censored")
    df[e] = ev
    n0 = len(df)
    df = df.dropna(subset=[t, e, g])
    df = df[df[t] >= 0].reset_index(drop=True)
    clean_steps = (f"Dropped {n0 - len(df)} row(s) with missing/invalid values ({n0} -> {len(df)} rows)."
                   if len(df) < n0 else "No cleaning applied; input used as-is.")
    df[g] = df[g].astype(str)
    lv = sorted(df[g].unique())
    pal = palette(len(lv), spec["appearance"].get("palette"))

    # The x axis runs to the last time, so the last event is on the figure, with
    # the same breaks and the same at-risk counts (time >= t) as the R engine.
    tmax = float(df[t].max())
    step = km_break_step(tmax)
    breaks = km_breaks(tmax, step)
    risk = km_at_risk(df[t], df[g], lv, breaks)
    leg_labs = wrap_labels(lv, LEGEND_WRAP)

    w, h = km_width(leg_labs, g), 5.8
    fig = plt.figure(figsize=(w, h))
    share = km_table_share(len(lv))
    gs = fig.add_gridspec(2, 1, height_ratios=[1 - share, share], hspace=0.32)
    ax = fig.add_subplot(gs[0])
    # No sharex: pandas, which lifelines plots through, hides the x labels of an
    # axes that shares x with an axes below it. The table gets the same limits.
    axr = fig.add_subplot(gs[1])
    rows = []
    for i, grp in enumerate(lv):
        sub = df[df[g] == grp]
        kmf = KaplanMeierFitter()
        kmf.fit(sub[t], sub[e], label=grp)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i],
                                   show_censors=True, censor_styles={"marker": "|", "ms": 6})
        med = float(kmf.median_survival_time_)
        rows.append({"group": grp, "n": int(len(sub)), "events": int(sub[e].sum()),
                     "median_survival": med if np.isfinite(med) else float("nan")})

    res = multivariate_logrank_test(df[t], df[g], df[e])
    lr_p = float(res.p_value)
    # The log-rank P goes in a subtitle, where no curve can cross it, as in R.
    ax.set_title("Log-rank P < 0.001" if lr_p < 0.001 else f"Log-rank P = {lr_p:.2g}",
                 fontsize=10, color="#4d4d4d")
    ax.set_xlabel(spec["appearance"].get("x_label") or "Time")
    ax.set_ylabel(spec["appearance"].get("y_label") or "Survival probability")
    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.05 * tmax, 1.05 * tmax)
    ax.set_xticks(breaks)
    # The legend sits right of the panel, where it never covers a curve.
    handles = [Line2D([0], [0], color=pal[i], lw=1.6) for i in range(len(lv))]
    leg = ax.legend(handles, leg_labs, title=g, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, fontsize=9, title_fontsize=10)
    leg._legend_box.align = "left"  # the title sits over the keys, as in ggplot2
    apply_pub_style(ax, spec["appearance"].get("theme"))
    _risk_table(axr, risk, lv, pal, breaks, tmax)
    fig.subplots_adjust(left=0.95 / w, right=1 - (w - 4.45) / w, top=1 - 0.35 / h,
                        bottom=0.12 / h)

    stats_df = pd.DataFrame(rows)
    stats_df["log_rank_p"] = lr_p
    lr_df = len(lv) - 1
    test_meta = {
        "name": "Kaplan-Meier estimate with log-rank test",
        "log_rank": {"chisq": float(res.test_statistic), "df": lr_df, "p_value": lr_p},
        "groups": rows,
        "at_risk": {"times": breaks,
                    "counts": {grp: [int(v) for v in risk[i]] for i, grp in enumerate(lv)}},
    }
    methods = _methods(rows, res, lr_df, lr_p, g)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "label": g, "width": w, "height": h,
        "build_script": functools.partial(_script, spec, t, e, g, lv, pal, step, w),
    }


def _median_text(m):
    """Format a median survival time, or "not reached".

    Args:
        m: The median, NaN when the curve stays above 0.5.

    Returns:
        str: The text.
    """
    return "not reached" if not np.isfinite(m) else f"{m:.3g}"


def _methods(rows, res, lr_df, lr_p, g):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        rows: The per-group summary.
        res: The log-rank test result.
        lr_df: The degrees of freedom of the log-rank test.
        lr_p: Its P value.
        g: The group column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    med = ", ".join(f"{r['group']} (median {_median_text(r['median_survival'])})" for r in rows)
    return (
        "Survival was estimated by the Kaplan-Meier method and compared between "
        f"groups of {g} with the log-rank test (chi-square({lr_df}) = {res.test_statistic:.2f}, "
        f"{p_phrase(lr_p)})."
        f" Median survival by group: {med}."
        " Censoring is marked with ticks and the number at risk is tabulated below the curves;"
        " it counts the subjects whose time is at or after each time point."
        f" Analyses were performed in Python {platform.python_version()} with lifelines {v('lifelines')}."
    )


def _script(spec, t, e, g, lv, pal, step, w, in_name, fig_stub):
    """Write the standalone script that redraws the curves and the risk table.

    Args:
        spec: The figure specification.
        t: The time column.
        e: The event column.
        g: The group column.
        lv: The groups.
        pal: Their colours.
        step: The step between the x breaks.
        w: The figure width.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

df = pd.read_csv("{in_name}"); df["{g}"] = df["{g}"].astype(str)
order, pal = {lv}, {pal}
tmax = df["{t}"].max()
breaks = [{step!r} * k for k in range(int(np.floor(tmax / {step!r} + 1e-9)) + 1)]
fig, (ax, axr) = plt.subplots(2, 1, figsize=({w:.2f}, 5.8),
                              gridspec_kw={{"height_ratios": [0.74, 0.26]}})
for i, grp in enumerate(order):
    sub = df[df["{g}"] == grp]
    kmf = KaplanMeierFitter().fit(sub["{t}"], sub["{e}"], label=grp)
    kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i], show_censors=True,
                               censor_styles={{"marker": "|", "ms": 6}})
    # At risk at t: the subjects whose time is t or later, as in R.
    for b in breaks:
        axr.text(b, len(order) - i, str(int((sub["{t}"] >= b).sum())), ha="center",
                 va="center", fontsize=9, color=pal[i])
res = multivariate_logrank_test(df["{t}"], df["{g}"], df["{e}"])
p = res.p_value
ax.set_title("Log-rank P < 0.001" if p < 0.001 else "Log-rank P = %.2g" % p, fontsize=10)
ax.set_xlabel("{spec['appearance'].get('x_label') or 'Time'}")
ax.set_ylabel("{spec['appearance'].get('y_label') or 'Survival probability'}")
ax.set_ylim(0, 1.02); ax.set_xlim(-0.05 * tmax, 1.05 * tmax); ax.set_xticks(breaks)
ax.legend(title="{g}", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
for s in ("top", "right"): ax.spines[s].set_visible(False)
axr.set_xlim(ax.get_xlim()); axr.set_ylim(0.5, len(order) + 0.5); axr.axis("off")
axr.set_title("Number at risk", loc="left", fontsize=10, fontweight="bold")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf", bbox_inches="tight")
print("Reproduced figure.")
'''
