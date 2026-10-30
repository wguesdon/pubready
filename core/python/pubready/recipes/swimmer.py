"""Recipe: swimmer plot of the time each patient spends on treatment or on study.

Mirrors the R recipe. Input is one row per patient: the patient (--id), the
duration (--time), the arm (--group, optional) and the status (--event,
optional: 1 = the course ended, 0 = ongoing). One horizontal bar per patient,
sorted from the longest course at the top, with an arrow at the end of each
ongoing course. Columns named in --markers hold the time of an event on the
course, and each one is drawn with its own shape. The plot is descriptive, so
there is no hypothesis test.
"""
import functools

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from ..bundle import slugify
from ..theme import apply_pub_style, legend_width, new_fig, palette, wrap_label, wrap_labels
from ..util import level_labels
from .survival_km import event_codes

#: The marker shapes, in the order of --markers, matched to the R engine.
SHAPES = ("o", "^", "D", "s", "v")

#: The legend labels wrap at this many characters, as in the R engine.
LEGEND_WRAP = 24


def swimmer_width(patients, arms, markers, title):
    """The figure width, as swimmer_width() of the R engine.

    3.8 in for the panel, the patient labels (when they are drawn, up to 40
    patients) at 0.055 in per character, and the legend of the arms and the
    markers, and never less than 6.4 in. A long label widens the figure instead
    of squeezing the panel.

    Args:
        patients: The patient labels.
        arms: The arms of the legend, or an empty list.
        markers: The marker columns.
        title: The legend title.

    Returns:
        float: The width in inches.
    """
    ychars = max([len(str(p)) for p in patients] + [0]) if len(patients) <= 40 else 0
    leg = wrap_labels(list(arms) + list(markers), LEGEND_WRAP)
    return max(6.4, 3.8 + 0.5 + 0.055 * ychars + legend_width(leg, title))


def arm_summary(d, g):
    """Summarise one arm.

    Args:
        d: The per-patient table.
        g: The arm.

    Returns:
        dict: The count of patients, the median duration and the count ongoing.
    """
    s = d[d["group"] == g]
    return {"n": int(len(s)), "median_duration": float(s["duration"].median()),
            "ongoing": int(s["ongoing"].sum())}


def marker_step(d, m):
    """A data log step for the values of one marker that fall outside the bar.

    A value after the end of the bar or below 0 is logged. The symbols stay
    where the data puts them. The R engine writes the same step.

    Args:
        d: The per-patient table.
        m: The marker column.

    Returns:
        str or None: The step, or None when every value is inside its bar.
    """
    v = d[m]
    out = (np.isfinite(v) & ((v > d["duration"]) | (v < 0))).to_numpy()
    if not out.any():
        return None
    vals = level_labels(v[out].reset_index(drop=True))
    bars = level_labels(d["duration"][out].reset_index(drop=True))
    items = [f"'{p}' ({a}, bar {b})" for p, a, b in zip(d["patient"][out], vals, bars)]
    shown = ", ".join(items[:6]) + (f" and {len(items) - 6} more" if len(items) > 6 else "")
    return (f"{int(out.sum())} value(s) of '{m}' fall after the end of the bar or below 0, and "
            f"the symbols are drawn where the data puts them: {shown}.")


def recipe_swimmer(df, spec):
    """Draw one bar per patient for the duration of the course.

    Args:
        df: The tidy table, one row per patient.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if the status is not coded 1/0, or
            if a patient has two rows.
    """
    dat = spec["data"]
    idc, tcol, ev, grp = dat.get("id"), dat.get("time"), dat.get("event"), dat.get("group")
    markers = list(dat.get("markers") or [])
    if not idc:
        raise ValueError("swimmer needs --id (the patient column)")
    if not tcol:
        raise ValueError("swimmer needs --time (the duration on treatment)")
    has_event, has_group = bool(ev), bool(grp)
    for c in [idc, tcol] + ([ev] if has_event else []) + ([grp] if has_group else []) + markers:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    status = (event_codes(df[ev], ev, "1 = the course ended, 0 = ongoing") if has_event
              else pd.Series(1.0, index=df.index))
    d = pd.DataFrame({
        "patient": df[idc].astype(str).where(df[idc].notna()),
        "duration": pd.to_numeric(df[tcol], errors="coerce"),
        "group": df[grp].astype(str).where(df[grp].notna()) if has_group else "all",
        "ongoing": (status == 0).where(status.notna()),
    })
    for m in markers:
        d[m] = pd.to_numeric(df[m], errors="coerce")
    # A row with a missing cell and a row with a negative duration are counted
    # on their own, so the data log gives the right reason for each drop.
    missing = (d["patient"].isna() | ~np.isfinite(d["duration"]) | d["group"].isna()
               | d["ongoing"].isna())
    negative = ~missing & (d["duration"] < 0)
    keep = ~missing & ~negative
    d = d[keep].reset_index(drop=True)
    d["ongoing"] = d["ongoing"].astype(bool)
    if d["patient"].duplicated().any():
        raise ValueError("swimmer needs one row per patient")
    if d.empty:
        raise ValueError("swimmer needs at least one patient with a duration")
    steps = []
    if missing.any():
        steps.append(f"Dropped {int(missing.sum())} row(s) with a missing patient, duration, "
                     "status or arm.")
    if negative.any():
        steps.append(f"Dropped {int(negative.sum())} row(s) with a negative duration.")
    steps.append(f"Kept {len(d)} of {n0} patient(s).")
    steps += [s for s in (marker_step(d, m) for m in markers) if s]

    # The arms are sorted, as the other recipes sort their groups; R sorts them
    # in the same C-locale order.
    g_lv = sorted(pd.unique(d["group"]))
    d = d.sort_values(["duration", "patient"], kind="mergesort").reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)
    ap = spec["appearance"]
    pal = palette(len(g_lv), ap.get("palette"))
    colour = {g: pal[i] for i, g in enumerate(g_lv)}
    span = float(d["duration"].max()) or 1.0

    w = swimmer_width(list(d["patient"]), g_lv if has_group else [], markers,
                      grp if has_group else "")
    h = max(3.6, min(10.0, 1.6 + 0.18 * len(d)))
    fig, ax = new_fig(w, h)
    ax.barh(d["rank"], d["duration"], height=0.6, color=[colour[g] for g in d["group"]])
    for _, r in d[d["ongoing"]].iterrows():
        ax.annotate("", xy=(r["duration"] + 0.04 * span, r["rank"]), xytext=(r["duration"], r["rank"]),
                    arrowprops={"arrowstyle": "-|>", "color": "0.12", "lw": 0.8})
    handles = ([Line2D([0], [0], color=colour[g], lw=6, label=wrap_label(g, LEGEND_WRAP))
                for g in g_lv] if has_group else [])
    for i, m in enumerate(markers):
        sel = np.isfinite(d[m])
        if sel.any():
            ax.scatter(d.loc[sel, m], d.loc[sel, "rank"], marker=SHAPES[i % len(SHAPES)], s=30,
                       facecolor="white", edgecolor="black", linewidths=0.9, zorder=3)
        handles.append(Line2D([0], [0], marker=SHAPES[i % len(SHAPES)], ls="", markersize=6,
                              markerfacecolor="white", markeredgecolor="black",
                              label=wrap_label(m, LEGEND_WRAP)))
    if len(d) <= 40:
        ax.set_yticks(d["rank"])
        ax.set_yticklabels(d["patient"])
    else:
        ax.set_yticks([])
    ax.set_ylim(0.2, len(d) + 0.8)
    ax.set_xlim(0, span * 1.08)
    if handles:
        ax.legend(handles=handles, title=grp if has_group else None, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9, title_fontsize=10)
    ax.set_title(f"{len(d)} patients; {int(d['ongoing'].sum())} ongoing (arrow)",
                 fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or tcol)
    ax.set_ylabel(ap.get("y_label") or "Patient")
    apply_pub_style(ax, ap.get("theme"))
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()

    stats_df = d.iloc[::-1][["patient", "group", "duration", "ongoing"] + markers].reset_index(drop=True)
    if not has_group:
        stats_df = stats_df.drop(columns="group")
    test_meta = {
        "name": "Swimmer plot (duration per patient)",
        "n_patients": int(len(d)), "n_ongoing": int(d["ongoing"].sum()),
        "arms": {g: arm_summary(d, g) for g in g_lv},
        "markers": {m: int(np.isfinite(d[m]).sum()) for m in markers},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(d, g_lv, markers, tcol, has_event),
        "df_used": d, "clean_steps": steps, "label": slugify(tcol), "width": w, "height": h,
        "build_script": functools.partial(_script, spec, has_event, has_group, w),
    }


def _methods(d, g_lv, markers, tcol, has_event):
    """Write the methods paragraph.

    Args:
        d: The per-patient table.
        g_lv: The arms, in order.
        markers: The marker columns.
        tcol: The duration column.
        has_event: Whether the input carries a status column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    arms = "; ".join(
        f"{'' if g == 'all' else g + ': '}{s['n']} patient(s), median duration "
        f"{s['median_duration']:.3g}, {s['ongoing']} ongoing"
        for g, s in ((g, arm_summary(d, g)) for g in g_lv))
    mk = ""
    if markers:
        mk = " Symbols mark " + ", ".join(
            f"{m} ({int(np.isfinite(d[m]).sum())} patient(s))" for m in markers) + "."
    return (
        f"Each of the {len(d)} patients is drawn as one horizontal bar whose length is the value of "
        f"{tcol}, sorted from the "
        "longest course at the top."
        + (" An arrow marks a course that was ongoing at the data cut-off." if has_event else "")
        + f"{mk} {arms}. The figure is descriptive, so no statistical test is applied. Drawn in "
        f"Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(spec, has_event, has_group, w, in_name, fig_stub):
    """Write the standalone script that redraws the bars.

    Args:
        spec: The figure specification.
        has_event: Whether the input carries a status column.
        has_group: Whether the input carries an arm column.
        w: The figure width in inches.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    dat = spec["data"]
    idc, tcol, ev, grp = dat["id"], dat["time"], dat.get("event"), dat.get("group")
    ongoing = f'raw["{ev}"] == 0' if has_event else "False"
    colour = f'raw["{grp}"].astype("category").cat.codes' if has_group else "0"
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}").sort_values("{tcol}").reset_index(drop=True)
raw["ongoing"] = {ongoing}
raw["colour"] = {colour}
fig, ax = plt.subplots(figsize=({w:.2f}, max(3.6, min(10, 1.6 + 0.18 * len(raw)))))
ax.barh(np.arange(1, len(raw) + 1), raw["{tcol}"], color=plt.cm.tab10(raw["colour"]))
ongoing = raw[raw["ongoing"]]
ax.scatter(ongoing["{tcol}"], ongoing.index + 1, marker=">", color="black")
ax.set_yticks(np.arange(1, len(raw) + 1)); ax.set_yticklabels(raw["{idc}"])
ax.set_xlabel("{tcol}"); ax.set_ylabel("Patient")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
