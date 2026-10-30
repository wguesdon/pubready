"""Recipe: spider plot of the change from baseline of each patient over time.

Mirrors the R recipe. One line per patient (--id), the time on --x and the
measured value on --y. Each value becomes a percentage change from that
patient's own baseline, which is the row with the smallest time. Two dashed
reference lines carry the RECIST 1.1 thresholds, +20% (--pd_threshold) and -30%
(--pr_threshold).

The plot is descriptive, so there is no hypothesis test. The stats table holds
one row per patient with the baseline value, the best and the worst change and
the category those two put the patient in. The best and the worst are measured
over the visits after the baseline, because every patient reads zero at the
baseline itself. That categorisation reads the change from baseline alone, so
the column is named for what it measures: a full RECIST 1.1 assessment also
reads the nadir, the non-target lesions and any new lesion.
"""
import functools

import numpy as np
import pandas as pd

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import cap_first, level_labels

#: The order a response legend reads in, from the best outcome to the worst.
CATEGORIES = ("partial response", "stable disease", "progressive disease")


def categorise(best, worst, pr_cut, pd_cut):
    """Put one patient in a response category from the change from baseline.

    Args:
        best: The smallest percentage change the patient reached.
        worst: The largest percentage change the patient reached.
        pr_cut: The change at or below which a patient is a partial response.
        pd_cut: The change at or above which a patient is progressive disease.

    Returns:
        The category name. A patient who reaches both thresholds is reported by
        the best value, which is how a best overall response is summarised.
    """
    if best <= pr_cut:
        return "partial response"
    if worst >= pd_cut:
        return "progressive disease"
    return "stable disease"


def check_one_row_per_time(ids, times, recipe, xcol):
    """Stop when a patient has two rows or more at one time point.

    A second baseline row would else become a later visit. The message names
    the patients and the times. spider_response and waterfall share it, and the
    R engine has the same rule.

    Args:
        ids: The patient of each row.
        times: The time of each row.
        recipe: The recipe name, for the message.
        xcol: The time column, for the message.

    Raises:
        ValueError: If a patient has two rows at one time.
    """
    tt = level_labels(pd.Series(list(times)))
    keys = pd.Series([f"{i}\r{t}" for i, t in zip(ids, tt)])
    dup = keys.duplicated().to_numpy()
    if not dup.any():
        return
    dup_ids = [i for i, flag in zip(ids, dup) if flag]
    dup_tt = [t for t, flag in zip(tt, dup) if flag]
    items, seen = [], set()
    for i, t in zip(dup_ids, dup_tt):
        if i not in seen:
            seen.add(i)
            items.append(f"'{i}' at {xcol} = {t}")
    shown = ", ".join(items[:6]) + (f" and {len(items) - 6} more" if len(items) > 6 else "")
    raise ValueError(f"{recipe} needs one row per patient and time point, and {len(items)} "
                     f"patient(s) have two rows or more at one time: {shown}; average or remove "
                     "the extra rows first")


def check_numeric_time(values, col):
    """Stop when the time column holds no number, as check_numeric_time() in R.

    The message says how many of the values are text, for example "week 0".
    spider_response and waterfall share it.

    Args:
        values: The time column.
        col: The column name, for the message.

    Raises:
        ValueError: If no value of the column reads as a number.
    """
    num = pd.to_numeric(values, errors="coerce")
    if num.notna().any():
        return
    txt = values.astype(str).str.strip()
    text = values.notna() & (txt != "")
    if text.any():
        raise ValueError(f"column '{col}' holds no number ({int(text.sum())} of {len(values)} "
                         f"values are text, for example '{txt[text].iloc[0]}'), and the time "
                         "has to be numeric")
    raise ValueError(f"column '{col}' holds no number, and the time has to be numeric")


def recipe_spider_response(df, spec):
    """Draw the change from baseline of every patient over time.

    Args:
        df: The tidy table, one row per measurement.
        spec: The figure specification.

    Returns:
        The bundle payload: the figure, the per patient table, the metadata,
        the methods paragraph and the standalone script.

    Raises:
        ValueError: If a column is missing, if the time holds no number, or if
            no patient carries a non-zero baseline.
    """
    x, y = spec["data"].get("x"), spec["data"].get("y")
    idc, grp = spec["data"].get("id"), spec["data"].get("group")
    if not x or not y:
        raise ValueError("spider_response needs --x (the time) and --y (the measured value)")
    if not idc:
        raise ValueError("spider_response needs --id (the patient column)")
    cols = [idc, x, y] + ([grp] if grp else [])
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    d = df[cols].copy()
    check_numeric_time(d[x], x)
    d[x] = pd.to_numeric(d[x], errors="coerce")
    d[y] = pd.to_numeric(d[y], errors="coerce")
    n_rows0 = len(d)
    d = d.dropna()
    n_missing = n_rows0 - len(d)
    d[idc] = d[idc].astype(str)
    check_one_row_per_time(d[idc], d[x], "spider_response", x)
    d = d.sort_values([idc, x], kind="mergesort").reset_index(drop=True)

    baseline = d.groupby(idc, sort=False)[y].first()
    usable = baseline[(baseline != 0) & np.isfinite(baseline)].index
    dropped = len(baseline) - len(usable)
    rows_zero = int((~d[idc].isin(usable)).sum())
    d = d[d[idc].isin(usable)].reset_index(drop=True)
    if len(d) < 2:
        raise ValueError("spider_response needs at least one patient with a non-zero baseline")
    d["change"] = 100 * (d[y] - d[idc].map(baseline)) / d[idc].map(baseline)

    ap = spec["appearance"]
    pd_cut = ap.get("pd_threshold")
    pd_cut = 20.0 if pd_cut is None else float(pd_cut)
    pr_cut = ap.get("pr_threshold")
    pr_cut = -30.0 if pr_cut is None else float(pr_cut)

    # The response is measured over the visits after the baseline. Every patient
    # reads zero at the baseline, so a patient who only grows would otherwise
    # carry a best change of zero rather than the smallest growth they reached.
    post = d[d[idc].duplicated()]
    no_follow_up = sorted(set(d[idc]) - set(post[idc]))
    rows_no_follow = int(d[idc].isin(no_follow_up).sum())
    if no_follow_up:
        d = d[~d[idc].isin(no_follow_up)].reset_index(drop=True)
        post = post[~post[idc].isin(no_follow_up)]
    if post.empty:
        raise ValueError("spider_response needs at least one patient with a visit "
                         "after the baseline")

    visits = d.groupby(idc, sort=False)["change"].size()
    per = post.groupby(idc, sort=False).agg(
        best_change_pct=("change", "min"), worst_change_pct=("change", "max"),
        last_change_pct=("change", "last")).reset_index()
    per.insert(1, "baseline", per[idc].map(baseline))
    per.insert(2, "visits", per[idc].map(visits))
    per["category_from_change"] = [
        categorise(row.best_change_pct, row.worst_change_pct, pr_cut, pd_cut)
        for row in per.itertuples()]

    steps = [f"Kept {len(per)} patient(s) across {len(d)} measurement(s)."]
    # Each reason is counted on its own: rows with a missing cell, then the
    # patients (and their rows) with a zero baseline or with no later visit.
    if n_missing:
        what = "patient, time, value or group" if grp else "patient, time or value"
        steps.append(f"Dropped {n_missing} row(s) with a missing {what}.")
    if dropped:
        steps.append(f"Dropped {dropped} patient(s) ({rows_zero} row(s)) whose baseline was "
                     "zero, because a percentage change needs a non-zero baseline.")
    if no_follow_up:
        steps.append(f"Dropped {len(no_follow_up)} patient(s) ({rows_no_follow} row(s)) with no "
                     "visit after the baseline, because a response needs a second "
                     "measurement.")

    if grp:
        colour_of = d[[idc, grp]].drop_duplicates(idc).set_index(idc)[grp].astype(str)
        legend_title = grp
        levels = list(dict.fromkeys(colour_of.loc[per[idc]]))
    else:
        colour_of = per.set_index(idc)["category_from_change"]
        legend_title = "Response"
        seen = set(colour_of)
        levels = [name for name in CATEGORIES if name in seen]
    pal = palette(len(levels), ap.get("palette"))
    colour_index = {name: pal[i] for i, name in enumerate(levels)}

    w, h = 5.8, 4.2
    fig, ax = new_fig(w, h)
    ax.axhline(0, color="grey", lw=0.8)
    ax.axhline(pd_cut, color="grey", lw=0.8, ls="--")
    ax.axhline(pr_cut, color="grey", lw=0.8, ls="--")
    for patient, rows in d.groupby(idc, sort=False):
        colour = colour_index[colour_of.loc[patient]]
        ax.plot(rows[x], rows["change"], color=colour, lw=1.1, alpha=0.9, zorder=2)
        ax.scatter(rows[x], rows["change"], color=colour, s=14, alpha=0.9, zorder=3)
    handles = [ax.plot([], [], color=colour_index[name], lw=1.6, label=name)[0]
               for name in levels]
    ax.legend(handles=handles, title=legend_title, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)

    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or "Change from baseline (%)")
    if ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    stats_df = per.rename(columns={idc: "patient"})
    if grp:
        stats_df["group"] = colour_of.loc[per[idc]].to_numpy()
    stats_df = stats_df.sort_values("best_change_pct").reset_index(drop=True)

    counts = stats_df["category_from_change"].value_counts()
    test_meta = {
        "name": "Spider plot (change from baseline)",
        "n_patients": int(len(per)),
        "n_measurements": int(len(d)),
        "n_dropped_patients": int(dropped),
        "n_without_follow_up": len(no_follow_up),
        "pd_threshold": pd_cut,
        "pr_threshold": pr_cut,
        "partial_response": int(counts.get("partial response", 0)),
        "stable_disease": int(counts.get("stable disease", 0)),
        "progressive_disease": int(counts.get("progressive disease", 0)),
    }
    methods = _methods(len(per), len(d), pd_cut, pr_cut, test_meta, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": d, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, pd_cut=pd_cut, pr_cut=pr_cut),
    }


def _methods(n_patients, n_rows, pd_cut, pr_cut, meta, spec):
    """Write the methods paragraph.

    Args:
        n_patients: The count of patients drawn.
        n_rows: The count of measurements drawn.
        pd_cut: The progressive disease threshold.
        pr_cut: The partial response threshold.
        meta: The metadata, which carries the category counts.
        spec: The figure specification.

    Returns:
        The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or spec["data"]["y"]
    return (
        f"{cap_first(y_lab)} was expressed as a percentage change from each patient's own "
        f"baseline, which is the first measurement in time, and drawn as one line per patient "
        f"(n = {n_patients} patients, {n_rows} measurements). Dashed reference lines mark "
        f"{pd_cut:+g}% and {pr_cut:+g}%, the RECIST 1.1 thresholds for progressive disease and "
        f"for a partial response. On the change from baseline alone, "
        f"{meta['partial_response']} patient(s) reached a partial response, "
        f"{meta['stable_disease']} were stable and {meta['progressive_disease']} reached "
        "progressive disease. That categorisation reads the target lesion change from baseline "
        "only. A full RECIST 1.1 assessment also reads the nadir, the non-target lesions and any "
        "new lesion. The figure is descriptive, so no statistical test is applied. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(spec, in_name, fig_stub, pd_cut, pr_cut):
    """Write the standalone script that redraws the figure.

    Args:
        spec: The figure specification.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.
        pd_cut: The progressive disease threshold.
        pr_cut: The partial response threshold.

    Returns:
        The script, as one string.
    """
    x, y, idc = spec["data"]["x"], spec["data"]["y"], spec["data"]["id"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

raw = pd.read_csv("{in_name}")
raw["{x}"] = pd.to_numeric(raw["{x}"], errors="coerce")
raw["{y}"] = pd.to_numeric(raw["{y}"], errors="coerce")
raw = raw.dropna(subset=["{idc}", "{x}", "{y}"]).sort_values(["{idc}", "{x}"])
base = raw.groupby("{idc}")["{y}"].first()
raw["change"] = 100 * (raw["{y}"] - raw["{idc}"].map(base)) / raw["{idc}"].map(base)

fig, ax = plt.subplots(figsize=({5.8}, {4.2}))
ax.axhline(0, color="grey", lw=0.8)
ax.axhline({pd_cut}, color="grey", lw=0.8, ls="--")
ax.axhline({pr_cut}, color="grey", lw=0.8, ls="--")
for _, rows in raw.groupby("{idc}"):
    ax.plot(rows["{x}"], rows["change"], lw=1.1)
    ax.scatter(rows["{x}"], rows["change"], s=14)
ax.set_xlabel("{x}"); ax.set_ylabel("Change from baseline (%)")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
