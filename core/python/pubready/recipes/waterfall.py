"""Recipe: waterfall plot of the best response of each patient. Mirrors the R recipe.

One bar per patient, sorted from the largest growth on the left to the deepest
shrinkage on the right, with dashed lines at the RECIST 1.1 thresholds (+20% and
-30% by default). With --x (the time), the data is one row per visit and each
patient's best change is the smallest percentage change from their own baseline
over the visits after it; the category follows spider_response. Without --x, the
data is one row per patient and --y is already the best percentage change. The
plot is descriptive, so there is no hypothesis test.
"""
import functools

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from .spider_response import CATEGORIES, categorise, check_numeric_time


def from_visits(d):
    """Best and worst change from baseline per patient from one row per visit.

    Args:
        d: The cleaned data with id, time, y and group columns.

    Returns:
        tuple: The per-patient table, the count of patients with a zero
        baseline, and the count without a visit after the baseline.
    """
    d = d.sort_values(["id", "time"], kind="mergesort").reset_index(drop=True)
    base = d.groupby("id", sort=True)["y"].first()
    usable = base.index[np.isfinite(base) & (base != 0)]
    dropped_zero = len(base) - len(usable)
    d = d[d["id"].isin(usable)].reset_index(drop=True)
    d["change"] = 100 * (d["y"] - d["id"].map(base)) / d["id"].map(base)
    post = d[d["id"].duplicated()]
    no_follow = len(set(d["id"]) - set(post["id"]))
    per = post.groupby("id", sort=True)["change"].agg(["min", "max"]).reset_index()
    per.columns = ["patient", "best_change_pct", "worst_change_pct"]
    per["group"] = per["patient"].map(d.groupby("id", sort=True)["group"].first())
    return per, dropped_zero, no_follow


def recipe_waterfall(df, spec):
    """Draw the best change from baseline of every patient as sorted bars.

    Args:
        df: The tidy table, one row per visit (with --x) or per patient.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, or if a patient has two rows when
            there is no time column.
    """
    dat = spec["data"]
    x, y, idc, grp = dat.get("x"), dat.get("y"), dat.get("id"), dat.get("group")
    if not y:
        raise ValueError("waterfall needs --y (the measured value, or the best percentage change)")
    if not idc:
        raise ValueError("waterfall needs --id (the patient column)")
    has_time, has_group = bool(x), bool(grp)
    for c in [idc, y] + ([x] if has_time else []) + ([grp] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")
    ap = spec["appearance"]
    pd_cut = 20.0 if ap.get("pd_threshold") is None else float(ap["pd_threshold"])
    pr_cut = -30.0 if ap.get("pr_threshold") is None else float(ap["pr_threshold"])

    if has_time:
        check_numeric_time(df[x], x)
    n0 = len(df)
    d = pd.DataFrame({"id": df[idc].astype(str).where(df[idc].notna()),
                      "time": pd.to_numeric(df[x], errors="coerce") if has_time else 0.0,
                      "y": pd.to_numeric(df[y], errors="coerce"),
                      "group": df[grp].astype(str).where(df[grp].notna()) if has_group else "all"})
    d = d.dropna().reset_index(drop=True)
    n_missing = n0 - len(d)
    extra = []
    if has_time:
        per, dropped_zero, no_follow = from_visits(d)
        per["category"] = [categorise(b, w, pr_cut, pd_cut)
                           for b, w in zip(per["best_change_pct"], per["worst_change_pct"])]
        if dropped_zero:
            extra.append(f"Dropped {dropped_zero} patient(s) whose baseline was zero, because a "
                         "percentage change needs a non-zero baseline.")
        if no_follow:
            extra.append(f"Dropped {no_follow} patient(s) with no visit after the baseline.")
    else:
        if d["id"].duplicated().any():
            raise ValueError("without --x, waterfall needs one row per patient; pass --x (the time) "
                             "for one row per visit")
        per = pd.DataFrame({"patient": d["id"], "best_change_pct": d["y"], "group": d["group"]})
        per["category"] = [categorise(b, b, pr_cut, pd_cut) for b in per["best_change_pct"]]
    if per.empty:
        raise ValueError("waterfall needs at least one patient with a response")
    steps = ([f"Dropped {n_missing} row(s) with a missing patient, time, value or group."]
             if n_missing else []) + extra + [f"Kept {len(per)} patient(s)."]

    per["neg"] = -per["best_change_pct"]
    per = per.sort_values(["neg", "patient"], kind="mergesort").drop(columns="neg")
    per = per.reset_index(drop=True)
    per["rank"] = np.arange(1, len(per) + 1)
    if has_group:
        levels, fill_of, legend_title = list(pd.unique(d["group"])), per["group"], grp
    else:
        levels = [c for c in CATEGORIES if c in set(per["category"])]
        fill_of, legend_title = per["category"], "Response"
    pal = palette(len(levels), ap.get("palette"))
    colour = {lv: pal[i] for i, lv in enumerate(levels)}
    n_pr = int((per["category"] == "partial response").sum())

    w, h = max(5.0, min(9.0, 3.2 + 0.16 * len(per))), 4.4
    fig, ax = new_fig(w, h)
    ax.bar(per["rank"], per["best_change_pct"], width=0.8,
           color=[colour[v] for v in fill_of])
    ax.axhline(0, color="0.12", lw=0.8)
    for cut in (pd_cut, pr_cut):
        ax.axhline(cut, color="0.45", lw=0.8, ls="--")
    if len(per) <= 30:
        ax.set_xticks(per["rank"])
        ax.set_xticklabels(per["patient"], rotation=90, fontsize=7)
    else:
        ax.set_xticks([])
    ax.set_xlim(0.4, len(per) + 0.6)
    handles = [Patch(facecolor=colour[lv], label=lv) for lv in levels]
    ax.legend(handles=handles, title=legend_title, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.set_title(f"{len(per)} patients; {n_pr} ({100 * n_pr / len(per):.0f}%) reach the "
                 f"{pr_cut:+g}% line", fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    ax.set_xlabel(ap.get("x_label") or "Patient")
    ax.set_ylabel(ap.get("y_label") or "Best change from baseline (%)")
    apply_pub_style(ax, ap.get("theme"))
    ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()

    cols =(["patient", "group", "best_change_pct"] + (["worst_change_pct"] if has_time else [])
            + ["category"])
    stats_df = per[cols].rename(columns={"category": "category_from_change"})
    if not has_group:
        stats_df = stats_df.drop(columns="group")
    counts = {}
    for g in pd.unique(per["group"]):
        vc = per.loc[per["group"] == g, "category"].value_counts()
        counts[g] = {c.replace(" ", "_"): int(vc.get(c, 0)) for c in CATEGORIES}
    test_meta = {"name": "Waterfall plot (best change from baseline)", "n_patients": int(len(per)),
                 "pd_threshold": pd_cut, "pr_threshold": pr_cut,
                 "input": "one row per visit" if has_time else "one row per patient",
                 "counts": counts}
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(per, has_time, pd_cut, pr_cut, y),
        "df_used": d, "clean_steps": steps, "label": slugify(y), "width": w, "height": h,
        "build_script": functools.partial(_script, spec, has_time, pd_cut, pr_cut),
    }


def _methods(per, has_time, pd_cut, pr_cut, y):
    """Write the methods paragraph.

    Args:
        per: The per-patient table.
        has_time: Whether the input is one row per visit.
        pd_cut: The progressive disease threshold.
        pr_cut: The partial response threshold.
        y: The value column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    vc = per["category"].value_counts()
    src = (f"The best response of each patient was the smallest percentage change in {y} from that "
           "patient's own baseline over the visits after it. " if has_time else
           f"The best percentage change of each patient was read from {y}. ")
    return (
        f"{src}The {len(per)} patients are drawn as bars sorted from the largest increase to the "
        f"deepest decrease, with dashed reference lines at {pd_cut:+g}% and {pr_cut:+g}%, the "
        "RECIST 1.1 thresholds for progressive disease and for a partial response. On the change "
        f"from baseline alone, {int(vc.get('partial response', 0))} patient(s) reached a partial "
        f"response, {int(vc.get('stable disease', 0))} were stable and "
        f"{int(vc.get('progressive disease', 0))} reached progressive disease. A full RECIST 1.1 "
        "assessment also reads the non-target lesions and any new lesion. The figure is "
        f"descriptive, so no statistical test is applied. Drawn in Python "
        f"{platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(spec, has_time, pd_cut, pr_cut, in_name, fig_stub):
    """Write the standalone script that redraws the waterfall.

    Args:
        spec: The figure specification.
        has_time: Whether the input is one row per visit.
        pd_cut: The progressive disease threshold.
        pr_cut: The partial response threshold.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    d = spec["data"]
    x, y, idc = d.get("x"), d["y"], d["id"]
    if has_time:
        body = f'''raw = raw.dropna(subset=["{idc}", "{x}", "{y}"]).sort_values(["{idc}", "{x}"])
base = raw.groupby("{idc}")["{y}"].first()
raw["change"] = 100 * (raw["{y}"] - raw["{idc}"].map(base)) / raw["{idc}"].map(base)
best = raw[raw["{idc}"].duplicated()].groupby("{idc}")["change"].min()'''
    else:
        body = f'best = raw.set_index("{idc}")["{y}"]'
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
{body}
best = best.sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(6, 4.4))
ax.bar(np.arange(1, len(best) + 1), best.to_numpy())
ax.axhline(0, color="k"); ax.axhline({pd_cut}, ls="--", color="grey"); ax.axhline({pr_cut}, ls="--", color="grey")
ax.set_xlabel("Patient"); ax.set_ylabel("Best change from baseline (%)")
for side in ("top", "right"): ax.spines[side].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print(best)
print("Reproduced figure.")
'''
