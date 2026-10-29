"""Recipe: Cox proportional-hazards model with a forest plot of hazard ratios
(lifelines). Mirrors the R cox_forest: HR with 95% CI and p per covariate, and a
proportional-hazards check.

The figure matches the R engine row for row: the label column on the left, one
row per model term (the first level of each factor gets its own "reference"
row), squares with 95% CI whiskers on a log axis, and the "HR (95% CI)" and "P"
columns on the right.
"""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from matplotlib import transforms
from matplotlib.ticker import NullLocator

from ..theme import apply_pub_style, palette, plain_log_axis
from ..util import format_p, p_phrase

# Tick sequences of the hazard-ratio axis, from coarse to fine (mantissas per decade).
TICK_SETS = [(1, 2, 5), (1, 1.5, 2, 3, 5, 7), (1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8)]


def forest_stars(p):
    """The significance stars of the forest plot, the same in both engines.

    The rule stops at three stars, as survminer::ggforest graded them, so a P
    below 1e-4 keeps "***".

    Args:
        p: The P value of one term.

    Returns:
        str: "***", "**", "*" or an empty string.
    """
    if p is None or np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def ticks_for(mant, lo, hi):
    """The values of one tick sequence in [lo, hi], extended across 1 when needed.

    On a side of 1 that the CIs reach but no tick covers, the next value of the
    sequence outward is added.

    Args:
        mant: The mantissas of the sequence in one decade, for example (1, 2, 5).
        lo: The smallest value on the axis, at most 1.
        hi: The largest value on the axis, at least 1.

    Returns:
        list: The tick positions in increasing order.
    """
    seq = sorted(m * 10.0 ** k for k in range(-4, 5) for m in mant)
    ticks = [v for v in seq if lo <= v <= hi]
    below = [v for v in seq if v < lo]
    above = [v for v in seq if v > hi]
    if lo < 1 and not any(v < 1 for v in ticks) and below:
        ticks = [below[-1]] + ticks
    if hi > 1 and not any(v > 1 for v in ticks) and above:
        ticks = ticks + [above[0]]
    return ticks


def forest_log_ticks(lo, hi):
    """Ticks of the hazard-ratio axis that cover the CIs on both sides of 1.

    The 1-2-5 sequence comes first. It is kept when it covers the CIs on both
    sides of 1 without widening the axis by more than 60%. A narrow range falls
    back to finer steps, and a range over many decades to the powers of ten. The
    R engine uses the same rule.

    Args:
        lo: The smallest lower CI bound.
        hi: The largest upper CI bound.

    Returns:
        list: The tick positions in increasing order, 1 included.
    """
    lo, hi = min(lo, 1.0), max(hi, 1.0)
    span = np.log(hi / lo)
    sets = TICK_SETS
    if len(ticks_for((1, 2, 5), lo, hi)) > 8:
        sets = [(1,)]
    ticks = [1.0]
    for mant in sets:
        ticks = ticks_for(mant, lo, hi)
        total = np.log(max([hi] + ticks) / min([lo] + ticks))
        if total <= 1.6 * span:
            return ticks
    return ticks


def label_inches(labels):
    """The width of the label column in inches, from the longest label.

    Args:
        labels: The row labels.

    Returns:
        float: The width, the same estimate as the R engine.
    """
    return 0.3 + 0.085 * max([len(s) for s in labels] + [4])


def design_and_rows(df, t, e, covs):
    """Build the numeric design matrix and the display rows of the forest plot.

    A categorical covariate is one-hot encoded with its first sorted level as the
    reference, and gives a "reference" row plus one row per other level. A
    numeric covariate gives one row.

    Args:
        df: The cleaned data.
        t: The time column.
        e: The event column.
        covs: The covariate columns.

    Returns:
        tuple: The design data frame and a list of row dicts with the keys
        ``label`` and ``term`` (None for a reference row).
    """
    design = df[[t, e]].copy()
    rows = []
    for c in covs:
        s = df[c]
        if s.dtype == object or str(s.dtype).startswith("category") or not np.issubdtype(s.dtype, np.number):
            cats = sorted(s.astype(str).unique())
            rows.append({"label": f"{c}: {cats[0]}", "term": None})
            for lvl in cats[1:]:
                design[f"{c}{lvl}"] = (s.astype(str) == lvl).astype(int)
                rows.append({"label": f"{c}: {lvl}", "term": f"{c}{lvl}"})
        else:
            design[c] = pd.to_numeric(s, errors="coerce")
            rows.append({"label": c, "term": c})
    return design, rows


def row_texts(row, stats):
    """The "HR (95% CI)" and "P" texts of one display row.

    Args:
        row: The row dict.
        stats: The stats table indexed by term.

    Returns:
        tuple: The HR text and the P text.
    """
    if row["term"] is None:
        return "reference", ""
    r = stats.loc[row["term"]]
    star = forest_stars(r["p_value"])
    p_txt = format_p(r["p_value"]) + (" " + star if star else "")
    return f"{r['HR']:.2f} ({r['ci_low']:.2f}–{r['ci_high']:.2f})", p_txt


def cox_subtitle(n, events, ph_p):
    """Write the subtitle: the counts and, when the PH check fails, a warning.

    The warning keeps the figure from showing the hazard ratios of a model
    that breaks the proportional-hazards assumption as valid.

    Args:
        n: The number of subjects.
        events: The number of events.
        ph_p: The P value of the proportional-hazards check, or NaN.

    Returns:
        str: The subtitle.
    """
    txt = f"{n} subjects, {events} events"
    if ph_p is not None and np.isfinite(ph_p) and ph_p < 0.05:
        txt += f"; PH assumption violated ({p_phrase(ph_p)})"
    return txt


def forest_plot(rows, stats_df, spec, n, events, ph_p=np.nan):
    """Draw the forest plot with the label column and the two text columns.

    Args:
        rows: The display rows from ``design_and_rows``.
        stats_df: The stats table.
        spec: The figure specification.
        n: The number of subjects.
        events: The number of events.
        ph_p: The P value of the proportional-hazards check, or NaN.

    Returns:
        dict: The figure, its ticks, the label width, the width and the height.
    """
    ap = spec["appearance"]
    title = ap.get("title")
    pal = palette(1, ap.get("palette"))
    stats = stats_df.set_index("term")
    k = len(rows)
    yy = np.arange(k, 0, -1).astype(float)
    est = [(y_, stats.loc[r["term"]]) for y_, r in zip(yy, rows) if r["term"] is not None]
    est = [(y_, s) for y_, s in est if np.isfinite(s["HR"])]
    lows = [s["ci_low"] for _, s in est if np.isfinite(s["ci_low"]) and s["ci_low"] > 0]
    highs = [s["ci_high"] for _, s in est if np.isfinite(s["ci_high"])]
    lo = min(lows) if lows else 1.0
    hi = max(highs) if highs else 1.0
    ticks = forest_log_ticks(lo, hi)
    xa, xb = min([lo] + ticks), max([hi] + ticks)
    r = np.log10(xb / xa) if xb > xa else np.log10(2.0)
    ylim = (0.4, k + 1.6)

    lab_in = label_inches([row["label"] for row in rows])
    w = lab_in + 3.0 + 2.5
    h = max(3.0, 0.3 * (k + 1) + 1.6 + (0.3 if title else 0.0))
    top_in = 0.45 + (0.3 if title else 0.0)
    bottom_in = 0.85
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([lab_in / w, bottom_in / h, 2.85 / w, (h - top_in - bottom_in) / h])
    axt = fig.add_axes([(lab_in + 3.0) / w, bottom_in / h, 2.3 / w, (h - top_in - bottom_in) / h])

    ax.axvline(1.0, ls=":", color="0.4", lw=0.9, zorder=1)
    for y_, s in est:
        ax.hlines(y_, s["ci_low"], s["ci_high"], color="0.15", lw=1.0, zorder=2)
    ax.scatter([s["HR"] for _, s in est], [y_ for y_, _ in est], marker="s", s=48,
               color=pal[0], zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(10 ** (np.log10(xa) - 0.05 * r), 10 ** (np.log10(xb) + 0.05 * r))
    ax.set_xticks(ticks)
    plain_log_axis(ax.xaxis)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(*ylim)
    ax.set_yticks([])
    ax.set_xlabel(ap.get("x_label") or "Hazard ratio (95% CI, log scale)")
    ax.set_title(cox_subtitle(n, events, ph_p), fontsize=9.5, color="#4d4d4d", loc="left")
    apply_pub_style(ax, ap.get("theme"))
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10.4)

    trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
    for y_, row in zip(yy, rows):
        fig.text(0.15 / w, y_, row["label"], transform=trans, ha="left", va="center",
                 fontsize=10.5)

    axt.set_xlim(0, 1)
    axt.set_ylim(*ylim)
    axt.axis("off")
    for y_, row in zip(yy, rows):
        hr_txt, p_txt = row_texts(row, stats)
        axt.text(0, y_, hr_txt, ha="left", va="center", fontsize=10)
        axt.text(0.66, y_, p_txt, ha="left", va="center", fontsize=10)
    for x_, head in ((0, "HR (95% CI)"), (0.66, "P")):
        axt.text(x_, k + 1, head, ha="left", va="center", fontsize=10, fontweight="bold")
    if title:
        fig.suptitle(title, fontweight="bold", fontsize=13)
    return {"fig": fig, "ticks": ticks, "label_in": lab_in, "width": w, "height": h}


def recipe_cox_forest(df, spec):
    t, e = spec["data"]["time"], spec["data"]["event"]
    covs = spec["data"].get("covariates")
    if not t or not e:
        raise ValueError("cox_forest needs --time and --event columns")
    if not covs:
        covs = [c for c in df.columns if c not in (t, e)]
    for c in [t, e] + covs:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found")

    df = df.copy()
    df[t] = pd.to_numeric(df[t], errors="coerce")
    ev = pd.to_numeric(df[e], errors="coerce")
    if not ev.dropna().isin([0, 1]).all():
        raise ValueError(f"event column '{e}' must be coded 1 = event, 0 = censored")
    df[e] = ev
    n0 = len(df)
    df = df.dropna(subset=[t, e] + covs)
    df = df[df[t] >= 0].reset_index(drop=True)
    clean_steps = (f"Dropped {n0 - len(df)} row(s) with missing/invalid values ({n0} -> {len(df)} rows)."
                   if len(df) < n0 else "No cleaning applied; input used as-is.")

    # Build a numeric design matrix; one-hot encode categoricals (drop first = reference).
    design, rows = design_and_rows(df, t, e, covs)

    cph = CoxPHFitter()
    # A tight precision makes lifelines converge to the estimate of survival::coxph;
    # the default stops about 2e-6 short of it. Terms are named as R names them.
    cph.fit(design, duration_col=t, event_col=e, fit_options={"precision": 1e-11})
    summ = cph.summary
    terms = list(summ.index)
    hr = summ["exp(coef)"].to_numpy()
    lo = summ["exp(coef) lower 95%"].to_numpy()
    hi = summ["exp(coef) upper 95%"].to_numpy()
    pv = summ["p"].to_numpy()

    ph_p = np.nan
    try:
        from lifelines.statistics import proportional_hazard_test
        ph = proportional_hazard_test(cph, design, time_transform="rank")
        ph_p = float(ph.summary["p"].min())
    except Exception:
        pass

    stats_df = pd.DataFrame({"term": terms, "HR": hr, "ci_low": lo, "ci_high": hi,
                             "z": summ["z"].to_numpy(), "p_value": pv})
    n_subj, n_events = int(len(df)), int(df[e].sum())
    fp = forest_plot(rows, stats_df, spec, n_subj, n_events, ph_p)

    test_meta = {
        "name": "Cox proportional-hazards model",
        "n": n_subj, "events": n_events, "covariates": covs,
        "proportional_hazards": {"test": "proportional_hazard_test (min p)", "p_value": ph_p},
        "hazard_ratios": [{"term": tt, "HR": float(h_), "ci_low": float(l_),
                           "ci_high": float(u_), "p_value": float(p_)}
                          for tt, h_, l_, u_, p_ in zip(terms, hr, lo, hi, pv)],
    }
    methods = _methods(stats_df, ph_p, covs, n_subj, n_events)
    return {
        "fig": fp["fig"], "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "label": "hazard_ratios",
        "width": fp["width"], "height": fp["height"],
        "build_script": functools.partial(_script, spec, t, e, covs, fp["ticks"],
                                          fp["label_in"], fp["width"], fp["height"]),
    }


def _methods(stats_df, ph_p, covs, n, events):
    from importlib.metadata import version as v
    import platform
    sig = stats_df[stats_df["p_value"] < 0.05]
    if len(sig):
        parts = "; ".join(f"{r.term} (HR = {r.HR:.2f}, 95% CI {r.ci_low:.2f}-{r.ci_high:.2f}, {p_phrase(r.p_value)})"
                          for r in sig.itertuples())
        sig_txt = f" Significant predictors: {parts}."
    else:
        sig_txt = " No individual predictor reached significance."
    ph_txt = "" if np.isnan(ph_p) else (
        f" The proportional-hazards assumption was {'violated' if ph_p < 0.05 else 'not violated'}"
        f" (minimum covariate cox.zph {p_phrase(ph_p)}).")
    return (
        f"Hazard ratios were estimated with a Cox proportional-hazards model ({n} subjects, {events} events) "
        f"including {', '.join(covs)}." + sig_txt + ph_txt +
        f" Analyses were performed in Python {platform.python_version()} with lifelines {v('lifelines')}"
        f", and the forest plot was drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, t, e, covs, ticks, lab_in, w, h, in_name, fig_stub):
    """Write the standalone script that refits the model and redraws the figure.

    Args:
        spec: The figure specification.
        t: The time column.
        e: The event column.
        covs: The covariate columns.
        ticks: The tick positions of the hazard-ratio axis.
        lab_in: The width of the label column in inches.
        w: The figure width in inches.
        h: The figure height in inches.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script text.
    """
    ap = spec["appearance"]
    colour = palette(1, ap.get("palette"))[0]
    xlab = ap.get("x_label") or "Hazard ratio (95% CI, log scale)"
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from matplotlib import transforms
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter

df = pd.read_csv("{in_name}")
covs = {covs!r}
design = df[["{t}", "{e}"]].copy()
rows = []  # (label, term); term None marks the reference level of a factor
for c in covs:
    s = df[c]
    if s.dtype == object or not np.issubdtype(s.dtype, np.number):
        cats = sorted(s.astype(str).unique())
        rows.append(("%s: %s" % (c, cats[0]), None))
        for lvl in cats[1:]:
            design["%s%s" % (c, lvl)] = (s.astype(str) == lvl).astype(int)
            rows.append(("%s: %s" % (c, lvl), "%s%s" % (c, lvl)))
    else:
        design[c] = pd.to_numeric(s, errors="coerce")
        rows.append((c, c))
cph = CoxPHFitter().fit(design, "{t}", "{e}", fit_options={{"precision": 1e-11}})
s = cph.summary
print(s)

k = len(rows); yy = np.arange(k, 0, -1).astype(float)
ticks = {[float(x) for x in ticks]!r}
w, h, lab_in = {w:.4f}, {h:.4f}, {lab_in:.4f}
fig = plt.figure(figsize=(w, h))
ax = fig.add_axes([lab_in / w, 0.85 / h, 2.85 / w, (h - 1.3) / h])
axt = fig.add_axes([(lab_in + 3.0) / w, 0.85 / h, 2.3 / w, (h - 1.3) / h])
ax.axvline(1.0, ls=":", color="0.4", lw=0.9)
trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
for y_, (label, term) in zip(yy, rows):
    fig.text(0.15 / w, y_, label, transform=trans, ha="left", va="center", fontsize=10.5)
    if term is None:
        axt.text(0, y_, "reference", ha="left", va="center", fontsize=10)
        continue
    hr_, lo_, hi_, p_ = s.loc[term, ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    ax.hlines(y_, lo_, hi_, color="0.15", lw=1.0)
    ax.scatter([hr_], [y_], marker="s", s=48, color="{colour}", zorder=3)
    star = "***" if p_ < 0.001 else "**" if p_ < 0.01 else "*" if p_ < 0.05 else ""
    p_txt = ("< 0.001" if p_ < 0.001 else "%.2g" % p_) + (" " + star if star else "")
    axt.text(0, y_, "%.2f (%.2f\\u2013%.2f)" % (hr_, lo_, hi_), ha="left", va="center", fontsize=10)
    axt.text(0.66, y_, p_txt, ha="left", va="center", fontsize=10)
axt.text(0, k + 1, "HR (95% CI)", ha="left", va="center", fontsize=10, fontweight="bold")
axt.text(0.66, k + 1, "P", ha="left", va="center", fontsize=10, fontweight="bold")
ax.set_xscale("log")
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FuncFormatter("{{:g}}".format))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_ylim(0.4, k + 1.6); axt.set_ylim(0.4, k + 1.6); axt.set_xlim(0, 1); axt.axis("off")
ax.set_yticks([]); ax.set_xlabel({xlab!r}, fontsize=13)
ax.set_title("%d subjects, %d events" % (len(df), int(df["{e}"].sum())), fontsize=9.5,
             color="#4d4d4d", loc="left")
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
fig.savefig("{fig_stub}.pdf", bbox_inches="tight")
print("Reproduced figure.")
'''
