"""Recipe: the composition of each sample, compared between conditions.

Mirrors core/r/R/composition.R. Input is long: one row per sample (--id) and
category (--x) with a count (--y), and an optional condition (--group). Without
--y each row is one observation. The counts are summed per sample and category
and turned into proportions within each sample.

The figure is one 100% stacked bar per sample, ordered by group and then by
name, one strip per group. With two groups each category is compared with the
Mann-Whitney U test on the per-sample proportions (the p-value rule of R's
wilcox.test), with three or more with the Kruskal-Wallis test, and the P values
are adjusted across the categories with the Benjamini-Hochberg method.
"""
import functools
import re
import textwrap

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style, palette
from ..util import mann_whitney_r
from .distribution import as_labels


def wrap_label(s, width):
    """Wrap a label on spaces, the same lines as R's strwrap.

    Args:
        s: The label.
        width: The longest line wanted, in characters.

    Returns:
        str: The wrapped label.
    """
    return textwrap.fill(str(s), width, break_on_hyphens=False, break_long_words=False)


def bh_adjust(p):
    """Benjamini-Hochberg adjustment with the rule of R's p.adjust.

    Missing P values stay missing and do not count in the number of tests.

    Args:
        p: The P values.

    Returns:
        numpy.ndarray: The adjusted P values.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = len(q)
    if n == 0:
        return out
    o = np.argsort(-q, kind="stable")
    i = np.arange(n, 0, -1)
    adj = np.minimum(1.0, np.minimum.accumulate(n / i * q[o]))
    res = np.empty(n)
    res[o] = adj
    out[ok] = res
    return out


def natural_key(s):
    """Sort key of a sample name: the text before a trailing number, the number,
    then the full name, so that S2 comes before S10, as in the R engine.

    Args:
        s: The sample name.

    Returns:
        tuple: The key.
    """
    m = re.search(r"[0-9]+$", s)
    prefix = s[:m.start()] if m else s
    return (prefix.encode("utf-8"), m is None, float(m.group()) if m else 0.0, s.encode("utf-8"))


def signif(p):
    """Significance stars of an adjusted P value, empty when it is 0.05 or more.

    Args:
        p: The adjusted P value.

    Returns:
        str: The stars.
    """
    if p is None or not np.isfinite(p):
        return ""
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def test_category(v, g, lv):
    """Compare the per-sample proportions of one category between the groups.

    Args:
        v: The proportion of the category in each sample.
        g: The group of each sample.
        lv: The group levels, in order.

    Returns:
        tuple: The statistic and the P value, NaN when every proportion is equal.
    """
    if len(np.unique(v)) < 2:
        return float("nan"), float("nan")
    if len(lv) == 2:
        return mann_whitney_r(v[g == lv[0]], v[g == lv[1]])
    res = stats.kruskal(*[v[g == gl] for gl in lv])
    return float(res.statistic), float(res.pvalue)


def legend_label(cat, p_adj):
    """The legend label of one category, wrapped at 14 characters with stars.

    Args:
        cat: The category.
        p_adj: Its adjusted P value.

    Returns:
        str: The label.
    """
    stars = signif(p_adj)
    return wrap_label(cat, 14) + (f" {stars}" if stars else "")


def dims(ids, g_of, lv, labels):
    """The figure size and label layout, with the rules of the R engine.

    Args:
        ids: The sample names, in plot order.
        g_of: The group of each sample.
        lv: The group levels.
        labels: The legend labels.

    Returns:
        dict: The width, the height, the rotation, the legend columns and the
        strip width of each group.
    """
    n = len(ids)
    slot = 0.3 if n <= 20 else 0.22
    longest = max(len(s) for s in ids)
    rotate = n > 6 or longest > 8
    lab_h = 0.075 * min(longest, 40) if rotate else 0.25
    leg_lines = sum(lab.count("\n") + 1 for lab in labels)
    leg_h = 0.19 * leg_lines + 0.08 * len(labels) + 0.4
    height = max(4.4, 3.4 + lab_h + (0.8 if n > 20 else 0.0))
    leg_cols = 2 if leg_h > height - 0.6 else 1
    if leg_cols == 2:
        height = max(height, leg_h / 2 + 1.0)
    width = max(4.8, 1.0 + slot * n + 0.12 * (len(lv) - 1) + leg_cols * 1.45 + 0.3)
    strip_chars = [max(8, int(np.floor(int(np.sum(g_of == gl)) * slot / 0.08))) for gl in lv]
    return {"width": width, "height": height, "rotate": rotate, "leg_cols": leg_cols,
            "small_x": n > 20, "strip_chars": strip_chars}


def recipe_composition(df, spec):
    """Describe and compare the composition of each sample.

    Args:
        df: The tidy data.
        spec: The figure specification.

    Returns:
        dict: The recipe result for the bundle writer.
    """
    data = spec["data"]
    id_col, cat_col, cnt, grp = data.get("id"), data.get("x"), data.get("y"), data.get("group")
    has_count, has_group = bool(cnt), bool(grp)
    if not id_col or not cat_col:
        raise ValueError("composition needs --id (the sample) and --x (the category); "
                         "--y (a count) and --group are optional")
    for c in [id_col, cat_col] + ([cnt] if has_count else []) + ([grp] if has_group else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    steps = []
    n0 = len(df)
    d = pd.DataFrame({"id": df[id_col], "category": df[cat_col],
                      "group": df[grp] if has_group else "all"})
    if has_count:
        yv = pd.to_numeric(df[cnt], errors="coerce")
        coerced = int((yv.isna() & df[cnt].notna()).sum())
        if coerced:
            steps.append(f"Coerced column '{cnt}' to numeric; {coerced} value(s) became NA.")
        d["count"] = yv.astype(float)
    else:
        d["count"] = 1.0
        steps.append("No count column was given, so each row counts as one observation.")
    keep = (d["id"].notna() & d["category"].notna() & d["group"].notna()
            & np.isfinite(d["count"].to_numpy(dtype=float)))
    if (~keep).sum():
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing sample, category, group "
                     f"or count ({n0} -> {int(keep.sum())} rows).")
    d = d[keep].reset_index(drop=True)
    for c in ("id", "category", "group"):
        d[c] = as_labels(d[c])
    if (d["count"] < 0).any():
        raise ValueError("composition needs counts of zero or more")
    if has_group:
        per = d[["id", "group"]].drop_duplicates()
        dup = per["id"].duplicated()
        if dup.any():
            raise ValueError(f"each sample must belong to one group; sample "
                             f"'{per['id'][dup].iloc[0]}' has more than one")

    cats = list(pd.unique(d["category"]))
    ids = list(pd.unique(d["id"]))
    m = (d.groupby(["id", "category"], sort=False)["count"].sum().unstack(fill_value=0.0)
         .reindex(index=ids, columns=cats, fill_value=0.0).to_numpy(dtype=float))
    totals = m.sum(axis=1)
    if (totals == 0).any():
        steps.append(f"Dropped {int((totals == 0).sum())} sample(s) with a total count of zero.")
        ids = [s for s, t in zip(ids, totals) if t > 0]
        m = m[totals > 0]
        totals = totals[totals > 0]
    if len(cats) < 2:
        raise ValueError("composition needs at least 2 categories")
    if not ids:
        raise ValueError("composition needs at least 1 sample with a count above zero")
    first_group = d.drop_duplicates("id").set_index("id")["group"]
    g_of = np.array([first_group[s] for s in ids])
    lv = sorted(set(g_of))
    order = []
    for gl in lv:
        in_g = [i for i in range(len(ids)) if g_of[i] == gl]
        order += sorted(in_g, key=functools.partial(_id_key, ids=ids))
    ids = [ids[i] for i in order]
    g_of, m, totals = g_of[order], m[order], totals[order]
    prop = m / totals[:, None]
    if not steps:
        steps = ["No cleaning applied; input used as-is."]

    n_by = [int(np.sum(g_of == gl)) for gl in lv]
    run_test = has_group and len(lv) >= 2 and all(n >= 2 for n in n_by)
    test_label = ("none (descriptive)" if not run_test
                  else "Mann-Whitney U test" if len(lv) == 2 else "Kruskal-Wallis test")
    stat = np.full(len(cats), np.nan)
    pval = np.full(len(cats), np.nan)
    if run_test:
        for j in range(len(cats)):
            stat[j], pval[j] = test_category(prop[:, j], g_of, lv)
    padj = bh_adjust(pval) if run_test else np.full(len(cats), np.nan)

    rows = []
    for j, c in enumerate(cats):
        for gl in lv:
            v = prop[g_of == gl, j]
            rows.append({"category": c, "group": gl, "n_samples": int(len(v)),
                         "mean_proportion": float(v.mean()),
                         "sd_proportion": float(np.std(v, ddof=1)) if len(v) >= 2 else float("nan"),
                         "median_proportion": float(np.median(v)),
                         "test": test_label, "statistic": stat[j], "p_value": pval[j],
                         "p_adj": padj[j], "significant": bool(padj[j] < 0.05)})
    stats_df = pd.DataFrame(rows)

    labels = [legend_label(c, padj[j]) for j, c in enumerate(cats)]
    lay = dims(ids, g_of, lv, labels)
    fig = _plot(prop, ids, g_of, lv, cats, labels, padj, run_test, has_group, test_label, lay,
                spec, id_col, cat_col)

    long = pd.DataFrame({"sample": np.repeat([ids], len(cats), axis=0).ravel(),
                         "group": np.repeat([g_of], len(cats), axis=0).ravel(),
                         "category": np.repeat(cats, len(ids)),
                         "count": m.T.ravel(), "proportion": prop.T.ravel()})
    name = (f"{test_label} per category with Benjamini-Hochberg adjustment" if run_test
            else "descriptive (no test)")
    test_meta = {
        "name": name, "n_samples": len(ids), "n_categories": len(cats),
        "groups": dict(zip(lv, n_by)) if has_group else None,
        "p_adjust": "BH" if run_test else "none",
        "n_significant": int(np.sum(padj < 0.05)),
        "categories": {c: {"statistic": stat[j], "p_value": pval[j], "p_adj": padj[j]}
                       for j, c in enumerate(cats)},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(ids, cats, lv, n_by, has_group, run_test, test_label, padj,
                            id_col, cat_col),
        "df_used": long, "clean_steps": steps, "label": slugify(cat_col),
        "width": lay["width"], "height": lay["height"],
        "build_script": functools.partial(_script, spec, has_count, has_group),
    }


def _id_key(i, ids):
    """Natural sort key of the sample at index ``i``.

    Args:
        i: The sample index.
        ids: The sample names.

    Returns:
        tuple: The key of ``natural_key``.
    """
    return natural_key(ids[i])


def _plot(prop, ids, g_of, lv, cats, labels, padj, run_test, has_group, test_label, lay, spec,
          id_col, cat_col):
    """Draw one 100% stacked bar per sample, one panel per group.

    Args:
        prop: The proportions, samples in rows and categories in columns.
        ids: The sample names, in plot order.
        g_of: The group of each sample.
        lv: The group levels.
        cats: The categories.
        labels: The legend labels.
        padj: The adjusted P value of each category.
        run_test: Whether the categories were tested.
        has_group: Whether there is a group column.
        test_label: The test name.
        lay: The layout from ``dims``.
        spec: The figure specification.
        id_col: The sample column.
        cat_col: The category column.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    ap = spec["appearance"]
    pal = palette(len(cats), ap.get("palette"))
    counts = [int(np.sum(g_of == gl)) for gl in lv]
    fig, axes = plt.subplots(1, len(lv), figsize=(lay["width"], lay["height"]), sharey=True,
                             squeeze=False, gridspec_kw={"width_ratios": counts})
    axes = axes[0]
    for k, (ax, gl) in enumerate(zip(axes, lv)):
        rows = np.where(g_of == gl)[0]
        xpos = np.arange(len(rows))
        # The first category sits on top of the stack, as in ggplot2, so the
        # stack reads in the order of the legend.
        bottom = np.zeros(len(rows))
        for j in range(len(cats) - 1, -1, -1):
            vals = prop[rows, j]
            ax.bar(xpos, vals, bottom=bottom, width=0.85, color=pal[j], edgecolor="white",
                   linewidth=0.3)
            bottom += vals
        ax.set_xticks(xpos)
        ax.set_xticklabels([ids[i] for i in rows],
                           rotation=90 if lay["rotate"] else 0,
                           fontsize=8.5 if lay["small_x"] else 10.4)
        ax.set_xlim(-0.6, len(rows) - 0.4)
        ax.set_ylim(0, 1.01)
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, decimals=0))
        if has_group:
            ax.set_title(wrap_label(gl, lay["strip_chars"][k]), fontsize=10)
        apply_pub_style(ax, ap.get("theme"))
        if lay["small_x"]:
            ax.tick_params(axis="x", labelsize=8.5)
        if k > 0:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", left=False)
    axes[0].set_ylabel(ap.get("y_label") or "Proportion")
    n_sig = int(np.sum(padj < 0.05))
    if run_test:
        sub = (f"{test_label} per category, BH-adjusted\n{n_sig} of {len(cats)} categories with "
               f"P adj < 0.05{' (* in the legend)' if n_sig else ''}")
    elif has_group:
        sub = f"Composition of {len(ids)} samples\nno test, because a group has fewer than 2 samples"
    else:
        sub = f"Composition of {len(ids)} samples (descriptive, no test)"
    handles = [Patch(facecolor=pal[j], edgecolor="none", label=labels[j]) for j in range(len(cats))]
    axes[-1].legend(handles=handles, title=cat_col, loc="center left", bbox_to_anchor=(1.04, 0.5),
                    frameon=False, fontsize=9, title_fontsize=10, ncol=lay["leg_cols"],
                    handlelength=1.1, handleheight=1.1, labelspacing=0.45, alignment="left")
    h = lay["height"]
    top_in = 0.06
    header = []
    if ap.get("title"):
        header.append(fig.text(0.5, 1 - (top_in + 0.12) / h, ap["title"], ha="center",
                               va="center", fontsize=13, fontweight="bold"))
        top_in += 0.3
    header.append(fig.text(0.5, 1 - top_in / h, sub, ha="center", va="top", fontsize=9.5,
                           color="#4d4d4d"))
    top = 1 - (top_in + 0.17 * (sub.count("\n") + 1) - 0.02) / h
    fig.tight_layout(rect=(0, 0.27 / h, 1, top), w_pad=0.4)
    # Centre the header and the x title on the panels, as ggplot2 does.
    x0, x1 = axes[0].get_position().x0, axes[-1].get_position().x1
    for txt in header:
        txt.set_x((x0 + x1) / 2)
    fig.text((x0 + x1) / 2, 0.03 / h, ap.get("x_label") or id_col, ha="center", va="bottom",
             fontsize=13)
    return fig


def _methods(ids, cats, lv, n_by, has_group, run_test, test_label, padj, id_col, cat_col):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        ids: The samples.
        cats: The categories.
        lv: The group levels.
        n_by: The sample count of each group.
        has_group: Whether there is a group column.
        run_test: Whether the categories were tested.
        test_label: The test name.
        padj: The adjusted P values.
        id_col: The sample column.
        cat_col: The category column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    txt = (f"The counts of {len(cats)} categories of {cat_col} were summed within each of "
           f"{len(ids)} samples ({id_col}) and converted to proportions of the sample total.")
    if has_group:
        txt += (f" The samples fall into {len(lv)} groups ("
                + ", ".join(f"{g}: n = {n}" for g, n in zip(lv, n_by)) + ").")
    if run_test:
        exact = (" (two-sided, with the exact null distribution when there are fewer than 50 samples "
                 "per group and no ties)" if len(lv) == 2 else "")
        txt += (" For each category, the per-sample proportions were compared between the groups "
                f"with the {test_label}{exact}, and the P values were adjusted across the "
                f"{len(cats)} categories with the Benjamini-Hochberg method; "
                f"{int(np.sum(padj < 0.05))} categories had an adjusted P below 0.05.")
    elif has_group:
        txt += " No test was run, because a group holds fewer than 2 samples."
    else:
        txt += " The composition is described without a test."
    if has_group:
        txt += (" Proportions are compositional: they sum to one within each sample, so a rise in "
                "one category forces the others down and the per-category tests are not "
                "independent. The tests are therefore exploratory; a compositional method, for "
                "example a log-ratio transform, is needed for a confirmatory claim.")
    return (txt + f" Computed in Python {platform.python_version()} with scipy {v('scipy')} and "
            f"drawn with matplotlib {v('matplotlib')}.")


def _script(spec, has_count, has_group, in_name, fig_stub):
    """Write the standalone reproduction script.

    Args:
        spec: The figure specification.
        has_count: Whether there is a count column.
        has_group: Whether there is a group column.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script.
    """
    data = spec["data"]
    cnt = f'pd.to_numeric(raw["{data["y"]}"], errors="coerce")' if has_count else "1.0"
    grp = f'raw["{data["group"]}"].astype(str)' if has_group else '"all"'
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("{in_name}")
d = pd.DataFrame({{"id": raw["{data["id"]}"].astype(str), "category": raw["{data["x"]}"].astype(str),
                  "group": {grp}, "count": {cnt}}}).dropna()
cats = list(pd.unique(d["category"]))
tab = d.pivot_table(index="id", columns="category", values="count", aggfunc="sum", fill_value=0)[cats]
tab = tab[tab.sum(axis=1) > 0]
prop = tab.div(tab.sum(axis=1), axis=0)
g = d.drop_duplicates("id").set_index("id").loc[prop.index, "group"].to_numpy()
lv = sorted(set(g))
for c in cats if len(lv) >= 2 else []:
    v = prop[c].to_numpy()
    if len(lv) == 2:
        # The p-value rule of R's wilcox.test: exact below 50 per group without ties.
        a, b = v[g == lv[0]], v[g == lv[1]]
        exact = len(a) < 50 and len(b) < 50 and len(np.unique(v)) == len(v)
        print(c, stats.mannwhitneyu(a, b, alternative="two-sided",
                                    method="exact" if exact else "asymptotic", use_continuity=True))
    else:
        print(c, stats.kruskal(*[v[g == x] for x in lv]))
order = np.argsort(g, kind="stable")
prop = prop.iloc[order]
fig, ax = plt.subplots(figsize=(max(4.8, 2.5 + 0.3 * len(prop)), 4.8))
bottom = np.zeros(len(prop))
for c in reversed(cats):
    ax.bar(np.arange(len(prop)), prop[c].to_numpy(), bottom=bottom, width=0.85, edgecolor="white", label=c)
    bottom += prop[c].to_numpy()
ax.set_xticks(np.arange(len(prop))); ax.set_xticklabels(list(prop.index), rotation=90)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0)); ax.set_ylabel("Proportion")
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1], title="{data["x"]}", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
