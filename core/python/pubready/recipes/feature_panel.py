"""Recipe: small multiples, one panel per feature, with a test in each panel.

Mirrors the R recipe. Each feature (a gene, a protein, a cytokine) gets its own
panel and a test across the groups, and the p values are adjusted across the
features with the Benjamini-Hochberg method. Input is long (--x group --y value
--facet feature) or wide (--x group, and every other numeric column is a
feature). One test family covers the whole panel: two groups take the Welch
t-test or the Mann-Whitney U test, three or more the one-way ANOVA or the
Kruskal-Wallis test. --test auto picks the parametric family only when the
Shapiro-Wilk test accepts normality for every group of every feature.
"""
import functools
import math
import os
import textwrap

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import (cap_first, is_zero_variance, mann_whitney_r, rank_biserial_unpaired,
                    rank_floor_text, rank_test_cannot_reach, rank_test_min_p, shapiro_normal)

MAX_FEATURES = 36
MAX_QC_PANELS = 48
MIN_N = 2
LABELS = {"welch_t": "Welch t-test", "wilcoxon": "Mann-Whitney U test",
          "anova": "one-way ANOVA", "kruskal": "Kruskal-Wallis test"}
EFFECTS = {"welch_t": "Cohen's d", "wilcoxon": "rank-biserial r",
           "anova": "eta squared", "kruskal": "eta-squared (H)"}
ALIASES = {"mann_whitney": "wilcoxon", "mwu": "wilcoxon", "wilcox": "wilcoxon",
           "kw": "kruskal", "kruskal_wallis": "kruskal"}


def signif_stars(p):
    """Significance stars with the cutpoints of rstatix::add_significance.

    Args:
        p: The adjusted p value, or NaN.

    Returns:
        str: The stars, "ns", or "NA" for a missing value.
    """
    if p is None or np.isnan(p):
        return "NA"
    for cut, star in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*")):
        if p < cut:
            return star
    return "ns"


def bh_adjust(p):
    """Benjamini-Hochberg adjustment with the rules of R's p.adjust.

    Missing values stay missing and do not count as comparisons.

    Args:
        p: The p values.

    Returns:
        numpy.ndarray: The adjusted p values.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(-pv, kind="stable")
    ranks = np.arange(n, 0, -1)
    adj = np.minimum.accumulate(pv[order] * n / ranks)
    res = np.empty(n)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def group_label(v):
    """Text of one group value, as R's as.character writes it.

    Args:
        v: The value.

    Returns:
        str: The label; a whole number has no decimal part.
    """
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        return str(int(v))
    return str(v)


def group_levels(g):
    """The group levels: numbers in numeric order, text in code-point order.

    Args:
        g: The group values, as a pandas Series without missing values.

    Returns:
        list: The level labels.
    """
    if pd.api.types.is_numeric_dtype(g):
        return [group_label(v) for v in sorted(pd.unique(g))]
    return sorted(pd.unique(g.astype(str)))


def long_table(df, spec):
    """Bring a long or a wide table to one row per value.

    Args:
        df: The input table.
        spec: The figure specification.

    Returns:
        dict: The long data, the features, the groups, the steps and the layout.

    Raises:
        ValueError: If a column is missing.
    """
    dat = spec["data"]
    x, fcol, y = dat.get("x"), dat.get("facet"), dat.get("y")
    if not x:
        raise ValueError("feature_panel needs --x (the group column)")
    if x not in df.columns:
        raise ValueError(f"column '{x}' not found in data")
    if fcol:
        if not y:
            raise ValueError("feature_panel in long form needs --y (the value column)")
        for c in (fcol, y):
            if c not in df.columns:
                raise ValueError(f"column '{c}' not found in data")
        d = pd.DataFrame({"feature": df[fcol].astype(str).where(df[fcol].notna()),
                          "group": df[x], "value": pd.to_numeric(df[y], errors="coerce")})
        layout = "long"
    else:
        skip = {x, dat.get("id"), dat.get("label")}
        cols = [c for c in df.columns
                if pd.api.types.is_numeric_dtype(df[c]) and c not in skip]
        if not cols:
            raise ValueError("feature_panel found no numeric feature column; pass --facet for long data")
        d = pd.concat([pd.DataFrame({"feature": c, "group": df[x],
                                     "value": pd.to_numeric(df[c], errors="coerce")})
                       for c in cols], ignore_index=True)
        layout = "wide"
    n0 = len(d)
    keep = (d["feature"].notna() & (d["feature"].astype(str).str.len() > 0)
            & d["group"].notna() & np.isfinite(d["value"]))
    d = d[keep].reset_index(drop=True)
    glv = group_levels(d["group"])
    d["group"] = [group_label(v) for v in d["group"]]
    feats = list(pd.unique(d["feature"]))
    steps = [f"Read {layout} data: {len(feats)} feature(s), {len(glv)} group(s)."]
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} value(s) with a missing feature, group or value "
                     f"({n0} -> {len(d)}).")
    return {"d": d, "features": feats, "groups": glv, "steps": steps, "layout": layout}


def resolve_test(method, n_groups, all_normal):
    """Pick the test of the whole panel.

    Args:
        method: The --test value.
        n_groups: The count of groups.
        all_normal: Whether every feature-by-group sample passed Shapiro-Wilk.

    Returns:
        str: One of welch_t, wilcoxon, anova or kruskal.

    Raises:
        ValueError: If the test is unknown or does not fit the count of groups.
    """
    method = ALIASES.get((method or "auto").lower(), (method or "auto").lower())
    if method == "auto":
        if all_normal:
            return "welch_t" if n_groups == 2 else "anova"
        return "wilcoxon" if n_groups == 2 else "kruskal"
    if method not in LABELS:
        raise ValueError(f"feature_panel: unknown test '{method}'; use auto, welch_t, wilcoxon, "
                         "anova or kruskal")
    if method in ("welch_t", "wilcoxon") and n_groups != 2:
        raise ValueError(f"feature_panel: --test {method} needs exactly 2 groups; found "
                         f"{n_groups}. Use anova or kruskal.")
    return method


def cohens_d_pooled(a, b):
    """Cohen's d of two groups with the pooled standard deviation.

    Args:
        a: Values of the first group.
        b: Values of the second group.

    Returns:
        float: The effect size, negative when the first group is lower.
    """
    na, nb = len(a), len(b)
    sp = math.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return float((np.mean(a) - np.mean(b)) / sp)


def run_test(groups, test):
    """Run the test of one feature.

    Args:
        groups: The values of each group, in level order.
        test: The test key.

    Returns:
        list: The statistic, df1, df2, p and the effect size.
    """
    if test == "welch_t":
        a, b = groups
        va, vb = np.var(a, ddof=1) / len(a), np.var(b, ddof=1) / len(b)
        if va + vb == 0:
            return [np.nan] * 5
        res = stats.ttest_ind(a, b, equal_var=False)
        df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        return [float(res.statistic), float(df), np.nan, float(res.pvalue), cohens_d_pooled(a, b)]
    if test == "wilcoxon":
        a, b = groups
        u, p = mann_whitney_r(a, b)
        return [u, np.nan, np.nan, p, rank_biserial_unpaired(u, len(a), len(b))]
    v = np.concatenate(groups)
    if np.all(v == v[0]):
        return [np.nan] * 5
    if test == "anova":
        res = stats.f_oneway(*groups)
        ssb = sum(len(g) * (np.mean(g) - np.mean(v)) ** 2 for g in groups)
        return [float(res.statistic), float(len(groups) - 1), float(len(v) - len(groups)),
                float(res.pvalue), float(ssb / np.sum((v - np.mean(v)) ** 2))]
    res = stats.kruskal(*groups)
    h = float(res.statistic)
    # eta-squared (H) of rstatix::kruskal_effsize, as multi_group_compare reports.
    k = len(groups)
    return [h, float(k - 1), np.nan, float(res.pvalue), (h - k + 1) / (len(v) - k)]


def test_feature(d, feature, glv, test):
    """Test one feature; a group with fewer than MIN_N values gives no test.

    Args:
        d: The long data.
        feature: The feature name.
        glv: The group levels.
        test: The test key.

    Returns:
        list: The statistic, df1, df2, p and the effect size.
    """
    sub = d[d["feature"] == feature]
    groups = [sub.loc[sub["group"] == g, "value"].to_numpy(dtype=float) for g in glv]
    if any(len(g) < MIN_N for g in groups):
        return [np.nan] * 5
    try:
        return run_test(groups, test)
    except Exception:
        return [np.nan] * 5


def wrap_label(s, width=12):
    """Wrap a long group label over several lines, as scales::label_wrap does.

    Args:
        s: The label.
        width: The line width in characters.

    Returns:
        str: The wrapped label.
    """
    return textwrap.fill(str(s), width=width - 1, break_long_words=False, break_on_hyphens=False)


def title_name(s):
    """Wrap a feature name for a panel title, as the R engine does.

    The name is wrapped at 24 characters (lines of 23 or fewer, as strwrap
    gives), and a line still longer than 26 characters is cut and ends in an
    ellipsis.

    Args:
        s: The feature name.

    Returns:
        str: The title text.
    """
    lines = textwrap.wrap(str(s), width=23, break_long_words=False, break_on_hyphens=False)
    return "\n".join(ln if len(ln) <= 26 else ln[:25] + "…" for ln in lines)


def n_range(values):
    """Describe the group sizes over the features.

    Args:
        values: The size of one group in each feature.

    Returns:
        str: One size, or the smallest and the largest.
    """
    lo, hi = int(np.min(values)), int(np.max(values))
    return str(lo) if lo == hi else f"{lo} to {hi}"


def recipe_feature_panel(df, spec):
    """Draw one panel per feature with the BH-adjusted test of each.

    Args:
        df: The long or the wide table.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, there are too many features or
            fewer than two groups.
    """
    lg = long_table(df, spec)
    d, feats, glv = lg["d"], lg["features"], lg["groups"]
    if len(feats) > MAX_FEATURES:
        raise ValueError(f"feature_panel draws at most {MAX_FEATURES} features and found "
                         f"{len(feats)}; use the heatmap recipe for a larger panel")
    if len(glv) < 2:
        raise ValueError("feature_panel needs at least 2 groups")

    # A constant feature (every value the same) is not tested: it stays out of
    # the normality check, of the tests and of the BH count.
    const = np.array([is_zero_variance(d.loc[d["feature"] == f, "value"]) for f in feats], dtype=bool)
    pairs = [(f, g) for f, c in zip(feats, const) if not c for g in glv]
    keys = [f"{f}: {g}" for f, g in pairs]
    samples = {k: d.loc[(d["feature"] == f) & (d["group"] == g), "value"].to_numpy(dtype=float)
               for k, (f, g) in zip(keys, pairs)}
    sw = [bool(shapiro_normal(v)[1]) for v in samples.values()]
    all_normal = all(sw)
    method = spec["test"].get("method") or "auto"
    forced = method.lower() != "auto"
    # When the normality check fails but the rank test cannot reach P < 0.05 at
    # the group sizes of the panel, auto keeps the parametric family.
    tested = [f for f, c in zip(feats, const) if not c]
    sizes = [max([int(((d["feature"] == f) & (d["group"] == g)).sum()) for f in tested], default=0)
             for g in glv]
    rank_min_p = rank_test_min_p(sizes, "mann_whitney" if len(glv) == 2 else "kruskal")
    kept = not forced and not all_normal and rank_test_cannot_reach(rank_min_p)
    test = resolve_test(method, len(glv), all_normal or kept)

    res = np.array([test_feature(d, f, glv, test) for f in feats], dtype=float).reshape(len(feats), 5)
    res[const, :] = np.nan
    p = res[:, 3]
    q = bh_adjust(p)
    stats_df = pd.DataFrame({
        "feature": feats, "test": LABELS[test], "statistic": res[:, 0], "df1": res[:, 1],
        "df2": res[:, 2], "p_value": p, "p_adj": q, "p_adj_signif": [signif_stars(v) for v in q],
        "effect_size": res[:, 4], "effect_size_name": EFFECTS[test],
        "note": ["not tested: constant" if c else "" for c in const]})
    for g in glv:
        stats_df[f"n_{g}"] = [int(((d["feature"] == f) & (d["group"] == g)).sum()) for f in feats]
    n_untested = int((np.isnan(p) & ~const).sum())
    n_tested = int((~np.isnan(p)).sum())
    n_sig = int(np.nansum(q < 0.05))

    ap = spec["appearance"]
    ncol = math.ceil(math.sqrt(len(feats)))
    nrow = math.ceil(len(feats) / ncol)
    pal = palette(len(glv), ap.get("palette"))
    panel_w = max(2.2, 0.55 * len(glv) + 1.0)
    w, h = 0.9 + panel_w * ncol, 1.2 + 2.45 * nrow
    fig, _ = new_fig(w, h)
    fig.clf()
    axes = fig.subplots(nrow, ncol, squeeze=False)
    y_lab = ap.get("y_label") or (spec["data"]["y"] if lg["layout"] == "long" else "Value")
    ticks = [wrap_label(g) for g in glv]
    for i, f in enumerate(feats):
        ax = axes[i // ncol][i % ncol]
        sub = d[d["feature"] == f]
        vals = [sub.loc[sub["group"] == g, "value"].to_numpy(dtype=float) for g in glv]
        present = [k for k, v in enumerate(vals) if len(v)]
        if present:
            bp = ax.boxplot([vals[k] for k in present], positions=present, widths=0.6,
                            patch_artist=True, showfliers=False, manage_ticks=False,
                            medianprops={"color": "black", "lw": 1.0},
                            whiskerprops={"lw": 0.8}, showcaps=False,
                            boxprops={"lw": 0.8})
            for patch, k in zip(bp["boxes"], present):
                patch.set_facecolor(pal[k])
                patch.set_alpha(0.85)
        for k, v in enumerate(vals):
            jit = np.random.RandomState(k + 1).uniform(-0.12, 0.12, size=len(v))
            ax.scatter(np.full(len(v), k) + jit, v, s=9, color=pal[k], alpha=0.8,
                       edgecolors="black", linewidths=0.25, zorder=3)
        qv = q[i]
        qtxt = ("not tested: constant" if const[i] else
                f"q = {'NA' if np.isnan(qv) else format(qv, '.2g')}")
        ax.set_title(f"{title_name(f)}\n{qtxt}", fontsize=10, linespacing=1.15)
        ax.set_xticks(range(len(glv)))
        ax.set_xticklabels(ticks, fontsize=9, linespacing=0.95)
        ax.set_xlim(-0.6, len(glv) - 0.4)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
        apply_pub_style(ax, ap.get("theme"))
        ax.tick_params(labelsize=9)
        # As facet_wrap does, the group names are written under the lowest
        # panel of each column only.
        if i + ncol < len(feats):
            ax.tick_params(axis="x", labelbottom=False)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for j in range(len(feats), nrow * ncol):
        axes[j // ncol][j % ncol].set_visible(False)
    fig.supxlabel(ap.get("x_label") or spec["data"]["x"], fontsize=13)
    fig.supylabel(y_lab, fontsize=13)
    sub_txt = (f"{LABELS[test]} per feature; BH-adjusted q across {n_tested} features; "
               f"{n_sig} with q < 0.05")
    fig.suptitle(sub_txt, fontsize=10, color="#4d4d4d")
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.12 / h))
    if ap.get("title"):
        fig.text(0.5, 1.0, ap["title"], ha="center", va="bottom", fontweight="bold", fontsize=13)

    qc = None
    if len(keys) <= MAX_QC_PANELS:
        qc = {"quantity": "group values of each feature", "panels": samples}
    test_meta = {
        "name": LABELS[test], "test": test, "p_adjust": "BH", "n_features": len(feats),
        "n_groups": len(glv), "groups": glv, "n_significant": n_sig, "n_untested": n_untested,
        "n_constant": int(const.sum()), "n_tested": n_tested, "forced": forced,
        "kept_parametric": kept, "rank_test_attainable_p": rank_min_p,
        "assumptions": {"normality": {"test": "shapiro-wilk", "all_normal": all_normal,
                                      "n_checked": len(keys), "n_non_normal": int(len(sw) - sum(sw))}},
    }
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(test_meta, stats_df, glv, y_lab, spec["data"]["x"]),
        "df_used": d[["feature", "group", "value"]], "clean_steps": lg["steps"], "qc": qc,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]), "width": w, "height": h,
        "build_script": functools.partial(_script, spec, lg["layout"], test, ncol),
    }


def _methods(tm, stats_df, glv, y_lab, x):
    """Write the methods paragraph, as the R engine words it.

    Args:
        tm: The test metadata.
        stats_df: The per-feature table.
        glv: The group levels.
        y_lab: The name of the measured value.
        x: The group column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    nc = tm["assumptions"]["normality"]["n_checked"]
    if tm["forced"]:
        floor_txt = (rank_floor_text(tm["rank_test_attainable_p"], False)
                     if tm["test"] in ("wilcoxon", "kruskal") else "")
        choice = (f"The {tm['name']} was set for every feature by the user. "
                  + (floor_txt + " " if floor_txt else ""))
    elif tm["kept_parametric"]:
        choice = ("One test family was chosen for the whole panel: the Shapiro-Wilk test was run on "
                  f"each of the {nc} feature-by-group samples, and "
                  f"{tm['assumptions']['normality']['n_non_normal']} of them departed from normality "
                  "(P <= 0.05, or fewer than 3 values). "
                  f"{rank_floor_text(tm['rank_test_attainable_p'], True)} The {tm['name']} was used "
                  "for every feature. ")
    elif tm["assumptions"]["normality"]["all_normal"]:
        choice = ("One test family was chosen for the whole panel: the Shapiro-Wilk test was run on "
                  f"each of the {nc} feature-by-group samples, and because none departed from "
                  f"normality (P > 0.05) the parametric {tm['name']} was used for every feature. ")
    else:
        choice = ("One test family was chosen for the whole panel: the Shapiro-Wilk test was run on "
                  f"each of the {nc} feature-by-group samples, and because "
                  f"{tm['assumptions']['normality']['n_non_normal']} of them departed from normality "
                  f"(P <= 0.05, or fewer than 3 values) the rank-based {tm['name']} was used for "
                  "every feature. ")
    nt = ", ".join(f"{g} (n = {n_range(stats_df[f'n_{g}'])})" for g in glv)
    untested = (f"{tm['n_untested']} feature(s) had a group with fewer than 2 values and were not "
                "tested. " if tm["n_untested"] else "")
    constant = (f"{tm['n_constant']} feature(s) were constant (all values equal) and were not "
                "tested. " if tm["n_constant"] else "")
    return (
        f"{cap_first(y_lab)} was compared between the levels of {x} ({', '.join(glv)}) separately "
        f"for each of {tm['n_features']} features. {choice}P values were adjusted across the "
        f"{tm['n_tested']} features with the Benjamini-Hochberg method; {tm['n_significant']} "
        f"features had an adjusted P (q) below 0.05. The effect size is {EFFECTS[tm['test']]}. "
        f"{constant}{untested}Group sizes per feature: {nt}. Drawn in Python {platform.python_version()} "
        f"with matplotlib {v('matplotlib')} and scipy {v('scipy')}."
    )


def _script(spec, layout, test, ncol, in_name, fig_stub):
    """Write the standalone script that redraws the panel.

    Args:
        spec: The figure specification.
        layout: "long" or "wide".
        test: The test key.
        ncol: The count of panel columns.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    dat = spec["data"]
    x = dat["x"]
    if layout == "long":
        read = (f'd = pd.DataFrame({{"feature": raw["{dat["facet"]}"].astype(str), '
                f'"group": raw["{x}"].astype(str), "value": pd.to_numeric(raw["{dat["y"]}"], '
                'errors="coerce")})')
    else:
        read = (f'cols = [c for c in raw.select_dtypes("number").columns if c != "{x}"]\n'
                f'd = pd.concat([pd.DataFrame({{"feature": c, "group": raw["{x}"].astype(str), '
                '"value": raw[c]}) for c in cols], ignore_index=True)')
    call = {"welch_t": "stats.ttest_ind(*g, equal_var=False).pvalue",
            "wilcoxon": "stats.mannwhitneyu(*g).pvalue",
            "anova": "stats.f_oneway(*g).pvalue", "kruskal": "stats.kruskal(*g).pvalue"}[test]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

raw = pd.read_csv("{in_name}")
{read}
d = d[np.isfinite(d["value"])]
feats = list(pd.unique(d["feature"]))
levels = sorted(pd.unique(d["group"]))
p = []
for f in feats:
    g = [d.loc[(d["feature"] == f) & (d["group"] == lv), "value"].to_numpy() for lv in levels]
    p.append({call})
q = multipletests(p, method="fdr_bh")[1]
nrow = math.ceil(len(feats) / {ncol})
fig, axes = plt.subplots(nrow, {ncol}, figsize=(3 * {ncol}, 2.6 * nrow), squeeze=False)
for i, f in enumerate(feats):
    ax = axes[i // {ncol}][i % {ncol}]
    g = [d.loc[(d["feature"] == f) & (d["group"] == lv), "value"].to_numpy() for lv in levels]
    ax.boxplot(g, showfliers=False)
    ax.set_xticks(range(1, len(levels) + 1)); ax.set_xticklabels(levels, fontsize=8)
    ax.set_title(f"{{f}}\\nq = {{q[i]:.2g}}", fontsize=10)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
