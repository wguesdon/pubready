"""Recipe: compare two groups on one continuous outcome. Mirrors the R recipe:
picks t-test or Mann-Whitney from the assumption checks, box/violin with points
and a significance bracket on top."""
import numpy as np
import pandas as pd
import pingouin as pg
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator

from ..clean import clean_xy
from ..theme import add_dist_geom, apply_pub_style, new_fig, palette
from ..util import (cap_first, mann_whitney_r, p_stars, rank_biserial_paired,
                    rank_biserial_unpaired, shapiro_normal, wilcoxon_signed_rank_r)


def _resolve(method, all_normal, equal_var, paired):
    method = (method or "auto").lower()
    if method == "auto":
        family = "t" if all_normal else "wilcoxon"
        var_equal = equal_var if family == "t" else None
    elif method in ("t", "t_test", "student_t", "welch_t"):
        family = "t"
        var_equal = True if method == "student_t" else (False if method == "welch_t" else equal_var)
    elif method in ("wilcoxon", "mann_whitney", "wilcox", "mwu"):
        family, var_equal = "wilcoxon", None
    else:
        raise ValueError(f"unknown test method '{method}'")
    if family == "t":
        label = ("paired t-test" if paired else
                 "Student's two-sample t-test" if var_equal else "Welch two-sample t-test")
    else:
        label = "Wilcoxon signed-rank test" if paired else "Mann-Whitney U test"
    return dict(family=family, var_equal=var_equal, label=label, paired=paired)


def recipe_two_group_compare(df, spec):
    x, y = spec["data"]["x"], spec["data"]["y"]
    df, steps = clean_xy(df, x, y)
    df[x] = df[x].astype(str)
    lv = sorted(df[x].unique())
    if len(lv) != 2:
        raise ValueError(f"two_group_compare needs exactly 2 groups in '{x}'; found {len(lv)}")
    paired = bool(spec["test"].get("paired"))

    a = df.loc[df[x] == lv[0], y].to_numpy()
    b = df.loc[df[x] == lv[1], y].to_numpy()
    p1, n1 = shapiro_normal(a); p2, n2 = shapiro_normal(b)
    all_normal = bool(n1) and bool(n2)
    try:
        var_p = float(stats.levene(a, b, center="median").pvalue)
    except Exception:
        var_p = np.nan
    equal_var = np.isnan(var_p) or var_p > 0.05

    r = _resolve(spec["test"].get("method"), all_normal, equal_var, paired)

    if r["family"] == "t":
        res = pg.ttest(a, b, paired=paired, correction=not r["var_equal"])
        statistic = float(res["T"].iloc[0]); pval = float(res["p_val"].iloc[0])
        dof = float(res["dof"].iloc[0]); eff = float(res["cohen_d"].iloc[0])
        eff_name = "Cohen's d"
    else:
        if paired:
            d = a - b
            statistic, pval = wilcoxon_signed_rank_r(d)
            eff = rank_biserial_paired(d)
        else:
            statistic, pval = mann_whitney_r(a, b)
            eff = rank_biserial_unpaired(statistic, len(a), len(b))
        dof = np.nan; eff_name = "rank-biserial r"

    ap = spec["appearance"]
    pal = palette(2, ap.get("palette"))
    geom = ap.get("geom") or "box"
    w, h = 3.8, 4.0
    fig, ax = new_fig(w, h)
    add_dist_geom(ax, df, x, y, lv, pal, geom, ap.get("show_points", True))

    blabel = (ap.get("bracket") or {}).get("label", "p.signif")
    ann = f"p = {pval:.2g}" if blabel == "p.format" else p_stars(pval)
    annot = Annotator(ax, [(lv[0], lv[1])], data=df, x=x, y=y, order=lv)
    annot.configure(line_width=1.0, fontsize=12)
    annot.set_custom_annotations([ann])
    annot.annotate()

    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    yl = ap.get("y_limits") or [None, None]
    if yl[0] is not None or yl[1] is not None:
        ax.set_ylim(yl[0], yl[1])
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    stats_df = pd.DataFrame([{
        "outcome": y, "group1": lv[0], "group2": lv[1], "n1": len(a), "n2": len(b),
        "test": r["label"], "statistic": statistic, "df": dof, "p_value": pval,
        "p_adjust_method": spec["test"].get("p_adjust", "none"),
        "effect_size_name": eff_name, "effect_size": eff,
        "shapiro_p_group1": p1, "shapiro_p_group2": p2, "var_equal_p": var_p,
    }])

    test_meta = {
        "name": r["label"], "statistic": statistic, "p_value": pval,
        "adjustment": spec["test"].get("p_adjust", "none"),
        "effect_size": {"name": eff_name, "value": eff},
        "assumptions": {
            "normality": {"test": "shapiro-wilk", "p_group1": p1, "p_group2": p2},
            "equal_variance": {"test": "Levene", "p": var_p},
        },
    }

    methods = _methods(r, stats_df, eff_name, eff, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": methods, "df_used": df, "clean_steps": steps,
        "width": w, "height": h,
        "qc": {"quantity": "group values", "panels": {lv[0]: a, lv[1]: b}},
        "build_script": lambda in_name, fig_stub: _script(spec, r, lv, pal, in_name, fig_stub, ann),
    }


def _methods(r, stats_df, eff_name, eff_value, spec):
    from importlib.metadata import version as v
    y_label = spec["appearance"].get("y_label") or spec["data"]["y"]
    var_txt = ""
    if r["family"] == "t":
        var_txt = (" Equal variances were assumed." if r["var_equal"]
                   else " Equal variances were not assumed, and the Welch correction was applied.")
    else:
        var_txt = " A non-parametric test was used because the normality assumption was not met."
    n1, n2 = int(stats_df["n1"].iloc[0]), int(stats_df["n2"].iloc[0])
    # The value is given, not only the name, so the paragraph can be read on
    # its own. The R engine writes the same sentence.
    eff_txt = (eff_name if eff_value is None or not np.isfinite(eff_value)
               else f"{eff_name} = {eff_value:.2f}")
    return (
        f"{cap_first(y_label)} was compared between the two groups with a {r['label']}."
        " Normality was assessed with the Shapiro-Wilk test." + var_txt +
        f" Effect size is reported as {eff_txt}."
        " Significance was set at P < 0.05."
        f" Analyses were performed in Python {_pyver()} with pingouin {v('pingouin')};"
        f" figures were produced with seaborn {v('seaborn')} and statannotations {v('statannotations')}."
        f" Group sizes were n = {n1} and n = {n2}."
    )


def _pyver():
    import platform
    return platform.python_version()


def _script(spec, r, lv, pal, in_name, fig_stub, ann):
    x, y = spec["data"]["x"], spec["data"]["y"]
    ap = spec["appearance"]
    geom = ap.get("geom") or "box"
    paired = "True" if r["paired"] else "False"
    if r["family"] == "t":
        test_line = (f'res = pg.ttest(a, b, paired={paired}, '
                     f'correction={not r["var_equal"]}); '
                     f'p = float(res["p_val"].iloc[0])')
    elif r["paired"]:
        # R drops the zero differences and takes the exact route only when no tie
        # and no zero remain, so the same rule is written out here.
        test_line = (
            'd = a - b; nz = d[d != 0]\n'
            'exact = len(nz) < 50 and len(nz) == len(d) and len(set(abs(nz))) == len(nz)\n'
            'res = (ss.wilcoxon(nz, alternative="two-sided", method="exact") if exact\n'
            '       else ss.wilcoxon(nz, alternative="two-sided", method="approx", correction=True))\n'
            'p = float(res.pvalue)')
    else:
        test_line = (
            'exact = len(a) < 50 and len(b) < 50 and len(set(list(a) + list(b))) == len(a) + len(b)\n'
            'res = (ss.mannwhitneyu(a, b, alternative="two-sided", method="exact") if exact\n'
            '       else ss.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic",\n'
            '                            use_continuity=True))\n'
            'p = float(res.pvalue)')
    geom_line = (f'sns.violinplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, hue_order=order, '
                 f'palette={pal}, ax=ax, legend=False, cut=0, inner=None)'
                 if geom == "violin" else
                 f'sns.boxplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, hue_order=order, '
                 f'palette={pal}, ax=ax, legend=False, width=0.6, fliersize=0)')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg
from scipy import stats as ss
from statannotations.Annotator import Annotator

df = pd.read_csv("{in_name}")
df["{x}"] = df["{x}"].astype(str)
order = {lv}
a = df.loc[df["{x}"] == order[0], "{y}"].to_numpy()
b = df.loc[df["{x}"] == order[1], "{y}"].to_numpy()
{test_line}

fig, ax = plt.subplots(figsize=(3.8, 4.0))
{geom_line}
sns.stripplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, hue_order=order, palette={pal},
              ax=ax, legend=False, size=4, alpha=0.75, edgecolor="black", linewidth=0.3, jitter=0.12)
annot = Annotator(ax, [(order[0], order[1])], data=df, x="{x}", y="{y}", order=order)
annot.configure(line_width=1.0, fontsize=12)
annot.set_custom_annotations(["{ann}"])
annot.annotate()
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
