"""Recipe: compare a continuous outcome across 3+ groups. Mirrors the R recipe:
one-way ANOVA / Welch / Kruskal-Wallis from the assumptions, matching post-hoc,
brackets for significant pairs only, omnibus test as a subtitle."""
import itertools

import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
import seaborn as sns
from scipy import stats
from statannotations.Annotator import Annotator

from ..clean import clean_xy
from ..theme import apply_pub_style, new_fig, palette
from ..util import cap_first, format_p, p_stars, shapiro_normal


def _resolve(method, all_normal, equal_var):
    method = (method or "auto").lower()
    if method == "auto":
        fam = "anova" if (all_normal and equal_var) else "welch" if all_normal else "kruskal"
    elif method in ("anova", "aov"):
        fam = "anova"
    elif method in ("welch", "welch_anova"):
        fam = "welch"
    elif method in ("kruskal", "kruskal_wallis", "kw"):
        fam = "kruskal"
    else:
        raise ValueError(f"unknown ANOVA method '{method}'")
    label = {"anova": "one-way ANOVA", "welch": "Welch's one-way ANOVA",
             "kruskal": "Kruskal-Wallis test"}[fam]
    posthoc = {"anova": "Tukey's HSD", "welch": "the Games-Howell test",
               "kruskal": "Dunn's test with Benjamini-Hochberg adjustment"}[fam]
    return fam, label, posthoc


def recipe_multi_group_compare(df, spec):
    x, y = spec["data"]["x"], spec["data"]["y"]
    df, steps = clean_xy(df, x, y)
    df[x] = df[x].astype(str)
    lv = sorted(df[x].unique())
    if len(lv) < 2:
        raise ValueError(f"multi_group_compare needs 2+ groups in '{x}'")

    groups = [df.loc[df[x] == g, y].to_numpy() for g in lv]
    all_normal = all(bool(shapiro_normal(g)[1]) for g in groups)
    try:
        lev_p = float(stats.levene(*groups, center="mean").pvalue)
    except Exception:
        lev_p = np.nan
    equal_var = np.isnan(lev_p) or lev_p > 0.05
    fam, label, posthoc_label = _resolve(spec["test"].get("method"), all_normal, equal_var)

    if fam == "anova":
        row = pg.anova(data=df, dv=y, between=x).iloc[0]
        omni = dict(stat=float(row.F), df1=row.ddof1, df2=row.ddof2, p=float(row.p_unc),
                    eff_name="partial eta-squared", eff=float(row.np2))
        ph = pg.pairwise_tukey(data=df, dv=y, between=x)
        pairs = [(str(r.A), str(r.B), float(r.p_tukey)) for r in ph.itertuples()]
    elif fam == "welch":
        row = pg.welch_anova(data=df, dv=y, between=x).iloc[0]
        omni = dict(stat=float(row.F), df1=row.ddof1, df2=row.ddof2, p=float(row.p_unc),
                    eff_name=None, eff=np.nan)
        ph = pg.pairwise_gameshowell(data=df, dv=y, between=x)
        pairs = [(str(r.A), str(r.B), float(r.pval)) for r in ph.itertuples()]
    else:
        row = pg.kruskal(data=df, dv=y, between=x).iloc[0]
        omni = dict(stat=float(row.H), df1=row.ddof1, df2=np.nan, p=float(row.p_unc),
                    eff_name=None, eff=np.nan)
        m = sp.posthoc_dunn(df, val_col=y, group_col=x, p_adjust="fdr_bh")
        pairs = [(str(i), str(j), float(m.loc[i, j])) for i, j in itertools.combinations(list(m.index), 2)]

    ap = spec["appearance"]
    pal = palette(len(lv), ap.get("palette"))
    geom = ap.get("geom") or "box"
    blabel = (ap.get("bracket") or {}).get("label", "p.signif")
    w, h = max(4.2, 1.1 * len(lv) + 1.0), 4.2
    fig, ax = new_fig(w, h)
    if geom == "violin":
        sns.violinplot(data=df, x=x, y=y, hue=x, order=lv, palette=pal, ax=ax,
                       legend=False, cut=0, inner=None, alpha=0.9)
    else:
        sns.boxplot(data=df, x=x, y=y, hue=x, order=lv, palette=pal, ax=ax,
                    legend=False, width=0.6, fliersize=0)
    if ap.get("show_points", True):
        sns.stripplot(data=df, x=x, y=y, hue=x, order=lv, palette=pal, ax=ax,
                      legend=False, size=3.5, alpha=0.7, edgecolor="black",
                      linewidth=0.3, jitter=0.12)

    sig = [(a, b, p) for a, b, p in pairs if p < 0.05]
    if sig:
        anns = [f"p = {p:.2g}" if blabel == "p.format" else p_stars(p) for _, _, p in sig]
        annot = Annotator(ax, [(a, b) for a, b, _ in sig], data=df, x=x, y=y, order=lv)
        annot.configure(line_width=0.9, fontsize=10)
        annot.set_custom_annotations(anns)
        annot.annotate()

    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    ax.set_title(f"{cap_first(label)}, P = {omni['p']:.2g}", fontsize=10, color="#4d4d4d")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    stats_df = pd.DataFrame([{
        "outcome": y, "group1": a, "group2": b, "post_hoc": posthoc_label,
        "p_adj": p, "p_adj_signif": p_stars(p), "displayed": p < 0.05,
    } for a, b, p in pairs])

    n_by = {g: int(len(gr)) for g, gr in zip(lv, groups)}
    test_meta = {
        "name": label,
        "omnibus": {"statistic": omni["stat"], "df1": omni["df1"], "df2": omni["df2"],
                    "p_value": omni["p"],
                    "effect_size": None if omni["eff_name"] is None else
                    {"name": omni["eff_name"], "value": omni["eff"]}},
        "post_hoc": posthoc_label,
        "n_comparisons": len(pairs), "n_significant": len(sig),
        "assumptions": {"normality": {"test": "shapiro-wilk", "all_normal": all_normal},
                        "equal_variance": {"test": "Levene", "p": lev_p}},
    }
    methods = _methods(label, posthoc_label, omni, n_by, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": steps, "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, fam, label, lv, pal, in_name, fig_stub, blabel),
    }


def _methods(label, posthoc_label, omni, n_by, spec):
    from importlib.metadata import version as v
    import platform
    y_label = spec["appearance"].get("y_label") or spec["data"]["y"]
    if label.startswith("Kruskal"):
        omni_txt = f"a Kruskal-Wallis test (H({int(omni['df1'])}) = {omni['stat']:.2f}, P = {omni['p']:.2g})"
    else:
        art = "a " if label == "one-way ANOVA" else ""
        omni_txt = f"{art}{label} (F({int(omni['df1'])}, {omni['df2']:.0f}) = {omni['stat']:.2f}, P = {omni['p']:.2g})"
    eff = f" The effect size was {omni['eff_name']} = {omni['eff']:.2f}." if omni["eff_name"] else ""
    ns = ", ".join(f"{g} (n = {n})" for g, n in n_by.items())
    return (
        f"{cap_first(y_label)} was compared across {len(n_by)} groups with {omni_txt}." + eff +
        f" Pairwise comparisons used {posthoc_label}; only comparisons reaching significance"
        " (adjusted P < 0.05) are annotated on the figure."
        f" Analyses were performed in Python {platform.python_version()} with pingouin {v('pingouin')}"
        f" and scikit-posthocs {v('scikit-posthocs')}; figures were produced with seaborn {v('seaborn')}."
        f" Group sizes: {ns}."
    )


def _script(spec, fam, label, lv, pal, in_name, fig_stub, blabel):
    x, y = spec["data"]["x"], spec["data"]["y"]
    ap = spec["appearance"]
    geom = ap.get("geom") or "box"
    if fam == "anova":
        omni_line = f'row = pg.anova(data=df, dv="{y}", between="{x}").iloc[0]; omni_p = float(row.p_unc)'
        ph_line = (f'ph = pg.pairwise_tukey(data=df, dv="{y}", between="{x}"); '
                   f'pairs = [(str(r.A), str(r.B), float(r.p_tukey)) for r in ph.itertuples()]')
    elif fam == "welch":
        omni_line = f'row = pg.welch_anova(data=df, dv="{y}", between="{x}").iloc[0]; omni_p = float(row.p_unc)'
        ph_line = (f'ph = pg.pairwise_gameshowell(data=df, dv="{y}", between="{x}"); '
                   f'pairs = [(str(r.A), str(r.B), float(r.pval)) for r in ph.itertuples()]')
    else:
        omni_line = f'row = pg.kruskal(data=df, dv="{y}", between="{x}").iloc[0]; omni_p = float(row.p_unc)'
        ph_line = (f'import itertools, scikit_posthocs as sp; '
                   f'm = sp.posthoc_dunn(df, val_col="{y}", group_col="{x}", p_adjust="fdr_bh"); '
                   f'pairs = [(str(i), str(j), float(m.loc[i, j])) for i, j in itertools.combinations(list(m.index), 2)]')
    geom_line = (f'sns.violinplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, palette={pal}, ax=ax, legend=False, cut=0, inner=None)'
                 if geom == "violin" else
                 f'sns.boxplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, palette={pal}, ax=ax, legend=False, width=0.6, fliersize=0)')
    star = 'lambda p: "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*"'
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, seaborn as sns, pingouin as pg
from statannotations.Annotator import Annotator

df = pd.read_csv("{in_name}"); df["{x}"] = df["{x}"].astype(str)
order = {lv}
{omni_line}
{ph_line}
stars = {star}
sig = [(a, b, p) for a, b, p in pairs if p < 0.05]

fig, ax = plt.subplots(figsize=(4.6, 4.2))
{geom_line}
sns.stripplot(data=df, x="{x}", y="{y}", hue="{x}", order=order, palette={pal}, ax=ax,
              legend=False, size=3.5, alpha=0.7, edgecolor="black", linewidth=0.3, jitter=0.12)
if sig:
    annot = Annotator(ax, [(a, b) for a, b, _ in sig], data=df, x="{x}", y="{y}", order=order)
    annot.configure(line_width=0.9, fontsize=10)
    annot.set_custom_annotations([stars(p) for _, _, p in sig])
    annot.annotate()
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
ax.set_title("{cap_first(label)}, P = %s" % format(omni_p, ".2g"), fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
