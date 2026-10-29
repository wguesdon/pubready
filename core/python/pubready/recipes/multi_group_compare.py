"""Recipe: compare a continuous outcome across 3+ groups. Mirrors the R recipe:
one-way ANOVA / Welch / Kruskal-Wallis from the assumptions, matching post-hoc,
brackets for significant pairs only, or a compact letter display when more than
six pairs differ, omnibus test as a subtitle."""
import functools
import itertools

import numpy as np
import pandas as pd
import pingouin as pg
import scikit_posthocs as sp
from scipy import stats

from ..clean import clean_xy
from ..theme import (add_dist_geom, apply_pub_style, category_width, geom_top, new_fig, palette,
                     plain_number_axis, script_geom_code, wrap_labels)
from ..util import cap_first, cld_letters, df_text, p_stars, shapiro_normal

MAX_BRACKETS = 6
"""More significant pairs than this switch the figure to a compact letter display."""


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
        lev_p = float(stats.levene(*groups, center="median").pvalue)
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
        # eta-squared (H) of rstatix::kruskal_effsize, as the R engine reports.
        omni = dict(stat=float(row.H), df1=row.ddof1, df2=np.nan, p=float(row.p_unc),
                    eff_name="eta-squared (H)",
                    eff=(float(row.H) - len(lv) + 1) / (len(df) - len(lv)))
        m = sp.posthoc_dunn(df, val_col=y, group_col=x, p_adjust="fdr_bh")
        pairs = [(str(i), str(j), float(m.loc[i, j])) for i, j in itertools.combinations(list(m.index), 2)]

    ap = spec["appearance"]
    pal = palette(len(lv), ap.get("palette"))
    geom = ap.get("geom") or "box"
    blabel = (ap.get("bracket") or {}).get("label", "p.signif")
    w, h = category_width(len(lv), wrap_labels(lv)), 4.2
    fig, ax = new_fig(w, h)
    add_dist_geom(ax, df, x, y, lv, pal, geom, ap.get("show_points", True))

    # Significant comparisons only. Up to MAX_BRACKETS of them get a bracket each,
    # on its own level; more would stack too many lines, so the figure shows a
    # compact letter display instead, as in the R engine.
    sig = [(a, b, p) for a, b, p in pairs if p < 0.05]
    display = "compact letter display" if len(sig) > MAX_BRACKETS else "brackets"
    # The brackets and the letters start above the drawn geom. The spacing is a
    # share of the span of the drawn axis (from 0 for a bar), so a label never
    # reaches the next line.
    tops = {g: geom_top(gr, geom) for g, gr in zip(lv, groups)}
    ylo = float(df[y].min())
    if geom == "bar":
        ylo = min(ylo, 0.0)
    top_all = max(tops.values())
    span = top_all - ylo
    if not np.isfinite(span) or span <= 0:
        span = max(abs(top_all), 1.0)
    letters = None
    if display == "brackets" and sig:
        _draw_brackets(ax, sig, lv, top_all, span, blabel)
    elif display == "compact letter display":
        letters = cld_letters(lv, [a for a, _, _ in pairs], [b for _, b, _ in pairs],
                              [p for _, _, p in pairs])
        for i, g in enumerate(lv):
            ax.text(i, tops[g] + 0.04 * span, letters[g], ha="center", va="bottom", fontsize=10)
        ax.set_ylim(top=top_all + 0.14 * span)

    ax.set_xticks(range(len(lv)))
    ax.set_xticklabels(wrap_labels(lv))
    plain_number_axis(ax.yaxis)
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    ax.set_title(f"{cap_first(label)}, P = {omni['p']:.2g}", fontsize=10, color="#4d4d4d")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    # The same columns as the R engine. The estimate is the mean of group2 minus
    # the mean of group1 for Tukey and Games-Howell; Dunn's test has none.
    sizes = df.groupby(x)[y].size()
    means = df.groupby(x)[y].mean()
    stats_df = pd.DataFrame([{
        "outcome": y, "group1": a, "group2": b, "contrast": f"{b} - {a}",
        "n1": int(sizes[a]), "n2": int(sizes[b]), "post_hoc": posthoc_label,
        "estimate": float(means[b] - means[a]) if fam != "kruskal" else np.nan,
        "p_adj": p, "p_adj_signif": p_stars(p), "displayed": p < 0.05,
        "omnibus_test": label, "omnibus_statistic": omni["stat"], "omnibus_df1": omni["df1"],
        "omnibus_df2": omni["df2"], "omnibus_p": omni["p"],
        "omnibus_effect_size_name": omni["eff_name"], "omnibus_effect_size": omni["eff"],
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
        "display": display, "letters": letters,
        "assumptions": {"normality": {"test": "shapiro-wilk", "all_normal": all_normal},
                        "equal_variance": {"test": "Levene", "p": lev_p}},
    }
    methods = _methods(label, posthoc_label, omni, n_by, spec, display, len(sig))
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": steps, "width": w, "height": h,
        "qc": {"quantity": "group values", "panels": {str(g): gr for g, gr in zip(lv, groups)}},
        "build_script": functools.partial(_script, spec, fam, label, lv, pal, letters, w,
                                          blabel=blabel),
    }


def _bracket_key(row, pos):
    """Sort key of a bracket: the narrowest pair first, then the leftmost.

    Args:
        row: The pair as (group1, group2, p).
        pos: The position of each group on the axis.

    Returns:
        tuple: The span of the pair and its left end.
    """
    i1, i2 = pos[row[0]], pos[row[1]]
    return (abs(i2 - i1), min(i1, i2))


def _draw_brackets(ax, sig, lv, top_all, span, blabel):
    """Draw one bracket per significant pair, each on its own level.

    The rule is the one of .mg_bracket_rows() in the R engine: the narrowest
    pairs lowest, from 8% of the span above the tallest geom, 13% of the span
    apart, so a label never reaches the next line.

    Args:
        ax: The axes.
        sig: The significant pairs as (group1, group2, p).
        lv: The groups, in the order of the axis.
        top_all: The top of the tallest geom.
        span: The span of the drawn axis.
        blabel: "p.signif" for stars or "p.format" for the P value.
    """
    pos = {g: i for i, g in enumerate(lv)}
    rows = sorted(sig, key=functools.partial(_bracket_key, pos=pos))
    levels = [top_all + span * (0.08 + 0.13 * k) for k in range(len(rows))]
    top = levels[-1] + 0.10 * span
    tip = 0.01 * (top - (top_all - span))
    for (a, b, p), yk in zip(rows, levels):
        x1, x2 = sorted((pos[a], pos[b]))
        ax.plot([x1, x1, x2, x2], [yk - tip, yk, yk, yk - tip], color="black", lw=0.8)
        lab = f"p = {p:.2g}" if blabel == "p.format" else p_stars(p)
        ax.text((x1 + x2) / 2, yk, lab, ha="center", va="bottom", fontsize=10)
    ax.set_ylim(top=top)


def _methods(label, posthoc_label, omni, n_by, spec, display, n_sig):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        label: The omnibus test.
        posthoc_label: The post hoc test.
        omni: The omnibus result.
        n_by: The size of each group.
        spec: The figure specification.
        display: "brackets" or "compact letter display".
        n_sig: The count of significant pairs.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_label = spec["appearance"].get("y_label") or spec["data"]["y"]
    if label.startswith("Kruskal"):
        omni_txt = f"a Kruskal-Wallis test (H({int(omni['df1'])}) = {omni['stat']:.2f}, P = {omni['p']:.2g})"
    else:
        art = "a " if label == "one-way ANOVA" else ""
        omni_txt = (f"{art}{label} (F({df_text(omni['df1'])}, {df_text(omni['df2'])}) = "
                    f"{omni['stat']:.2f}, P = {omni['p']:.2g})")
    eff = f" The effect size was {omni['eff_name']} = {omni['eff']:.2f}." if omni["eff_name"] else ""
    ns = ", ".join(f"{g} (n = {n})" for g, n in n_by.items())
    if display == "compact letter display":
        ph_txt = (f" Pairwise comparisons used {posthoc_label}. Because {n_sig} pairs differed"
                  " (adjusted P < 0.05), the figure shows a compact letter display (the"
                  " insert-and-absorb algorithm of Piepho, 2004) instead of brackets: groups that"
                  " share a letter do not differ significantly.")
    else:
        ph_txt = (f" Pairwise comparisons used {posthoc_label}; only comparisons reaching"
                  " significance (adjusted P < 0.05) are annotated on the figure with brackets.")
    return (
        f"{cap_first(y_label)} was compared across {len(n_by)} groups with {omni_txt}." + eff +
        ph_txt +
        f" Analyses were performed in Python {platform.python_version()} with pingouin {v('pingouin')}"
        f" and scikit-posthocs {v('scikit-posthocs')}; figures were produced with seaborn {v('seaborn')}."
        f" Group sizes: {ns}."
    )


def _script(spec, fam, label, lv, pal, letters, w, in_name, fig_stub, blabel="p.signif"):
    """Write the standalone script that redraws the figure.

    Args:
        spec: The figure specification.
        fam: The test family: anova, welch or kruskal.
        label: The omnibus test.
        lv: The groups, in the order of the axis.
        pal: Their colours.
        letters: The compact letter display, or None for brackets.
        w: The figure width.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.
        blabel: "p.signif" for stars or "p.format" for the P value.

    Returns:
        str: The script.
    """
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
    geom_code = script_geom_code(x, y, pal, geom, ap.get("show_points", True))
    lab_expr = 'f"p = {{p:.2g}}"' if blabel == "p.format" else "stars(p)"  # escaped for format()
    if letters is not None:
        # With more than six significant pairs the figure shows letters.
        annot = (f'letters = {letters!r}  # compact letter display (Piepho 2004)\n'
                 f'span = df["{y}"].max() - df["{y}"].min()\n'
                 'for i, g in enumerate(order):\n'
                 f'    top = df.loc[df["{x}"] == g, "{y}"].max()\n'
                 '    ax.text(i, top + 0.04 * span, letters[g], ha="center", va="bottom", fontsize=10)')
    else:
        annot = ('span = df["{y}"].max() - df["{y}"].min(); top = df["{y}"].max()\n'
                 'pos = {{g: i for i, g in enumerate(order)}}\n'
                 'sig.sort(key=bracket_key)\n'
                 'for k, (a, b, p) in enumerate(sig):\n'
                 '    yk = top + span * (0.08 + 0.13 * k); x1, x2 = sorted((pos[a], pos[b]))\n'
                 '    ax.plot([x1, x1, x2, x2], [yk - 0.01 * span, yk, yk, yk - 0.01 * span], color="black", lw=0.8)\n'
                 f'    ax.text((x1 + x2) / 2, yk, {lab_expr}, ha="center", va="bottom", fontsize=10)').format(y=y)
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd, seaborn as sns, pingouin as pg


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*"


def bracket_key(row):
    i1, i2 = order.index(row[0]), order.index(row[1])
    return (abs(i2 - i1), min(i1, i2))


df = pd.read_csv("{in_name}"); df["{x}"] = df["{x}"].astype(str)
order = {lv}
{omni_line}
{ph_line}
sig = [(a, b, p) for a, b, p in pairs if p < 0.05]

fig, ax = plt.subplots(figsize=({w:.1f}, 4.2))
{geom_code}
{annot}
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
ax.set_title("{cap_first(label)}, P = %s" % format(omni_p, ".2g"), fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
