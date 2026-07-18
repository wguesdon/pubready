"""Recipe: two-way / three-way factorial ANOVA (parametric, statsmodels Type II).
The Python engine has no standard aligned rank transform, so unlike the R engine
this path is parametric only. Grouped box, faceted by the third factor, with each
effect's p-value shown on top."""
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf

from ..theme import apply_pub_style, palette
from ..util import cap_first


def _plab(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "P = NA"
    return "P < 0.001" if p < 0.001 else f"P = {p:.2g}"


def recipe_factorial_anova(df, spec):
    y = spec["data"]["y"]
    f1, f2, f3 = spec["data"]["x"], spec["data"]["fill"], spec["data"].get("facet")
    if not f2:
        raise ValueError("factorial_anova needs at least two factors: pass --x and --fill")
    factors = [f1, f2] + ([f3] if f3 else [])
    for f in factors + [y]:
        if f not in df.columns:
            raise ValueError(f"column '{f}' not found")
    df = df.copy()
    for f in factors:
        df[f] = df[f].astype(str)
    df[y] = pd.to_numeric(df[y], errors="coerce")
    n0 = len(df)
    df = df.dropna(subset=[y] + factors).reset_index(drop=True)
    clean_steps = (f"Dropped {n0 - len(df)} row(s) with missing values ({n0} -> {len(df)} rows)."
                   if len(df) < n0 else "No cleaning applied; input used as-is.")

    ren = {y: "_y_"}
    ren.update({f: f"_f{i}_" for i, f in enumerate(factors)})
    d2 = df.rename(columns=ren)
    rhs = "*".join(f"C(_f{i}_)" for i in range(len(factors)))
    model = smf.ols(f"_y_ ~ {rhs}", data=d2).fit()
    aov = sm.stats.anova_lm(model, typ=2)
    total_ss = float(aov["sum_sq"].sum())

    def clean_term(t):
        for i, f in enumerate(factors):
            t = t.replace(f"C(_f{i}_)", f)
        return t

    rows = []
    for term, r in aov.iterrows():
        if term == "Residual":
            continue
        rows.append({"effect": clean_term(term), "df1": float(r["df"]),
                     "df2": float(aov.loc["Residual", "df"]), "statistic": float(r["F"]),
                     "p_value": float(r["PR(>F)"]), "effect_size": float(r["sum_sq"]) / total_ss,
                     "significant": float(r["PR(>F)"]) < 0.05})
    eff = pd.DataFrame(rows)

    sig = eff[eff["significant"]]
    subtitle = ("\n".join(f"{r.effect.replace(':', ' × ')}: {_plab(r.p_value)}"
                          for r in sig.itertuples())
                if len(sig) else "No significant effects (P ≥ 0.05)")

    ap = spec["appearance"]
    f1_lv = sorted(df[f1].unique()); f2_lv = sorted(df[f2].unique())
    pal = palette(len(f2_lv), ap.get("palette"))
    facet_lv = sorted(df[f3].unique()) if f3 else [None]
    w = (6.8 if f3 else 5.2); h = 4.4
    fig, axes = plt.subplots(1, len(facet_lv), figsize=(w, h), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, fl in zip(axes, facet_lv):
        sub = df if fl is None else df[df[f3] == fl]
        sns.boxplot(data=sub, x=f1, y=y, hue=f2, order=f1_lv, hue_order=f2_lv,
                    palette=pal, ax=ax, width=0.7, fliersize=0, legend=False)
        if ap.get("show_points", True):
            sns.stripplot(data=sub, x=f1, y=y, hue=f2, order=f1_lv, hue_order=f2_lv,
                          dodge=True, palette=pal, ax=ax, size=3, alpha=0.6,
                          edgecolor="black", linewidth=0.3, legend=False)
        if fl is not None:
            ax.set_title(str(fl), fontsize=11)
        ax.set_xlabel(ap.get("x_label") or f1)
        ax.set_ylabel(ap.get("y_label") or y if ax is axes[0] else "")
        apply_pub_style(ax)
    handles = [mpatches.Patch(color=pal[i], label=lvl) for i, lvl in enumerate(f2_lv)]
    axes[-1].legend(handles=handles, title=f2, loc="center left",
                    bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    fig.suptitle(subtitle, fontsize=9.5, color="#4d4d4d")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    test_meta = {
        "name": f"{len(factors)}-way ANOVA", "type": "type_II_anova",
        "factors": factors,
        "effects": [{"effect": r.effect, "df1": r.df1, "df2": r.df2,
                     "statistic": r.statistic, "p_value": r.p_value,
                     "effect_size": r.effect_size} for r in eff.itertuples()],
    }
    methods = _methods(eff, factors, spec)
    return {
        "fig": fig, "stats": eff, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, factors, f1_lv, f2_lv, facet_lv, pal, in_name, fig_stub),
    }


def _methods(eff, factors, spec):
    from importlib.metadata import version as v
    import platform
    y_label = spec["appearance"].get("y_label") or spec["data"]["y"]
    parts = "; ".join(
        f"{r.effect.replace(':', ' x ')} (F({int(r.df1)}, {int(r.df2)}) = {r.statistic:.2f}, {_plab(r.p_value)})"
        for r in eff.itertuples())
    return (
        f"{cap_first(y_label)} was analyzed by a {len(factors)}-way ANOVA (Type II sums of squares) "
        f"with factors {', '.join(factors)}. Effects tested: {parts}."
        " Significance was set at P < 0.05."
        f" Analyses were performed in Python {platform.python_version()} with statsmodels {v('statsmodels')};"
        f" figures were produced with seaborn {v('seaborn')}."
        " (The Python engine uses parametric ANOVA only; the aligned rank transform is R-only.)"
    )


def _script(spec, factors, f1_lv, f2_lv, facet_lv, pal, in_name, fig_stub):
    y = spec["data"]["y"]; f1, f2 = factors[0], factors[1]
    f3 = factors[2] if len(factors) >= 3 else None
    ap = spec["appearance"]
    facet_line = (f'facet_lv = sorted(df["{f3}"].unique())' if f3 else "facet_lv = [None]")
    return f'''#!/usr/bin/env python3
# Standalone reproduction (parametric ANOVA). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import pandas as pd, seaborn as sns
import statsmodels.api as sm, statsmodels.formula.api as smf

df = pd.read_csv("{in_name}")
factors = {factors}
for f in factors: df[f] = df[f].astype(str)
ren = {{"{y}": "_y_"}}; ren.update({{f: "_f%d_" % i for i, f in enumerate(factors)}})
d2 = df.rename(columns=ren)
rhs = "*".join("C(_f%d_)" % i for i in range(len(factors)))
aov = sm.stats.anova_lm(smf.ols("_y_ ~ " + rhs, data=d2).fit(), typ=2)
def clean(t):
    for i, f in enumerate(factors): t = t.replace("C(_f%d_)" % i, f)
    return t
def plab(p): return "P < 0.001" if p < 0.001 else "P = %.2g" % p
sig = [(clean(t), r["PR(>F)"]) for t, r in aov.iterrows() if t != "Residual" and r["PR(>F)"] < 0.05]
subtitle = "\\n".join("%s: %s" % (t.replace(":", " × "), plab(p)) for t, p in sig) or "No significant effects (P ≥ 0.05)"

f1_lv, f2_lv = {f1_lv}, {f2_lv}
{facet_line}
fig, axes = plt.subplots(1, len(facet_lv), figsize=({6.8 if f3 else 5.2}, 4.4), sharey=True, squeeze=False)
axes = axes[0]
for ax, fl in zip(axes, facet_lv):
    sub = df if fl is None else df[df["{f3}"] == fl]
    sns.boxplot(data=sub, x="{f1}", y="{y}", hue="{f2}", order=f1_lv, hue_order=f2_lv, palette={pal}, ax=ax, width=0.7, fliersize=0, legend=False)
    sns.stripplot(data=sub, x="{f1}", y="{y}", hue="{f2}", order=f1_lv, hue_order=f2_lv, dodge=True, palette={pal}, ax=ax, size=3, alpha=0.6, edgecolor="black", linewidth=0.3, legend=False)
    if fl is not None: ax.set_title(str(fl))
    ax.set_xlabel("{ap.get('x_label') or f1}"); ax.set_ylabel("{ap.get('y_label') or y}" if ax is axes[0] else "")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1].legend(handles=[mpatches.Patch(color={pal}[i], label=l) for i, l in enumerate(f2_lv)], title="{f2}", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
fig.suptitle(subtitle, fontsize=9.5, color="#4d4d4d")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
