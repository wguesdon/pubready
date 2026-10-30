"""Recipe: two-way / three-way factorial ANOVA (parametric, statsmodels Type II).
The Python engine has no standard aligned rank transform, so unlike the R engine
this path is parametric only. Grouped box, faceted by the third factor, with each
effect's p-value shown on top. ``--geom interaction`` draws the interaction plot
of the same ANOVA instead: the mean +/- SEM of each cell, one line per --fill
level."""
import functools
import textwrap

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from ..theme import apply_pub_style, palette, rain_offsets
from ..util import (SCRIPT_LABEL_FN, cap_first, check_spread, check_test_method, level_labels,
                    level_order)

# The dodge width of the interaction plot, as position_dodge(0.25) in R.
INT_DODGE = 0.25
# The dodge width of the bar geom, as position_dodge(0.8) in R.
BAR_DODGE = 0.8
# The dodge width of the box geom: seaborn dodges the boxes over their width of 0.7.
BOX_DODGE = 0.7


TESTS = ("auto", "anova", "aov", "art", "aligned_rank")
"""The --test values that factorial_anova knows (F63)."""


def _plab(p):
    """Format a P value for the subtitle and the methods.

    Args:
        p: The P value.

    Returns:
        str: The label.
    """
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "P = NA"
    return "P < 0.001" if p < 0.001 else f"P = {p:.2g}"


def clean_step_lines(y, coerced, n_missing, n0, n1):
    """The data log of the cleaning, as .fa_clean_steps() of the R engine writes it.

    Args:
        y: The outcome column.
        coerced: The count of outcome values that became NA as numbers.
        n_missing: The count of missing values of each column, in order.
        n0: The count of rows before the cleaning.
        n1: The count of rows after it.

    Returns:
        list: The steps: the coerced values, the rows with no value in each
        column (a row with two missing values counts in both), and the rows
        dropped in all.
    """
    steps = []
    if coerced:
        steps.append(f"Coerced column '{y}' to numeric; {coerced} value(s) became NA.")
    for col, n in n_missing.items():
        if n:
            steps.append(f"{n} row(s) have no value in '{col}'.")
    if n1 < n0:
        steps.append(f"Dropped {n0 - n1} row(s) with a missing value ({n0} -> {n1} rows).")
    return steps or ["No cleaning applied; input used as-is."]


def _clean_term(t, factors):
    """Turn a statsmodels term such as C(_f0_):C(_f1_) back into column names.

    Args:
        t: The term.
        factors: The factor columns, in model order.

    Returns:
        str: The term with the column names.
    """
    for i, f in enumerate(factors):
        t = t.replace(f"C(_f{i}_)", f)
    return t


def wrap_x(label):
    """Wrap an --x level on spaces at 12 characters, the lines of the R engine.

    Args:
        label: The level.

    Returns:
        str: The wrapped level.
    """
    return textwrap.fill(str(label), 12, break_on_hyphens=False, break_long_words=False)


def needs_wrap(levels):
    """Whether any --x level is longer than 12 characters.

    Args:
        levels: The levels.

    Returns:
        bool: True when the labels are wrapped.
    """
    return any(len(str(lv)) > 12 for lv in levels)


def fig_size(x_lv, fill_lv, fill_name, n_facet, has_facet):
    """The figure size, with the rule of the R engine.

    The width holds every --x level times every --fill level in each panel and
    the legend, and is never below 5.2 (6.8 with a facet) inches. A wrapped
    label adds a line of height.

    Args:
        x_lv: The --x levels.
        fill_lv: The --fill levels.
        fill_name: The --fill column, the legend title.
        n_facet: The number of panels.
        has_facet: Whether there is a --facet column.

    Returns:
        tuple: The width and the height in inches.
    """
    slot = 0.3 * len(fill_lv) + 0.35
    legend_w = 0.5 + 0.08 * max(len(str(s)) for s in list(fill_lv) + [fill_name])
    w = max(6.8 if has_facet else 5.2, 0.9 + n_facet * len(x_lv) * slot + legend_w)
    lines = max(wrap_x(s).count("\n") + 1 for s in x_lv)
    h = 4.4 + (0.18 * (lines - 1) if needs_wrap(x_lv) else 0.0)
    return w, h


def _interaction_panel(ax, sub, y, f1, f2, f1_lv, f2_lv, pal):
    """Draw the interaction plot of one panel: mean +/- SEM per cell.

    Args:
        ax: The axes.
        sub: The rows of the panel.
        y: The outcome column.
        f1: The --x column.
        f2: The --fill column.
        f1_lv: The --x levels.
        f2_lv: The --fill levels.
        pal: The colour of each --fill level.
    """
    k = len(f2_lv)
    base = np.arange(len(f1_lv), dtype=float)
    for i, lv in enumerate(f2_lv):
        g = sub[sub[f2] == lv].groupby(f1)[y]
        means = g.mean().reindex(f1_lv).to_numpy(dtype=float)
        sems = (g.std(ddof=1) / np.sqrt(g.size())).reindex(f1_lv).to_numpy(dtype=float)
        # The offset of position_dodge(width): (i + 0.5 - k / 2) * width / k.
        x = base + (i + 0.5 - k / 2) * INT_DODGE / k
        ok = np.isfinite(means)
        ax.plot(x[ok], means[ok], color=pal[i], lw=1.6, zorder=2)
        ax.errorbar(x[ok], means[ok], yerr=np.nan_to_num(sems[ok]), fmt="none", ecolor=pal[i],
                    elinewidth=1.1, capsize=3, zorder=2)
        ax.scatter(x[ok], means[ok], s=36, color=pal[i], zorder=3)
    ax.set_xticks(base)
    ax.set_xlim(-0.6, len(f1_lv) - 0.4)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 3, 5, 10]))


def _bar_panel(ax, sub, y, f1, f2, f1_lv, f2_lv, pal):
    """Draw the bar geom of one panel: mean bar +/- SEM per cell, as R's stat_summary.

    The bars are dodged as position_dodge(0.8) with a bar width of 0.7, so the
    layout matches the R engine.

    Args:
        ax: The axes.
        sub: The rows of the panel.
        y: The outcome column.
        f1: The --x column.
        f2: The --fill column.
        f1_lv: The --x levels.
        f2_lv: The --fill levels.
        pal: The colour of each --fill level.

    Returns:
        list: The x position of each --fill level's bars.
    """
    k = len(f2_lv)
    base = np.arange(len(f1_lv), dtype=float)
    bar_w = 0.7 / k
    centres = []
    for i, lv in enumerate(f2_lv):
        g = sub[sub[f2] == lv].groupby(f1)[y]
        means = g.mean().reindex(f1_lv).to_numpy(dtype=float)
        sems = (g.std(ddof=1) / np.sqrt(g.size())).reindex(f1_lv).to_numpy(dtype=float)
        x = base + (i + 0.5 - k / 2) * BAR_DODGE / k
        ok = np.isfinite(means)
        ax.bar(x[ok], means[ok], width=bar_w, color=pal[i], alpha=0.9, edgecolor="black",
               linewidth=0.3, zorder=1)
        ax.errorbar(x[ok], means[ok], yerr=np.nan_to_num(sems[ok]), fmt="none", ecolor="black",
                    elinewidth=0.8, capsize=3, zorder=2)
        centres.append(x)
    ax.set_xticks(base)
    ax.set_xlim(-0.6, len(f1_lv) - 0.4)
    return centres


def _cell_points(ax, sub, y, f1, f2, f1_lv, f2_lv, pal, dodge):
    """Draw the points of each cell at fixed offsets around the cell centre.

    The centre is the dodged position of the box or the bar of the cell. Each
    point gets a golden-ratio offset by its row in the cell, within
    +/- 0.15 / (k + 2) for k --fill levels, the jitter of
    position_jitterdodge(jitter.width = 0.15) in the R engine. The points keep
    their places in every run, as .fa_point_x() does in R.

    Args:
        ax: The axes.
        sub: The rows of the panel.
        y: The outcome column.
        f1: The --x column.
        f2: The --fill column.
        f1_lv: The --x levels.
        f2_lv: The --fill levels.
        pal: The colour of each --fill level.
        dodge: The dodge width of the boxes or the bars.
    """
    k = len(f2_lv)
    a = 0.15 / (k + 2)
    for i, lv in enumerate(f2_lv):
        for j, v in enumerate(f1_lv):
            vals = sub.loc[(sub[f2] == lv) & (sub[f1] == v), y].to_numpy(dtype=float)
            if len(vals) == 0:
                continue
            x = j + (i + 0.5 - k / 2) * dodge / k + rain_offsets(len(vals), -a, a)
            ax.scatter(x, vals, s=10, color=pal[i], alpha=0.7, edgecolors="black",
                       linewidths=0.3, zorder=3)


def recipe_factorial_anova(df, spec):
    """Run a two-way or three-way Type II ANOVA and draw the grouped figure.

    Args:
        df: The tidy data.
        spec: The figure specification.

    Returns:
        dict: The recipe result for the bundle writer.
    """
    check_test_method("factorial_anova", spec["test"].get("method"), TESTS)
    y = spec["data"]["y"]
    f1, f2, f3 = spec["data"]["x"], spec["data"]["fill"], spec["data"].get("facet")
    if not f2:
        raise ValueError("factorial_anova needs at least two factors: pass --x and --fill")
    method = (spec["test"].get("method") or "auto").lower()
    if method in ("art", "aligned_rank"):
        raise ValueError("the aligned rank transform (non-parametric factorial) is R-only; "
                         "rerun with --engine r")
    factors = [f1, f2] + ([f3] if f3 else [])
    for f in factors + [y]:
        if f not in df.columns:
            raise ValueError(f"column '{f}' not found")
    df = df.copy()
    for f in factors:
        df[f] = level_labels(df[f])
    # Each reason for a dropped row is counted on its own: a text value of the
    # outcome that becomes NA, and a missing value in each column.
    y_raw = df[y]
    df[y] = pd.to_numeric(df[y], errors="coerce")
    coerced = int((df[y].isna() & y_raw.notna()).sum())
    n_missing = {c: int(df[c].isna().sum()) for c in [y] + factors}
    n0 = len(df)
    df = df.dropna(subset=[y] + factors).reset_index(drop=True)
    clean_steps = clean_step_lines(y, coerced, n_missing, n0, len(df))
    check_spread(df[y], "factorial_anova", f"the outcome '{y}'")

    ren = {y: "_y_"}
    ren.update({f: f"_f{i}_" for i, f in enumerate(factors)})
    d2 = df.rename(columns=ren)
    rhs = "*".join(f"C(_f{i}_)" for i in range(len(factors)))
    model = smf.ols(f"_y_ ~ {rhs}", data=d2).fit()
    aov = sm.stats.anova_lm(model, typ=2)
    # Partial eta squared, SS_effect / (SS_effect + SS_residual), as the R engine reports.
    ss_res = float(aov.loc["Residual", "sum_sq"])

    rows = []
    for term, r in aov.iterrows():
        if term == "Residual":
            continue
        rows.append({"effect": _clean_term(term, factors), "df1": int(round(r["df"])),
                     "df2": int(round(aov.loc["Residual", "df"])), "statistic": float(r["F"]),
                     "p_value": float(r["PR(>F)"]),
                     "effect_size": float(r["sum_sq"]) / (float(r["sum_sq"]) + ss_res),
                     "effect_size_name": "partial eta squared",
                     "significant": float(r["PR(>F)"]) < 0.05})
    eff = pd.DataFrame(rows)

    sig = eff[eff["significant"]]
    subtitle = ("\n".join(f"{r.effect.replace(':', ' × ')}: {_plab(r.p_value)}"
                          for r in sig.itertuples())
                if len(sig) else "No significant effects (P ≥ 0.05)")

    ap = spec["appearance"]
    # --levels orders the --x levels; --fill and --facet take the rule without a
    # list: numeric order when every label is a number, else code-point order.
    f1_lv = level_order(df[f1], spec["data"].get("levels"), "level")
    f2_lv = level_order(df[f2], None, "level")
    pal = palette(len(f2_lv), ap.get("palette"))
    facet_lv = level_order(df[f3], None, "level") if f3 else [None]
    geom = ap.get("geom") or "box"
    w, h = fig_size(f1_lv, f2_lv, f2, len(facet_lv), bool(f3))
    wrap = needs_wrap(f1_lv)
    fig, axes = plt.subplots(1, len(facet_lv), figsize=(w, h), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, fl in zip(axes, facet_lv):
        sub = df if fl is None else df[df[f3] == fl]
        if geom == "interaction":
            _interaction_panel(ax, sub, y, f1, f2, f1_lv, f2_lv, pal)
        elif geom == "bar":
            _bar_panel(ax, sub, y, f1, f2, f1_lv, f2_lv, pal)
            if ap.get("show_points", True):
                _cell_points(ax, sub, y, f1, f2, f1_lv, f2_lv, pal, BAR_DODGE)
        else:
            sns.boxplot(data=sub, x=f1, y=y, hue=f2, order=f1_lv, hue_order=f2_lv,
                        palette=pal, ax=ax, width=BOX_DODGE, fliersize=0, legend=False)
            if ap.get("show_points", True):
                _cell_points(ax, sub, y, f1, f2, f1_lv, f2_lv, pal, BOX_DODGE)
        if wrap or geom in ("interaction", "bar"):
            ax.set_xticklabels([wrap_x(s) if wrap else str(s) for s in f1_lv])
        if fl is not None:
            ax.set_title(str(fl), fontsize=11)
        ax.set_xlabel(ap.get("x_label") or f1)
        ax.set_ylabel(ap.get("y_label") or y if ax is axes[0] else "")
        apply_pub_style(ax, ap.get("theme"))
    if geom == "interaction":
        handles = [Line2D([0], [0], color=pal[i], lw=1.6, marker="o", ms=6, label=lvl)
                   for i, lvl in enumerate(f2_lv)]
    else:
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
        "qc": {"quantity": "model residuals", "panels": {"residuals": model.resid.to_numpy()}},
        "build_script": functools.partial(_script, spec, factors, f1_lv, f2_lv, facet_lv, pal, w, h),
    }


def _methods(eff, factors, spec):
    """Write the methods paragraph.

    Args:
        eff: The effects table.
        factors: The factor columns.
        spec: The figure specification.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    y_label = spec["appearance"].get("y_label") or spec["data"]["y"]
    parts = "; ".join(
        f"{r.effect.replace(':', ' x ')} (F({int(r.df1)}, {int(r.df2)}) = {r.statistic:.2f}, {_plab(r.p_value)})"
        for r in eff.itertuples())
    return (
        f"{cap_first(y_label)} was analyzed by a {len(factors)}-way ANOVA (Type II sums of squares) "
        f"with factors {', '.join(factors)}. Effects tested: {parts}. Effect sizes are partial eta squared."
        " Significance was set at P < 0.05."
        f" Analyses were performed in Python {platform.python_version()} with statsmodels {v('statsmodels')};"
        f" figures were produced with seaborn {v('seaborn')}."
        " (The Python engine uses parametric ANOVA only; the aligned rank transform is R-only.)"
    )


def _script_points(f1, f2, y, pal, dodge):
    """Write the code that draws the points of each cell in an emitted script.

    The code places the points as _cell_points does: the dodged cell centre plus
    a golden-ratio offset by row in the cell. It runs inside the panel loop.

    Args:
        f1: The --x column.
        f2: The --fill column.
        y: The outcome column.
        pal: The colour of each --fill level.
        dodge: The dodge width of the boxes or the bars.

    Returns:
        str: The code lines, indented for the panel loop.
    """
    return (
        "    # Each point sits at its dodged cell centre plus a fixed golden-ratio offset.\n"
        "    k = len(f2_lv); a = 0.15 / (k + 2)\n"
        "    for i, lv in enumerate(f2_lv):\n"
        "        for j, v in enumerate(f1_lv):\n"
        f'            vals = sub.loc[(sub["{f2}"] == lv) & (sub["{f1}"] == v), "{y}"].to_numpy(dtype=float)\n'
        "            n = np.arange(1, len(vals) + 1)\n"
        f"            x = j + (i + 0.5 - k / 2) * {dodge} / k - a + 2 * a * np.mod(n * 0.6180339887498949, 1.0)\n"
        f'            ax.scatter(x, vals, s=10, color={pal}[i], alpha=0.7, edgecolors="black", '
        "linewidths=0.3, zorder=3)"
    )


def _script(spec, factors, f1_lv, f2_lv, facet_lv, pal, w, h, in_name, fig_stub):
    """Write the standalone reproduction script.

    Args:
        spec: The figure specification.
        factors: The factor columns.
        f1_lv: The --x levels.
        f2_lv: The --fill levels.
        facet_lv: The --facet levels, or [None].
        pal: The colour of each --fill level.
        w: The figure width in inches.
        h: The figure height in inches.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script.
    """
    y = spec["data"]["y"]; f1, f2 = factors[0], factors[1]
    f3 = factors[2] if len(factors) >= 3 else None
    ap = spec["appearance"]
    facet_line = f"facet_lv = {facet_lv!r}"
    wrap_line = ('\n    ax.set_xticks(range(len(f1_lv))); ax.set_xticklabels([textwrap.fill(s, 12, break_on_hyphens=False) for s in f1_lv])'
                 if needs_wrap(f1_lv) else "")
    if (ap.get("geom") or "box") == "interaction":
        draw = f'''    k = len(f2_lv)
    for i, lv in enumerate(f2_lv):
        g = sub[sub["{f2}"] == lv].groupby("{f1}")["{y}"]
        m = g.mean().reindex(f1_lv); se = (g.std() / g.size() ** 0.5).reindex(f1_lv)
        x = [j + (i + 0.5 - k / 2) * 0.25 / k for j in range(len(f1_lv))]
        ax.errorbar(x, m, yerr=se, color={pal}[i], marker="o", capsize=3, label=lv)
    ax.set_xticks(range(len(f1_lv))); ax.set_xticklabels(f1_lv)'''
        legend = f'axes[-1].legend(title="{f2}", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)'
    elif (ap.get("geom") or "box") == "bar":
        draw = f'''    k = len(f2_lv)
    for i, lv in enumerate(f2_lv):
        g = sub[sub["{f2}"] == lv].groupby("{f1}")["{y}"]
        m = g.mean().reindex(f1_lv); se = (g.std() / g.size() ** 0.5).reindex(f1_lv)
        x = [j + (i + 0.5 - k / 2) * 0.8 / k for j in range(len(f1_lv))]
        ax.bar(x, m, width=0.7 / k, color={pal}[i], alpha=0.9, edgecolor="black", linewidth=0.3)
        ax.errorbar(x, m, yerr=se, fmt="none", ecolor="black", capsize=3)
    ax.set_xticks(range(len(f1_lv))); ax.set_xticklabels(f1_lv)'''
        if ap.get("show_points", True):
            draw += "\n" + _script_points(f1, f2, y, pal, BAR_DODGE)
        legend = f'axes[-1].legend(handles=[mpatches.Patch(color={pal}[i], label=l) for i, l in enumerate(f2_lv)], title="{f2}", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)'
    else:
        draw = f'''    sns.boxplot(data=sub, x="{f1}", y="{y}", hue="{f2}", order=f1_lv, hue_order=f2_lv, palette={pal}, ax=ax, width={BOX_DODGE}, fliersize=0, legend=False)'''
        if ap.get("show_points", True):
            draw += "\n" + _script_points(f1, f2, y, pal, BOX_DODGE)
        legend = f'axes[-1].legend(handles=[mpatches.Patch(color={pal}[i], label=l) for i, l in enumerate(f2_lv)], title="{f2}", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)'
    imports = "import textwrap\n" if needs_wrap(f1_lv) else ""
    return f'''#!/usr/bin/env python3
# Standalone reproduction (parametric ANOVA). Run inside the pinned container.
{imports}import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import numpy as np
import pandas as pd, seaborn as sns
import statsmodels.api as sm, statsmodels.formula.api as smf


{SCRIPT_LABEL_FN}

df = pd.read_csv("{in_name}")
factors = {factors}
for f in factors: df[f] = df[f].map(label)
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
fig, axes = plt.subplots(1, len(facet_lv), figsize=({w:g}, {h:g}), sharey=True, squeeze=False)
axes = axes[0]
for ax, fl in zip(axes, facet_lv):
    sub = df if fl is None else df[df["{f3}"] == fl]
{draw}{wrap_line}
    if fl is not None: ax.set_title(str(fl))
    ax.set_xlabel("{ap.get('x_label') or f1}"); ax.set_ylabel("{ap.get('y_label') or y}" if ax is axes[0] else "")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
{legend}
fig.suptitle(subtitle, fontsize=9.5, color="#4d4d4d")
fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
