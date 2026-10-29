"""Recipe: compare proportions of a categorical outcome across groups. Mirrors
the R recipe: a chi-squared test of independence, or Fisher's exact when any
expected cell is below 5 (or forced); an odds ratio for a 2x2 table. Figure: a
100% stacked bar of the outcome proportions within each group. Fisher's test
follows stats::fisher.test (pubready.fisher) on any table."""
import functools
import inspect
import os

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
from scipy import stats as ss

from .. import fisher as fisher_module
from ..bundle import slugify
from ..fisher import fisher_2x2, fisher_rxc_p
from ..theme import apply_pub_style, category_width, legend_width, palette, wrap_labels


def group_order(g, levels=None):
    """Group order: the order of --levels when it is given, else first appearance.

    Args:
        g: The group values as text.
        levels: The groups in the order that --levels gives, or None. A listed
            group that no row holds is left out.

    Returns:
        list: The groups in drawing order.

    Raises:
        ValueError: If a group of the data is missing from levels.
    """
    u = list(pd.unique(g))
    if not levels:
        return u
    levels = list(dict.fromkeys(str(v).strip() for v in levels))
    miss = [v for v in u if v not in levels]
    if miss:
        raise ValueError("--levels does not list the group(s) "
                         + ", ".join(f"'{v}'" for v in miss) + " of the data")
    return [v for v in levels if v in u]


def recipe_proportions(df, spec):
    x, y = spec["data"]["x"], spec["data"]["y"]
    if not x or not y:
        raise ValueError("proportions needs --x (group) and --y (outcome)")
    for c in (x, y):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    d = df[[x, y]].copy()
    n0 = len(d)
    d = d.dropna()
    d[x] = d[x].astype(str); d[y] = d[y].astype(str)
    groups = group_order(d[x], spec["data"].get("levels"))
    outcomes = list(pd.unique(d[y]))
    steps = ([f"Dropped {n0 - len(d)} row(s) with a missing group or outcome ({n0} -> {len(d)})."]
             if len(d) < n0 else ["No cleaning applied; input used as-is."])
    tab = pd.crosstab(d[x], d[y]).reindex(index=groups, columns=outcomes).fillna(0)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        raise ValueError("proportions needs at least 2 groups and 2 outcomes")

    method = (spec["test"].get("method") or "auto").lower()
    chi2, p_chi, dof, expected = ss.chi2_contingency(tab.to_numpy(), correction=False)
    small = bool((expected < 5).any())
    use_fisher = small if method in ("auto", "") else method in ("fisher", "fisher_exact")
    if method in ("chisq", "chi", "chi_square", "chisquare"):
        use_fisher = False
    is2x2 = tab.shape == (2, 2)

    orr, or_ci = np.nan, (np.nan, np.nan)
    if use_fisher:
        # stats::fisher.test: the conditional MLE odds ratio and its exact CI for
        # a 2 x 2 table, the network-algorithm P value for a larger table.
        stat, dof_out, test_label = np.nan, np.nan, "Fisher's exact test"
        counts = tab.to_numpy().astype(int)
        if is2x2:
            ft = fisher_2x2(counts)
            pval, orr, or_ci = ft["p_value"], ft["odds_ratio"], (ft["ci_low"], ft["ci_high"])
        else:
            pval = fisher_rxc_p(counts)
    else:
        pval, stat, dof_out, test_label = p_chi, chi2, dof, "Pearson's chi-squared test"
        if is2x2:
            a, b, c2, d2 = tab.to_numpy().ravel()
            orr = (a * d2) / (b * c2) if (b * c2) else np.nan

    prop = tab.div(tab.sum(axis=1), axis=0)
    ap = spec["appearance"]
    pal = palette(len(outcomes), ap.get("palette"))
    # The axis, a slot per group, and the legend, as in the R engine.
    w = max(3.6, category_width(len(groups), wrap_labels(groups), 0) + legend_width(outcomes, y))
    h = 4.2
    fig, ax = plt.subplots(figsize=(w, h))
    xpos = np.arange(len(groups))
    bottom = np.zeros(len(groups))
    # The first outcome goes on top of the stack and the last at the bottom, as
    # ggplot2 stacks them in R, so the stack reads in the order of the legend.
    handles = {}
    for j in reversed(range(len(outcomes))):
        oc = outcomes[j]
        vals = prop[oc].to_numpy()
        handles[oc] = ax.bar(xpos, vals, bottom=bottom, width=0.7, color=pal[j],
                             edgecolor="white", label=oc)
        for i, v in enumerate(vals):
            if v >= 0.04:
                ax.text(xpos[i], bottom[i] + v / 2, f"{v * 100:.0f}%",
                        ha="center", va="center", color="white", fontsize=9)
        bottom += vals
    ax.set_xticks(xpos); ax.set_xticklabels(wrap_labels(groups))
    ax.set_ylim(0, 1); ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or "proportion")
    # The pad keeps the subtitle clear of the 100% tick label.
    ax.set_title(f"{test_label}, P = {pval:.2g}", fontsize=10, color="#4d4d4d", pad=12)
    apply_pub_style(ax, ap.get("theme"))
    ax.legend([handles[oc] for oc in outcomes], outcomes, title=y, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    fig.tight_layout()

    stats_df = pd.DataFrame([
        {"group": g, "outcome": oc, "n": int(tab.loc[g, oc]), "proportion": float(prop.loc[g, oc]),
         "omnibus_test": test_label, "omnibus_statistic": stat, "omnibus_df1": dof_out,
         "omnibus_df2": np.nan, "omnibus_p": float(pval), "odds_ratio": orr}
        for g in groups for oc in outcomes])
    test_meta = {"name": test_label, "statistic": None if np.isnan(stat) else float(stat),
                 "df": None if dof_out is None or (isinstance(dof_out, float) and np.isnan(dof_out)) else int(dof_out),
                 "p_value": float(pval), "odds_ratio": None if np.isnan(orr) else float(orr),
                 "or_conf_low": None if np.isnan(or_ci[0]) else or_ci[0],
                 "or_conf_high": None if np.isnan(or_ci[1]) else or_ci[1],
                 "n": int(tab.to_numpy().sum()), "any_expected_below_5": small}
    methods = _methods(test_label, tab.shape[0], tab.shape[1], int(tab.to_numpy().sum()), pval, orr,
                       is2x2, stat, dof_out, list(tab.index), list(tab.columns))
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": d.reset_index(drop=True), "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, use_fisher, is2x2),
    }


def _methods(test_label, n_grp, n_out, n, pval, orr, is2x2, stat, dof, groups, outcomes):
    from importlib.metadata import version as v
    import platform
    # The odds ratio is the odds of the first outcome in the first group over the
    # odds in the second group, so the text names the outcome and both groups.
    or_txt = (f" The odds ratio of {outcomes[0]} for {groups[0]} versus {groups[1]} was {orr:.2f}."
              if (is2x2 and np.isfinite(orr)) else "")
    p_txt = f"P = {pval:.2g}"
    if stat is not None and np.isfinite(stat):
        p_txt = f"chi-square({int(dof)}) = {stat:.2f}, {p_txt}"
    return (
        f"The distribution of the outcome ({n_out} categories) across {n_grp} groups (n = {n}) "
        f"was tested with {test_label} ({p_txt}).{or_txt} Proportions are shown as a 100% "
        "stacked bar per group. "
        f"Analyses were performed in Python {platform.python_version()} with scipy {v('scipy')}; "
        f"the figure was produced with matplotlib {v('matplotlib')}."
    )


def _script(spec, fisher, is2x2, in_name, fig_stub):
    x, y = spec["data"]["x"], spec["data"]["y"]
    if fisher:
        # The Fisher test of stats::fisher.test, carried whole into the script.
        call = "fisher_2x2" if is2x2 else "fisher_rxc_p"
        test_line = (inspect.getsource(fisher_module)
                     + f"\nprint({call}(tab.to_numpy().astype(int)))")
    else:
        test_line = 'print(ss.chi2_contingency(tab.to_numpy(), correction=False))'
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np, pandas as pd
from scipy import stats as ss

raw = pd.read_csv("{in_name}").dropna(subset=["{x}", "{y}"])
tab = pd.crosstab(raw["{x}"].astype(str), raw["{y}"].astype(str))
{test_line}
prop = tab.div(tab.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(4.2, 4.2))
xpos = np.arange(tab.shape[0]); bottom = np.zeros(tab.shape[0])
for oc in prop.columns:
    ax.bar(xpos, prop[oc].to_numpy(), bottom=bottom, width=0.7, edgecolor="white", label=oc)
    bottom += prop[oc].to_numpy()
ax.set_xticks(xpos); ax.set_xticklabels(list(prop.index))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0)); ax.set_ylabel("proportion")
ax.legend(title="{y}", frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
