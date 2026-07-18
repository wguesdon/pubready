"""Recipe: Kaplan-Meier survival curves compared between groups (lifelines).
Step curves with censoring ticks, the log-rank p-value on the plot, and a
number-at-risk table below. Mirrors the R survival_km."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test

from ..theme import apply_pub_style, palette
from ..util import format_p


def recipe_survival_km(df, spec):
    t, e, g = spec["data"]["time"], spec["data"]["event"], spec["data"]["x"]
    if not t or not e:
        raise ValueError("survival_km needs --time and --event columns")
    if not g:
        raise ValueError("survival_km needs a grouping column via --x")
    for c in (t, e, g):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found")

    df = df.copy()
    df[t] = pd.to_numeric(df[t], errors="coerce")
    ev = pd.to_numeric(df[e], errors="coerce")
    if not ev.dropna().isin([0, 1]).all():
        raise ValueError(f"event column '{e}' must be coded 1 = event, 0 = censored")
    df[e] = ev
    n0 = len(df)
    df = df.dropna(subset=[t, e, g])
    df = df[df[t] >= 0].reset_index(drop=True)
    clean_steps = (f"Dropped {n0 - len(df)} row(s) with missing/invalid values ({n0} -> {len(df)} rows)."
                   if len(df) < n0 else "No cleaning applied; input used as-is.")
    df[g] = df[g].astype(str)
    lv = sorted(df[g].unique())
    pal = palette(len(lv), spec["appearance"].get("palette"))

    w, h = 5.6, 5.8
    fig, ax = plt.subplots(figsize=(w, h))
    kmfs, rows = [], []
    for i, grp in enumerate(lv):
        sub = df[df[g] == grp]
        kmf = KaplanMeierFitter()
        kmf.fit(sub[t], sub[e], label=grp)
        kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i],
                                   show_censors=True, censor_styles={"marker": "|", "ms": 6})
        kmfs.append(kmf)
        rows.append({"group": grp, "n": int(len(sub)), "events": int(sub[e].sum()),
                     "median_survival": float(kmf.median_survival_time_)})

    res = multivariate_logrank_test(df[t], df[g], df[e])
    lr_p = float(res.p_value)
    ax.text(0.04, 0.08, f"log-rank p {'< 0.001' if lr_p < 0.001 else '= ' + format(lr_p, '.2g')}",
            transform=ax.transAxes, fontsize=11)
    ax.set_xlabel(spec["appearance"].get("x_label") or "Time")
    ax.set_ylabel(spec["appearance"].get("y_label") or "Survival probability")
    ax.set_ylim(0, 1.02)
    ax.legend(title=g, frameon=False, fontsize=10)
    apply_pub_style(ax)
    add_at_risk_counts(*kmfs, ax=ax)
    fig.tight_layout()

    stats_df = pd.DataFrame(rows)
    stats_df["log_rank_p"] = lr_p
    lr_df = len(lv) - 1
    test_meta = {
        "name": "Kaplan-Meier estimate with log-rank test",
        "log_rank": {"chisq": float(res.test_statistic), "df": lr_df, "p_value": lr_p},
        "groups": rows,
    }
    methods = _methods(rows, res, lr_df, lr_p, g, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "label": g, "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, t, e, g, lv, pal, in_name, fig_stub),
    }


def _methods(rows, res, lr_df, lr_p, g, spec):
    from importlib.metadata import version as v
    import platform
    med = ", ".join(f"{r['group']} (median {r['median_survival']:.3g})" for r in rows)
    return (
        "Survival was estimated by the Kaplan-Meier method and compared between "
        f"groups of {g} with the log-rank test (chi-square({lr_df}) = {res.test_statistic:.2f}, "
        f"P = {format_p(lr_p)})."
        f" Median survival by group: {med}."
        " Censoring is marked with ticks and the number at risk is tabulated below the curves."
        f" Analyses were performed in Python {platform.python_version()} with lifelines {v('lifelines')}."
    )


def _script(spec, t, e, g, lv, pal, in_name, fig_stub):
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test

df = pd.read_csv("{in_name}"); df["{g}"] = df["{g}"].astype(str)
order, pal = {lv}, {pal}
fig, ax = plt.subplots(figsize=(5.6, 5.8))
kmfs = []
for i, grp in enumerate(order):
    sub = df[df["{g}"] == grp]
    kmf = KaplanMeierFitter().fit(sub["{t}"], sub["{e}"], label=grp)
    kmf.plot_survival_function(ax=ax, ci_show=False, color=pal[i], show_censors=True, censor_styles={{"marker": "|", "ms": 6}})
    kmfs.append(kmf)
res = multivariate_logrank_test(df["{t}"], df["{g}"], df["{e}"])
p = res.p_value
ax.text(0.04, 0.08, "log-rank p " + ("< 0.001" if p < 0.001 else "= %.2g" % p), transform=ax.transAxes, fontsize=11)
ax.set_xlabel("{spec['appearance'].get('x_label') or 'Time'}"); ax.set_ylabel("{spec['appearance'].get('y_label') or 'Survival probability'}")
ax.set_ylim(0, 1.02); ax.legend(title="{g}", frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
add_at_risk_counts(*kmfs, ax=ax)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
