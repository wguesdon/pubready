"""Recipe: correlation between two continuous variables. Mirrors the R recipe:
picks Pearson or Spearman from the normality of each variable (unless forced),
scatter with a linear fit and CI band, coefficient and p annotated."""
import numpy as np
import pandas as pd
import pingouin as pg
import seaborn as sns

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import cap_first, format_p, shapiro_normal


def _resolve(method, x_normal, y_normal):
    method = (method or "auto").lower()
    if method == "auto":
        m = "pearson" if (x_normal and y_normal) else "spearman"
    elif method in ("pearson", "spearman", "kendall"):
        m = method
    else:
        raise ValueError(f"unknown correlation method '{method}' "
                         "(use pearson, spearman, or kendall)")
    label = {"pearson": "Pearson correlation",
             "spearman": "Spearman rank correlation",
             "kendall": "Kendall rank correlation"}[m]
    est_name = {"pearson": "r", "spearman": "rho", "kendall": "tau"}[m]
    return dict(method=m, label=label, estimate_name=est_name)


def recipe_correlation(df, spec):
    x, y = spec["data"]["x"], spec["data"]["y"]
    if not x or not y:
        raise ValueError("correlation needs --x and --y (two numeric columns)")
    for c in (x, y):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    xv = pd.to_numeric(df[x], errors="coerce")
    yv = pd.to_numeric(df[y], errors="coerce")
    mask = np.isfinite(xv) & np.isfinite(yv)
    xv, yv = xv[mask].to_numpy(), yv[mask].to_numpy()
    df_used = pd.DataFrame({x: xv, y: yv})
    dropped = int((~mask).sum())
    steps = ([f"Dropped {dropped} row(s) with a non-numeric or missing value in "
              f"'{x}' or '{y}' ({n0} -> {len(xv)})."] if dropped
             else ["No cleaning applied; input used as-is."])
    if len(xv) < 3:
        raise ValueError("correlation needs at least 3 complete pairs")

    px, xn = shapiro_normal(xv)
    py, yn = shapiro_normal(yv)
    r = _resolve(spec["test"].get("method"), xn, yn)

    res = pg.corr(xv, yv, method=r["method"])
    estimate = float(res["r"].iloc[0])
    pval = float(res["p_val"].iloc[0])
    if "CI95" in res:
        ci = res["CI95"].iloc[0]
        ci_low, ci_high = float(ci[0]), float(ci[1])
    else:
        ci_low, ci_high = np.nan, np.nan
    if r["method"] == "pearson":
        dof = len(xv) - 2
        statistic = (estimate * np.sqrt(dof / (1 - estimate ** 2))
                     if abs(estimate) < 1 else np.inf)
    else:
        dof, statistic = np.nan, np.nan

    ap = spec["appearance"]
    pal = palette(1, ap.get("palette"))
    w, h = 4.6, 4.2
    fig, ax = new_fig(w, h)
    sns.regplot(x=xv, y=yv, ax=ax, color=pal[0],
                scatter_kws={"s": 24, "alpha": 0.8, "edgecolor": "none"},
                line_kws={"color": "#333333", "linewidth": 1.4})
    pv = format_p(pval)
    est_txt = f"{r['estimate_name']} = {estimate:.2f}, p {'' if pv.startswith('<') else '= '}{pv}"
    ax.annotate(est_txt, xy=(0.03, 0.96), xycoords="axes fraction",
                ha="left", va="top", fontsize=11)
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
        "x": x, "y": y, "n": len(xv), "method": r["label"],
        "estimate_name": r["estimate_name"], "estimate": estimate,
        "statistic": statistic, "df": dof, "p_value": pval,
        "conf_low": ci_low, "conf_high": ci_high,
        "shapiro_p_x": px, "shapiro_p_y": py,
    }])

    test_meta = {
        "name": r["label"],
        "estimate": {"name": r["estimate_name"], "value": estimate},
        "df": None if np.isnan(dof) else int(dof),
        "p_value": pval,
        "conf_int": {"low": ci_low, "high": ci_high},
        "assumptions": {"normality": {"test": "shapiro-wilk", "p_x": px, "p_y": py}},
    }
    methods = _methods(r, stats_df, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": steps,
        "label": slugify(f"{y}_vs_{x}"), "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, r, pal[0], in_name, fig_stub),
    }


def _methods(r, stats_df, spec):
    from importlib.metadata import version as v
    import platform
    x_lab = spec["appearance"].get("x_label") or stats_df["x"].iloc[0]
    y_lab = spec["appearance"].get("y_label") or stats_df["y"].iloc[0]
    auto_txt = (" The method was chosen from the Shapiro-Wilk normality of each variable."
                if (spec["test"].get("method") or "auto").lower() == "auto" else "")
    est = f"{stats_df['estimate'].iloc[0]:.2f}"
    pv = format_p(float(stats_df["p_value"].iloc[0]))
    return (
        f"The association between {cap_first(x_lab)} and {y_lab} was assessed with a "
        f"{r['label']} (n = {int(stats_df['n'].iloc[0])})." + auto_txt +
        f" The correlation was {r['estimate_name']} = {est} "
        f"(P {'' if pv.startswith('<') else '= '}{pv}), with a linear fit and 95% "
        "confidence band shown. Significance was set at P < 0.05."
        f" Analyses were performed in Python {platform.python_version()} with pingouin {v('pingouin')};"
        f" the figure was produced with seaborn {v('seaborn')}."
    )


def _script(spec, r, color, in_name, fig_stub):
    x, y = spec["data"]["x"], spec["data"]["y"]
    ap = spec["appearance"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, seaborn as sns, pingouin as pg

raw = pd.read_csv("{in_name}")
x = pd.to_numeric(raw["{x}"], errors="coerce")
y = pd.to_numeric(raw["{y}"], errors="coerce")
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask].to_numpy(), y[mask].to_numpy()
print(pg.corr(x, y, method="{r['method']}"))

fig, ax = plt.subplots(figsize=(4.6, 4.2))
sns.regplot(x=x, y=y, ax=ax, color="{color}",
            scatter_kws={{"s": 24, "alpha": 0.8, "edgecolor": "none"}},
            line_kws={{"color": "#333333", "linewidth": 1.4}})
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
