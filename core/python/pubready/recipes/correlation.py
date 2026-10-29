"""Recipe: correlation between two continuous variables. Mirrors the R recipe:
picks Pearson or Spearman from the normality of each variable (unless forced),
scatter with a linear fit and CI band, method, coefficient and P in a subtitle."""
import functools
import inspect

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

from ..bundle import slugify
from ..theme import R_POINT_PT, apply_pub_style, new_fig, palette
from ..util import cap_first, shapiro_normal


def lm_band(x, y, n=80):
    """The least-squares line and its 95% confidence band over the range of x.

    The band is the t-based interval of the fitted mean, the band that
    geom_smooth(method = "lm") draws in the R engine on 80 points. It has no
    random part, so the figure is the same in every run.

    Args:
        x: The x values.
        y: The y values.
        n: The count of grid points.

    Returns:
        tuple: The grid, the fitted line, the lower bound and the upper bound.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    grid = np.linspace(x.min(), x.max(), n)
    fit = intercept + slope * grid
    dof = len(x) - 2
    s = np.sqrt(np.sum((y - (intercept + slope * x)) ** 2) / dof)
    se = s * np.sqrt(1 / len(x) + (grid - x.mean()) ** 2 / np.sum((x - x.mean()) ** 2))
    t = stats.t.ppf(0.975, dof)
    return grid, fit, fit - t * se, fit + t * se


def cor_label(method_label, estimate_name, estimate, p):
    """Write the subtitle, as cor_label() of the R engine does.

    Args:
        method_label: The method, for example "Pearson correlation".
        estimate_name: The coefficient symbol: r, rho or tau.
        estimate: The coefficient.
        p: The P value.

    Returns:
        str: For example "Pearson correlation, r = 0.71, P = 2.1e-10".
    """
    return f"{method_label}, {estimate_name} = {estimate:.2f}, P = {p:.2g}"


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
    # The points, then the line and its analytic 95% band on top, as R draws
    # geom_point(size = 2) and geom_smooth(method = "lm", fill = "grey75").
    ax.scatter(xv, yv, s=(2 * R_POINT_PT) ** 2, color=pal[0], alpha=0.8, edgecolors="none")
    grid, fit, lo, hi = lm_band(xv, yv)
    ax.fill_between(grid, lo, hi, color="#BFBFBF", alpha=0.4, linewidth=0)
    ax.plot(grid, fit, color="#333333", linewidth=1.4)
    # The method, the coefficient and the P value go in a subtitle, in the style
    # of the other recipes and in the same words as the R engine.
    ax.set_title(cor_label(r["label"], r["estimate_name"], estimate, pval),
                 fontsize=10, color="#4d4d4d")
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
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
        "qc": {"quantity": "variable values", "panels": {x: xv, y: yv}},
        "build_script": functools.partial(_script, spec, r, pal[0]),
    }


def _methods(r, stats_df, spec):
    from importlib.metadata import version as v
    import platform
    x_lab = spec["appearance"].get("x_label") or stats_df["x"].iloc[0]
    y_lab = spec["appearance"].get("y_label") or stats_df["y"].iloc[0]
    auto_txt = (" The method was chosen from the Shapiro-Wilk normality of each variable."
                if (spec["test"].get("method") or "auto").lower() == "auto" else "")
    est = f"{stats_df['estimate'].iloc[0]:.2f}"
    pv = f"{float(stats_df['p_value'].iloc[0]):.2g}"
    return (
        f"The association between {cap_first(x_lab)} and {y_lab} was assessed with a "
        f"{r['label']} (n = {int(stats_df['n'].iloc[0])})." + auto_txt +
        f" The correlation was {r['estimate_name']} = {est} (P = {pv}), with a linear fit and "
        "95% confidence band shown. Significance was set at P < 0.05."
        f" Analyses were performed in Python {platform.python_version()} with pingouin {v('pingouin')};"
        f" the figure was produced with matplotlib {v('matplotlib')}."
    )


def _script(spec, r, color, in_name, fig_stub):
    x, y = spec["data"]["x"], spec["data"]["y"]
    ap = spec["appearance"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, pingouin as pg
from scipy import stats


{inspect.getsource(lm_band)}

raw = pd.read_csv("{in_name}")
x = pd.to_numeric(raw["{x}"], errors="coerce")
y = pd.to_numeric(raw["{y}"], errors="coerce")
mask = np.isfinite(x) & np.isfinite(y)
x, y = x[mask].to_numpy(), y[mask].to_numpy()
res = pg.corr(x, y, method="{r['method']}"); print(res)
sub = "{r['label']}, {r['estimate_name']} = %.2f, P = %.2g" % (res["r"].iloc[0], res["p_val"].iloc[0])

fig, ax = plt.subplots(figsize=(4.6, 4.2))
ax.scatter(x, y, s={(2 * R_POINT_PT) ** 2:.2f}, color="{color}", alpha=0.8, edgecolors="none")
grid, fit, lo, hi = lm_band(x, y)
ax.fill_between(grid, lo, hi, color="#BFBFBF", alpha=0.4, linewidth=0)
ax.plot(grid, fit, color="#333333", linewidth=1.4)
ax.set_xlabel("{ap.get('x_label') or x}"); ax.set_ylabel("{ap.get('y_label') or y}")
ax.set_title(sub, fontsize=10, color="#4d4d4d")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
