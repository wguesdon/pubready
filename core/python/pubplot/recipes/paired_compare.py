"""Recipe: paired before/after comparison. Mirrors the R recipe: each subject
(--id) is measured under two conditions (--x), and a line connects the subject's
two points. The paired t-test or the Wilcoxon signed-rank test is chosen from the
normality of the within-subject differences."""
import numpy as np
import pandas as pd
from scipy import stats as ss

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig, palette
from ..util import cap_first, p_stars, shapiro_normal


def recipe_paired_compare(df, spec):
    x, y, idc = spec["data"]["x"], spec["data"]["y"], spec["data"].get("id")
    if not x or not y:
        raise ValueError("paired_compare needs --x (condition) and --y (value)")
    if not idc:
        raise ValueError("paired_compare needs --id (the subject column)")
    for c in (idc, x, y):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    d = df[[idc, x, y]].copy()
    d[y] = pd.to_numeric(d[y], errors="coerce")
    d = d.dropna()
    d[x] = d[x].astype(str)
    lv = list(pd.unique(d[x]))  # preserve appearance order (pre before post)
    if len(lv) != 2:
        raise ValueError(f"paired_compare needs exactly 2 conditions in '{x}'; found {len(lv)}")

    a = d[d[x] == lv[0]][[idc, y]].drop_duplicates(subset=idc).rename(columns={y: "v1"})
    b = d[d[x] == lv[1]][[idc, y]].drop_duplicates(subset=idc).rename(columns={y: "v2"})
    m = a.merge(b, on=idc)
    dropped = d[idc].nunique() - len(m)
    if len(m) < 2:
        raise ValueError("paired_compare needs at least 2 complete pairs")
    steps = [f"Kept {len(m)} complete pair(s) across the two conditions."]
    if dropped:
        steps.append(f"Dropped {dropped} subject(s) missing one condition.")

    diffs = (m["v2"] - m["v1"]).to_numpy()
    pdiff, dnorm = shapiro_normal(diffs)
    method = (spec["test"].get("method") or "auto").lower()
    if method in ("auto", ""):
        family = "t" if dnorm else "wilcoxon"
    elif method in ("t", "t_test", "paired_t", "student_t"):
        family = "t"
    elif method in ("wilcoxon", "wilcox", "signed_rank", "mwu"):
        family = "wilcoxon"
    else:
        raise ValueError(f"unknown paired test '{method}' (use t or wilcoxon)")

    if family == "t":
        res = ss.ttest_rel(m["v2"], m["v1"])
        stat = float(res.statistic); pval = float(res.pvalue); dof = len(m) - 1
        eff = float(np.mean(diffs) / np.std(diffs, ddof=1)); eff_name = "Cohen's dz"
        label = "paired t-test"
    else:
        nz = diffs[diffs != 0]
        res = ss.wilcoxon(m["v2"], m["v1"])
        stat = float(res.statistic); pval = float(res.pvalue); dof = np.nan
        rk = ss.rankdata(np.abs(nz))
        eff = float((rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()) if len(nz) else np.nan
        eff_name = "rank-biserial r"; label = "Wilcoxon signed-rank test"

    ap = spec["appearance"]
    pal = palette(2, ap.get("palette"))
    w, h = 3.8, 4.2
    fig, ax = new_fig(w, h)
    for _, r in m.iterrows():
        ax.plot([0, 1], [r["v1"], r["v2"]], color="grey", lw=0.6, alpha=0.7, zorder=1)
    ax.scatter([0] * len(m), m["v1"], color=pal[0], s=45, zorder=2, alpha=0.9)
    ax.scatter([1] * len(m), m["v2"], color=pal[1], s=45, zorder=2, alpha=0.9)
    ax.set_xticks([0, 1]); ax.set_xticklabels(lv); ax.set_xlim(-0.4, 1.4)

    allv = np.r_[m["v1"].to_numpy(), m["v2"].to_numpy()]
    yr = float(np.nanmax(allv) - np.nanmin(allv)) or 1.0
    ytop = float(np.nanmax(allv)) + 0.08 * yr
    blabel = (ap.get("bracket") or {}).get("label", "p.signif")
    brk = f"p = {pval:.2g}" if blabel == "p.format" else p_stars(pval)
    ax.plot([0, 0, 1, 1], [ytop, ytop + 0.02 * yr, ytop + 0.02 * yr, ytop], color="black", lw=1.0)
    ax.text(0.5, ytop + 0.03 * yr, brk, ha="center", va="bottom", fontsize=12)

    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_ylabel(ap.get("y_label") or y)
    if ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    stats_df = pd.DataFrame([{
        "outcome": y, "condition1": lv[0], "condition2": lv[1], "n_pairs": len(m),
        "test": label, "statistic": stat, "df": dof, "p_value": pval,
        "mean_diff": float(np.mean(diffs)), "effect_size_name": eff_name,
        "effect_size": eff, "shapiro_p_diff": pdiff,
    }])
    test_meta = {
        "name": label, "statistic": stat, "p_value": pval, "n_pairs": len(m),
        "mean_difference": float(np.mean(diffs)),
        "effect_size": {"name": eff_name, "value": eff},
        "assumptions": {"normality_of_differences": {"test": "shapiro-wilk", "p": pdiff}},
    }
    methods = _methods(label, lv, len(m), eff_name, eff, spec)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": m, "clean_steps": steps, "label": slugify(y),
        "width": w, "height": h,
        "qc": {"quantity": "within-subject differences", "panels": {"differences": diffs}},
        "build_script": lambda in_name, fig_stub: _script(spec, family, lv, in_name, fig_stub),
    }


def _methods(label, lv, n, eff_name, eff, spec):
    from importlib.metadata import version as v
    import platform
    y_lab = spec["appearance"].get("y_label") or spec["data"]["y"]
    return (
        f"{cap_first(y_lab)} was compared between {lv[0]} and {lv[1]} within each subject with a "
        f"{label} (n = {n} pairs). The test was chosen from the Shapiro-Wilk normality of the "
        f"within-subject differences. Effect size is reported as {eff_name} = {eff:.2f}. "
        "Significance was set at P < 0.05. "
        f"Analyses were performed in Python {platform.python_version()} with scipy {v('scipy')}; "
        f"the figure was produced with matplotlib {v('matplotlib')}."
    )


def _script(spec, family, lv, in_name, fig_stub):
    x, y, idc = spec["data"]["x"], spec["data"]["y"], spec["data"]["id"]
    test_line = ('res = ss.ttest_rel(m["v2"], m["v1"])' if family == "t"
                 else 'res = ss.wilcoxon(m["v2"], m["v1"])')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubplot container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats as ss

raw = pd.read_csv("{in_name}")
raw["{y}"] = pd.to_numeric(raw["{y}"], errors="coerce")
raw = raw.dropna(subset=["{idc}", "{x}", "{y}"])
raw["{x}"] = raw["{x}"].astype(str)
lv = {lv}
a = raw[raw["{x}"] == lv[0]][["{idc}", "{y}"]].drop_duplicates("{idc}").rename(columns={{"{y}": "v1"}})
b = raw[raw["{x}"] == lv[1]][["{idc}", "{y}"]].drop_duplicates("{idc}").rename(columns={{"{y}": "v2"}})
m = a.merge(b, on="{idc}")
{test_line}
print(res)
fig, ax = plt.subplots(figsize=(3.8, 4.2))
for _, r in m.iterrows(): ax.plot([0, 1], [r["v1"], r["v2"]], color="grey", lw=0.6, alpha=0.7)
ax.scatter([0]*len(m), m["v1"], s=45); ax.scatter([1]*len(m), m["v2"], s=45)
ax.set_xticks([0, 1]); ax.set_xticklabels(lv); ax.set_xlabel("{x}"); ax.set_ylabel("{y}")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
