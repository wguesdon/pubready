"""Recipe: clustered correlation heatmap among the numeric columns of a table.
Mirrors the R recipe: correlation matrix drawn as a seaborn clustermap with
dendrograms, each cell showing the coefficient and BH-adjusted stars. Pearson by
default; --test spearman|kendall."""
import inspect
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as ss
from scipy.cluster import hierarchy
from statsmodels.stats.multitest import multipletests

from ..bundle import slugify
from ..theme import diverging_cmap
from ..util import cap_first
from .heatmap import ordered_linkage


def recipe_correlation_heatmap(df, spec):
    num = df.select_dtypes(include="number")
    if num.shape[1] < 2:
        raise ValueError("correlation_heatmap needs at least 2 numeric columns")
    n0 = num.shape[1]
    keep = [c for c in num.columns if num[c].notna().sum() >= 3 and num[c].std() > 0]
    num = num[keep]
    if num.shape[1] < 2:
        raise ValueError("correlation_heatmap needs at least 2 numeric columns with variance")
    clean_steps = ([f"Dropped {n0 - num.shape[1]} non-numeric or zero-variance column(s) "
                    f"({n0} -> {num.shape[1]} variables)."] if num.shape[1] < n0
                   else ["No cleaning applied; numeric columns used as-is."])

    method = (spec["test"].get("method") or "auto").lower()
    if method in ("auto", ""):
        method = "pearson"
    if method not in ("pearson", "spearman", "kendall"):
        raise ValueError(f"unknown correlation method '{method}' "
                         "(use pearson, spearman, or kendall)")

    M = num.corr(method=method)
    varlist = list(M.columns)
    k = len(varlist)
    fn = {"pearson": ss.pearsonr, "spearman": ss.spearmanr, "kendall": ss.kendalltau}[method]

    P = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            if i == j:
                continue
            a, b = num[varlist[i]], num[varlist[j]]
            m2 = a.notna() & b.notna()
            try:
                P[i, j] = float(fn(a[m2], b[m2])[1])
            except Exception:
                P[i, j] = np.nan

    iu = np.triu_indices(k, 1)
    p_ut = P[iu]
    padj_ut = multipletests(p_ut, method="fdr_bh")[1]
    Padj = np.zeros((k, k))
    Padj[iu] = padj_ut
    Padj = Padj + Padj.T

    def star(p):
        if np.isnan(p):
            return ""
        return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

    labels = np.empty((k, k), dtype=object)
    for i in range(k):
        for j in range(k):
            labels[i, j] = "1" if i == j else f"{M.iloc[i, j]:.2f}{star(Padj[i, j])}"

    cluster = spec["appearance"].get("cluster") or "both"
    row_c = cluster in ("both", "rows")
    col_c = cluster in ("both", "columns")
    w = max(4.8, 1.4 + 0.55 * k)
    h = max(4.5, 1.2 + 0.55 * k)
    # The leaf rule of heatmap (ordered_linkage): complete linkage on the Euclidean
    # distance of the matrix rows, and at each merge the branch with the higher
    # mean correlation first; on a tie, the branch with the earlier input variable.
    mv = M.to_numpy()
    row_z = ordered_linkage(mv, mv.mean(axis=1)) if row_c and k > 1 else None
    col_z = ordered_linkage(mv.T, mv.mean(axis=0)) if col_c and k > 1 else None
    row_order = [varlist[i] for i in hierarchy.leaves_list(row_z)] if row_z is not None else varlist
    col_order = [varlist[i] for i in hierarchy.leaves_list(col_z)] if col_z is not None else varlist

    grid = sns.clustermap(
        M, cmap=diverging_cmap(), center=0, vmin=-1, vmax=1, method="complete",
        row_cluster=row_z is not None, col_cluster=col_z is not None,
        row_linkage=row_z, col_linkage=col_z,
        annot=labels, fmt="", annot_kws={"size": 7},
        figsize=(w, h), cbar_kws={"label": "correlation"},
        linewidths=0.5, linecolor="white")
    fig = grid.figure

    stats_df = pd.DataFrame([{
        "var1": varlist[a], "var2": varlist[b], "r": float(M.iloc[a, b]),
        "p_value": float(P[a, b]), "p_adj": float(Padj[a, b]),
        "significant": bool(Padj[a, b] < 0.05),
    } for a, b in zip(*iu)])

    n_sig = int((padj_ut < 0.05).sum())
    test_meta = {
        "name": f"Clustered {method} correlation heatmap (seaborn.clustermap)",
        "method": method, "clustering": cluster,
        "n_variables": k, "n_pairs": len(p_ut), "n_significant": n_sig,
        "row_order": [str(v) for v in row_order], "column_order": [str(v) for v in col_order],
    }
    methods = _methods(method, cluster, k, n_sig)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": num.reset_index(drop=True), "clean_steps": clean_steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(in_name, fig_stub, method, row_c, col_c),
    }


def _methods(method, cluster, k, n_sig):
    from importlib.metadata import version as v
    import platform
    from math import comb
    clust_txt = {"both": "Variables were hierarchically clustered on both axes.",
                 "rows": "Rows were hierarchically clustered.",
                 "columns": "Columns were hierarchically clustered.",
                 "none": "No clustering was applied."}.get(cluster, "")
    return (
        f"Pairwise {cap_first(method)} correlations among {k} variables were computed and drawn as a "
        f"clustered heatmap with seaborn. {clust_txt} Each cell shows the coefficient with "
        f"Benjamini-Hochberg adjusted significance stars ({n_sig} of the {comb(k, 2)} pairs "
        "were significant at adjusted P < 0.05)."
        f" Rendered in Python {platform.python_version()} with seaborn {v('seaborn')}."
    )


def _script(in_name, fig_stub, method, row_c, col_c):
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd, seaborn as sns
from scipy.cluster import hierarchy


{inspect.getsource(ordered_linkage)}

raw = pd.read_csv("{in_name}")
num = raw.select_dtypes(include="number")
M = num.corr(method="{method}")
# The leaf order of the figure: complete linkage on the Euclidean distance of
# the matrix rows, and at each merge the branch with the higher mean correlation
# first. On a tie, the branch with the earlier input variable.
mv = M.to_numpy()
row_z = ordered_linkage(mv, mv.mean(axis=1)) if {row_c} and len(mv) > 1 else None
col_z = ordered_linkage(mv.T, mv.mean(axis=0)) if {col_c} and len(mv) > 1 else None
g = sns.clustermap(M, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                   row_cluster=row_z is not None, col_cluster=col_z is not None,
                   row_linkage=row_z, col_linkage=col_z,
                   annot=True, fmt=".2f", annot_kws={{"size": 7}},
                   cbar_kws={{"label": "correlation"}}, linewidths=0.5, linecolor="white")
g.figure.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
