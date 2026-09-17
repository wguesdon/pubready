"""Recipe: clustered heatmap of a features-by-samples matrix. Tries
PyComplexHeatmap (to mirror the R ComplexHeatmap output); falls back to a seaborn
clustermap if PyComplexHeatmap raises. Row z-score and clustering by default;
a per-feature Welch t-test (BH) is computed when the annotation has two groups."""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..bundle import slugify
from ..io import read_tidy
from ..theme import PALETTE


def recipe_heatmap(df, spec):
    feat = df.iloc[:, 0].astype(str)
    mat = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    mat.index = list(feat)
    n0 = len(mat)
    keep = mat.notna().all(axis=1) & (mat.std(axis=1) > 0)
    mat = mat[keep]
    clean_steps = (f"Dropped {n0 - len(mat)} feature(s) with missing values or zero variance "
                   f"({n0} -> {len(mat)})." if len(mat) < n0 else "No cleaning applied; input used as-is.")
    df_used = df[df.iloc[:, 0].astype(str).isin(mat.index)].reset_index(drop=True)

    ap = spec["appearance"]
    scale_mode = ap.get("scale") or "row"
    cluster = ap.get("cluster") or "both"
    if scale_mode == "row":
        disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)
    elif scale_mode == "column":
        disp = (mat - mat.mean()) / mat.std()
    else:
        disp = mat.copy()
    disp = disp.fillna(0.0)
    legend_name = "value" if scale_mode == "none" else "z-score"
    cluster_rows = cluster in ("both", "rows")
    cluster_cols = cluster in ("both", "columns")

    ann, group_col = None, None
    if spec["data"].get("annotation"):
        araw = read_tidy(spec["data"]["annotation"])
        araw = araw.set_index(araw.columns[0])
        ann = araw.reindex(mat.columns)
        ann = ann.astype(str)
        group_col = ann.columns[0]

    # Per-feature differential test when the primary annotation has two groups.
    if ann is not None and ann[group_col].nunique() == 2:
        g = ann[group_col]
        lv = sorted(g.unique())
        ca = [c for c in mat.columns if g.loc[c] == lv[0]]
        cb = [c for c in mat.columns if g.loc[c] == lv[1]]
        pvals = mat.apply(lambda r: stats.ttest_ind(r[ca], r[cb], equal_var=False).pvalue, axis=1).to_numpy()
        padj = multipletests(pvals, method="fdr_bh")[1]
        stats_df = pd.DataFrame({
            "feature": mat.index, f"mean_{lv[0]}": mat[ca].mean(axis=1).to_numpy(),
            f"mean_{lv[1]}": mat[cb].mean(axis=1).to_numpy(),
            "p_value": pvals, "p_adj": padj, "significant": padj < 0.05,
        })
        test_label = f"per-feature Welch t-test ({lv[0]} vs {lv[1]}), BH-adjusted"
        n_sig = int((padj < 0.05).sum())
    else:
        stats_df = pd.DataFrame({"feature": mat.index, "mean": mat.mean(axis=1).to_numpy(),
                                 "sd": mat.std(axis=1).to_numpy()})
        test_label, n_sig = None, None

    w = max(5.5, 2.2 + 0.24 * mat.shape[1])
    h = max(5.0, 1.5 + 0.16 * mat.shape[0])

    used_pch = True
    try:
        import PyComplexHeatmap as pch
        col_ha = None
        if ann is not None:
            colors = {lvl: PALETTE[i % len(PALETTE)]
                      for i, lvl in enumerate(sorted(ann[group_col].unique()))}
            col_ha = pch.HeatmapAnnotation(
                **{group_col: pch.anno_simple(ann[group_col], colors=colors, add_text=False)},
                axis=1, verbose=0)
        plt.figure(figsize=(w, h))
        pch.ClusterMapPlotter(
            data=disp, top_annotation=col_ha,
            row_cluster=cluster_rows, col_cluster=cluster_cols,
            cmap="RdBu_r", label=legend_name,
            show_rownames=True, show_colnames=True,
            row_names_side="right", col_names_side="bottom", verbose=0)
        fig = plt.gcf()
    except Exception:
        plt.close("all")
        used_pch = False
        col_colors = None
        if ann is not None:
            colors = {lvl: PALETTE[i % len(PALETTE)]
                      for i, lvl in enumerate(sorted(ann[group_col].unique()))}
            col_colors = ann[group_col].map(colors)
        grid = sns.clustermap(
            disp, cmap="RdBu_r", center=0, row_cluster=cluster_rows, col_cluster=cluster_cols,
            col_colors=col_colors, figsize=(w, h), cbar_kws={"label": legend_name},
            xticklabels=True, yticklabels=True)
        fig = grid.figure

    test_meta = {
        "name": "Clustered heatmap", "engine_lib": "PyComplexHeatmap" if used_pch else "seaborn.clustermap",
        "scale": scale_mode, "clustering": cluster,
        "n_features": int(mat.shape[0]), "n_samples": int(mat.shape[1]),
        "per_feature_test": test_label, "n_significant": n_sig,
    }
    methods = _methods(scale_mode, cluster, mat.shape[0], mat.shape[1], test_label, n_sig, used_pch)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": clean_steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "extra_inputs": [spec["data"].get("annotation")],
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, scale_mode, cluster_rows, cluster_cols, legend_name, in_name, fig_stub),
    }


def _methods(scale_mode, cluster, n_feat, n_samp, test_label, n_sig, used_pch):
    from importlib.metadata import version as v
    import platform
    lib = "PyComplexHeatmap" if used_pch else "seaborn.clustermap"
    lib_ver = v("PyComplexHeatmap") if used_pch else v("seaborn")
    scale_txt = {"row": "Values were z-scored across samples within each feature (row scaling).",
                 "column": "Values were z-scored within each sample (column scaling).",
                 "none": "Raw values were shown without scaling."}.get(scale_mode)
    clust_txt = {"both": "Rows and columns were hierarchically clustered.",
                 "rows": "Rows were hierarchically clustered.",
                 "columns": "Columns were hierarchically clustered.",
                 "none": "No clustering was applied."}.get(cluster)
    test_txt = (f" A {test_label} was applied per feature; {n_sig} feature(s) were significant "
                "(adjusted P < 0.05)." if test_label else "")
    return (
        f"A heatmap of {n_feat} features across {n_samp} samples was drawn with {lib}. "
        + scale_txt + " " + clust_txt + test_txt +
        f" Rendered in Python {platform.python_version()} with {lib} {lib_ver}."
    )


def _script(spec, scale_mode, cluster_rows, cluster_cols, legend_name, in_name, fig_stub):
    ann = spec["data"].get("annotation")
    ann_name = os.path.basename(ann) if ann else None
    disp_line = {"row": "disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0).fillna(0)",
                 "column": "disp = ((mat - mat.mean()) / mat.std()).fillna(0)",
                 "none": "disp = mat"}.get(scale_mode)
    ann_lines = (f'ann = pd.read_csv("{ann_name}", index_col=0).astype(str).reindex(mat.columns)\n'
                 f'gc = ann.columns[0]\n'
                 f'colors = {{lvl: c for lvl, c in zip(sorted(ann[gc].unique()), ["#3B6DB3", "#C1432B", "#2E8B57"])}}\n'
                 f'col_colors = ann[gc].map(colors)') if ann_name else "col_colors = None"
    return f'''#!/usr/bin/env python3
# Standalone reproduction (seaborn clustermap). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import pandas as pd, seaborn as sns

raw = pd.read_csv("{in_name}", index_col=0)
mat = raw.apply(pd.to_numeric, errors="coerce")
mat = mat[mat.notna().all(axis=1) & (mat.std(axis=1) > 0)]
{disp_line}
{ann_lines}
g = sns.clustermap(disp, cmap="RdBu_r", center=0, row_cluster={cluster_rows}, col_cluster={cluster_cols},
                   col_colors=col_colors, cbar_kws={{"label": "{legend_name}"}},
                   xticklabels=True, yticklabels=True)
g.figure.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
