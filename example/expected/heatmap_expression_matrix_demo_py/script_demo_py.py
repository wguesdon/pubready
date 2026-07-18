#!/usr/bin/env python3
# Standalone reproduction (seaborn clustermap). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import pandas as pd, seaborn as sns

raw = pd.read_csv("input_expression_matrix.csv", index_col=0)
mat = raw.apply(pd.to_numeric, errors="coerce")
mat = mat[mat.notna().all(axis=1) & (mat.std(axis=1) > 0)]
disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0).fillna(0)
ann = pd.read_csv("expression_annotation.csv", index_col=0).astype(str).reindex(mat.columns)
gc = ann.columns[0]
colors = {lvl: c for lvl, c in zip(sorted(ann[gc].unique()), ["#3B6DB3", "#C1432B", "#2E8B57"])}
col_colors = ann[gc].map(colors)
g = sns.clustermap(disp, cmap="RdBu_r", center=0, row_cluster=True, col_cluster=True,
                   col_colors=col_colors, cbar_kws={"label": "z-score"},
                   xticklabels=True, yticklabels=True)
g.figure.savefig("figure_expression_matrix_demo_py.pdf")
print("Reproduced figure.")
