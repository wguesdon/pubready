#!/usr/bin/env python3
# Standalone reproduction (seaborn clustermap). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import numpy as np
import pandas as pd, seaborn as sns
from scipy.cluster import hierarchy


def ordered_linkage(values, weights):
    """Build the tree of one heatmap axis with the leaf order of the R engine.

    The tree is complete linkage on the Euclidean distance of the rows of values.
    At each merge the branch with the higher mean weight comes first, and on a
    tie the branch that holds the earlier input index. scipy draws the first
    column of a linkage row on the left, so the swap sets the drawn order, as
    .hm_order_tree() does on the merge matrix in R.

    Args:
        values: The drawn values, one row per leaf.
        weights: The weight of each leaf, the mean of its unscaled values.

    Returns:
        numpy.ndarray: The linkage matrix.
    """
    z = hierarchy.linkage(np.asarray(values, dtype=float), method="complete", metric="euclidean")
    n = len(weights)
    sum_w = np.concatenate([np.asarray(weights, dtype=float), np.zeros(n - 1)])
    size = np.concatenate([np.ones(n), np.zeros(n - 1)])
    first = np.concatenate([np.arange(n, dtype=float), np.zeros(n - 1)])
    for i in range(n - 1):
        a, b = int(z[i, 0]), int(z[i, 1])
        ma, mb = sum_w[a] / size[a], sum_w[b] / size[b]
        if mb > ma or (mb == ma and first[b] < first[a]):
            z[i, 0], z[i, 1] = b, a
        sum_w[n + i] = sum_w[a] + sum_w[b]
        size[n + i] = size[a] + size[b]
        first[n + i] = min(first[a], first[b])
    return z


raw = pd.read_csv("input_expression_matrix.csv", index_col=0)
mat = raw.apply(pd.to_numeric, errors="coerce")
mat = mat[mat.notna().all(axis=1) & (mat.std(axis=1) > 0)]
disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0).fillna(0)
ann = pd.read_csv("expression_annotation.csv", index_col=0).astype(str).reindex(mat.columns)
gc = ann.columns[0]
colors = {lvl: c for lvl, c in zip(sorted(ann[gc].unique()), ["#3B6DB3", "#C1432B", "#2E8B57"])}
col_colors = ann[gc].map(colors)
row_z = ordered_linkage(disp.to_numpy(), mat.mean(axis=1).to_numpy()) if len(disp) > 1 else None
col_z = ordered_linkage(disp.to_numpy().T, mat.mean(axis=0).to_numpy()) if disp.shape[1] > 1 else None
g = sns.clustermap(disp, cmap="RdBu_r", center=0, row_cluster=row_z is not None,
                   col_cluster=col_z is not None, row_linkage=row_z, col_linkage=col_z,
                   col_colors=col_colors, cbar_kws={"label": "z-score"},
                   xticklabels=True, yticklabels=True)
g.figure.savefig("figure_expression_matrix_demo_py.pdf")
print("Reproduced figure.")
