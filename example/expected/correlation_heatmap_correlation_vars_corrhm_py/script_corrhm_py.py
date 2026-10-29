#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
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


raw = pd.read_csv("input_correlation_vars.csv")
num = raw.select_dtypes(include="number")
M = num.corr(method="pearson")
# The leaf order of the figure: complete linkage on the Euclidean distance of
# the matrix rows, and at each merge the branch with the higher mean correlation
# first. On a tie, the branch with the earlier input variable.
mv = M.to_numpy()
row_z = ordered_linkage(mv, mv.mean(axis=1)) if True and len(mv) > 1 else None
col_z = ordered_linkage(mv.T, mv.mean(axis=0)) if True and len(mv) > 1 else None
g = sns.clustermap(M, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                   row_cluster=row_z is not None, col_cluster=col_z is not None,
                   row_linkage=row_z, col_linkage=col_z,
                   annot=True, fmt=".2f", annot_kws={"size": 7},
                   cbar_kws={"label": "correlation"}, linewidths=0.5, linecolor="white")
g.figure.savefig("figure_correlation_vars_corrhm_py.pdf")
print("Reproduced figure.")
