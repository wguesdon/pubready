#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import pandas as pd, seaborn as sns

raw = pd.read_csv("input_correlation_vars.csv")
num = raw.select_dtypes(include="number")
M = num.corr(method="pearson")
g = sns.clustermap(M, cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                   row_cluster=True, col_cluster=True,
                   annot=True, fmt=".2f", annot_kws={"size": 7},
                   cbar_kws={"label": "correlation"}, linewidths=0.5, linecolor="white")
g.figure.savefig("figure_correlation_vars_corrhm_py.pdf")
print("Reproduced figure.")
