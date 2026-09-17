#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, upsetplot

raw = pd.read_csv("input_set_membership.csv", index_col=0)
set_cols = ['apoptosis', 'inflammation', 'metabolism', 'cell_cycle']
bin_df = raw[set_cols].apply(lambda s: pd.to_numeric(s, errors="coerce") == 1)
bin_df = bin_df[bin_df.any(axis=1)].reset_index(drop=True)
data = upsetplot.from_indicators(set_cols, bin_df)
fig = plt.figure(figsize=(7, 4.5))
upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", show_counts=True).plot(fig=fig)
fig.savefig("figure_set_membership_upset_py.pdf")
print("Reproduced figure.")
