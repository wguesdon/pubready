#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_set_membership.csv", index_col=0)
set_cols = ['apoptosis', 'inflammation', 'metabolism', 'cell_cycle']
# 1, TRUE, T, yes or y is a member; any other code or an empty cell is not.
codes = raw[set_cols].astype(str)
for c in set_cols:
    codes[c] = codes[c].str.strip().str.lower().str.replace(r"\.0$", "", regex=True)
bin_df = codes.isin(["1", "true", "t", "yes", "y"])
bin_df = bin_df[bin_df.any(axis=1)].reset_index(drop=True)
combo = bin_df.groupby(set_cols, sort=False).size().reset_index(name="size")
combo["sets"] = [" & ".join(c for c in set_cols if r[c]) for _, r in combo.iterrows()]
combo["degree"] = combo[set_cols].sum(axis=1)
combo = combo.sort_values(["size", "degree", "sets"], ascending=[False, True, True], kind="stable")
shown = combo.head(15).reset_index(drop=True)
totals = bin_df.sum().sort_values(ascending=False, kind="stable")
order = list(totals.index); k, n = len(order), len(shown)

fig = plt.figure(figsize=(max(6.0, 0.32 * n + 3.5), max(3.8, 0.3 * k + 2.9)))
gs = fig.add_gridspec(2, 2, width_ratios=[max(n, 1), 3], height_ratios=[2, max(k, 1) * 0.6],
                      wspace=0.05, hspace=0.05)
ax_t, ax_m, ax_s = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
for j in range(n):
    member = shown.loc[j, order].to_numpy(dtype=bool)
    ys = np.arange(k)[::-1]
    ax_m.scatter([j] * k, ys, s=60, color=np.where(member, "black", "#CCCCCC"), lw=0)
    if member.sum() > 1:
        ax_m.vlines(j, ys[member].min(), ys[member].max(), color="black", lw=1.6)
ax_m.set_yticks(np.arange(k)[::-1]); ax_m.set_yticklabels(order); ax_m.set_xticks([])
ax_m.set_xlim(-0.5, n - 0.5); ax_m.set_ylim(-0.5, k - 0.5)
ax_t.bar(range(n), shown["size"], width=0.7, color="black"); ax_t.set_xlim(-0.5, n - 0.5)
for j, s in enumerate(shown["size"]):
    ax_t.text(j, s, str(s), ha="center", va="bottom", fontsize=8)
ax_t.set_xticks([]); ax_t.set_ylabel("Intersection size")
ax_s.barh(np.arange(k)[::-1], totals[order], height=0.6, color="black")
for y, v in zip(np.arange(k)[::-1], totals[order]):
    ax_s.text(v, y, " %d" % v, va="center", fontsize=8)
ax_s.set_yticks([]); ax_s.set_ylim(-0.5, k - 0.5); ax_s.set_xlabel("Set size")
fig.savefig("figure_set_membership_upset_py.pdf", bbox_inches="tight")
print("Reproduced figure.")
