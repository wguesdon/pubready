#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The AUC is
# the Mann-Whitney probability that a case outranks a control.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("input_status_marker_a_marker_b.csv")
markers = ['marker_a', 'marker_b']
raw = raw.dropna(subset=["status"] + markers)
is_case = raw["status"].astype(str).str.replace(r"\.0$", "", regex=True) == "disease"
fig, ax = plt.subplots(figsize=(4.8, 5.0))
ax.plot([0, 1], [0, 1], ls="--", color="0.6")
for m in markers:
    cases, controls = raw.loc[is_case, m].to_numpy(), raw.loc[~is_case, m].to_numpy()
    sign = 1 if np.median(controls) <= np.median(cases) else -1
    gt = (sign * cases[:, None] > sign * controls[None, :]).mean()
    eq = (cases[:, None] == controls[None, :]).mean()
    auc = gt + 0.5 * eq
    thr = np.sort(np.unique(sign * raw[m].to_numpy()))
    sens = [(sign * cases >= t).mean() for t in thr]
    fpr = [(sign * controls >= t).mean() for t in thr]
    ax.plot(fpr, sens, label=f"{m}: AUC {auc:.3f}")
    print(m, "AUC", auc)
ax.set_xlabel("1 - Specificity"); ax.set_ylabel("Sensitivity"); ax.set_aspect("equal")
for side in ("top", "right"): ax.spines[side].set_visible(False)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig("figure_status_marker_a_marker_b_roc_py.pdf")
print("Reproduced figure (healthy = control, disease = case).")
