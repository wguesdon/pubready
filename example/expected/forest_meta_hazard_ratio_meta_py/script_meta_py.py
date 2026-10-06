#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import numpy as np, pandas as pd
from statsmodels.stats.meta_analysis import combine_effects

raw = pd.read_csv("input_hazard_ratio.csv")
yi = np.log(raw["hazard_ratio"].to_numpy())
sei = (np.log(raw["ci_high"]) - np.log(raw["ci_low"])) / (2 * 1.959963984540054)
res = combine_effects(yi, np.asarray(sei) ** 2, method_re="dl")
print(res.summary_frame())
print("Q", res.q, "I2", res.i2, "tau2", res.tau2)
fig = res.plot_forest()
fig.savefig("figure_hazard_ratio_meta_py.pdf", bbox_inches="tight")
print("Reproduced figure.")
