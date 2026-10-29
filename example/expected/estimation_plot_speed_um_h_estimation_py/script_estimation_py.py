#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. It
# recomputes the BCa bootstrap intervals with the same Park-Miller resamples and
# draws a simple version of the estimation plot.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

M, A, SEED, B = 2147483647, 48271, 20261029, 5000


def uniforms(n):
    out, s = np.empty(n), SEED
    for i in range(n):
        s = (s * A) % M
        out[i] = s
    return out / M


def bca(r, g, paired):
    per = len(r) if paired else len(r) + len(g)
    u = uniforms(B * per).reshape(B, per)
    if paired:
        i = np.floor(u * len(r)).astype(int)
        boot = (g[i] - r[i]).sum(axis=1) / len(r)
        dd = g - r
        jack = (dd.sum() - dd) / (len(dd) - 1)
    else:
        ir = np.floor(u[:, :len(r)] * len(r)).astype(int)
        ig = np.floor(u[:, len(r):] * len(g)).astype(int)
        boot = g[ig].sum(axis=1) / len(g) - r[ir].sum(axis=1) / len(r)
        jack = np.concatenate([g.mean() - (r.sum() - r) / (len(r) - 1),
                               (g.sum() - g) / (len(g) - 1) - r.mean()])
    est = g.sum() / len(g) - r.sum() / len(r)
    p0 = np.sum(boot < est - 1e-11 * np.abs(np.concatenate([r, g])).max()) / B
    z0 = stats.norm.ppf(min(max(p0, 0.5 / B), 1 - 0.5 / B))
    dj = jack.mean() - jack
    acc = np.sum(dj ** 3) / (6 * np.sum(dj ** 2) ** 1.5)
    zq = stats.norm.ppf([0.025, 0.975])
    return est, *np.quantile(boot, stats.norm.cdf(z0 + (z0 + zq) / (1 - acc * (z0 + zq))))


d = pd.read_csv("input_speed_um_h.csv")
d = d[d["treatment"].notna() & np.isfinite(pd.to_numeric(d["speed_um_h"], errors="coerce"))]
d["treatment"] = d["treatment"].astype(str)
lv = ['vehicle', 'compound_A', 'compound_B', 'compound_C']
res = []
for g in lv[1:]:
    r = d.loc[d["treatment"] == lv[0], "speed_um_h"].to_numpy(float)
    v = d.loc[d["treatment"] == g, "speed_um_h"].to_numpy(float)
    est, lo, hi = bca(r, v, False)
    p = (stats.ttest_rel(v, r) if False else stats.ttest_ind(v, r, equal_var=False)).pvalue
    res.append((g, est, lo, hi, p))
    print(f"{g}: mean difference {est:.6g} (95% CI {lo:.6g} to {hi:.6g}), P = {p:.4g}")
fig, (a1, a2) = plt.subplots(2, 1, figsize=(5, 6))
for i, g in enumerate(lv):
    vals = d.loc[d["treatment"] == g, "speed_um_h"]
    a1.scatter(np.full(len(vals), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(vals)), vals)
a1.set_xticks(range(len(lv))); a1.set_xticklabels(lv); a1.set_ylabel("Migration speed (um/h)")
for i, (g, est, lo, hi, p) in enumerate(res, start=1):
    a2.errorbar([i], [est], yerr=[[est - lo], [hi - est]], fmt="o", color="black")
a2.axhline(0, ls="--", color="grey"); a2.set_xlim(a1.get_xlim())
a2.set_xticks(range(len(lv))); a2.set_xticklabels([""] + lv[1:])
a2.set_ylabel("Mean difference vs vehicle")
fig.tight_layout(); fig.savefig("figure_speed_um_h_estimation_py.pdf")
print("Reproduced figure.")
