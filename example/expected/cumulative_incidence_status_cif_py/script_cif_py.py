#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The
# functions below restate the cinc and crstm routines of the R package cmprsk.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

def cif_curve(time, any_event, cause_event):
    """The Aalen-Johansen CIF of one cause in one group, as cmprsk's ``cinc``.

    Args:
        time: The follow-up times, sorted in increasing order.
        any_event: 1 when the subject failed from any cause, else 0.
        cause_event: 1 when the subject failed from the cause, else 0.

    Returns:
        tuple: The times and the estimates at the step corners. The curve starts
        at (0, 0), each jump time appears twice (before and after the jump), and
        the last corner is the largest follow-up time.
    """
    n = len(time)
    fk = 1.0
    rs = float(n)
    xs, fs = [0.0], [0.0]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and time[j + 1] == time[i]:
            j += 1
        nd1 = int(cause_event[i:j + 1].sum())
        nd = int(any_event[i:j + 1].sum())
        if nd > 0:
            fkn = fk * (rs - nd) / rs
            if nd1 > 0:
                xs += [float(time[i]), float(time[i])]
                fs += [fs[-1], fs[-1] + fk * nd1 / rs]
            fk = fkn
        rs = float(n - (j + 1))
        i = j + 1
    xs.append(float(time[n - 1]))
    fs.append(fs[-1])
    return np.array(xs), np.array(fs)


def gray_test(time, code, group, ng, rho=0.0):
    """Gray's k-sample test of one cause, restated from cmprsk's ``crst``.

    Args:
        time: The follow-up times.
        code: 0 = censored, 1 = failed from the cause, 2 = failed from another cause.
        group: The group of each subject, coded 0 to ng - 1.
        ng: The count of groups.
        rho: The power of the weight function.

    Returns:
        tuple: The scores of the first ng - 1 groups and their variance matrix.
    """
    order = np.argsort(time, kind="stable")
    y = np.asarray(time, dtype=float)[order]
    m = np.asarray(code, dtype=int)[order]
    ig = np.asarray(group, dtype=int)[order]
    n = len(y)
    ng1 = ng - 1
    rs = np.bincount(ig, minlength=ng).astype(float)
    s = np.zeros(ng1)
    v = np.zeros((ng1, ng1))
    f1m, f1 = np.zeros(ng), np.zeros(ng)
    skmm, skm = np.ones(ng), np.ones(ng)
    v3 = np.zeros(ng)
    v2 = np.zeros((ng1, ng))
    c = np.zeros((ng, ng))
    a = np.zeros((ng, ng))
    fm = 0.0
    f = 0.0
    ll = 0
    while True:
        lu = ll
        while lu + 1 < n and not y[lu + 1] > y[ll]:
            lu += 1
        d = np.zeros((3, ng), dtype=int)
        for i in range(ll, lu + 1):
            d[m[i], ig[i]] += 1
        nd1 = int(d[1].sum())
        nd2 = int(d[2].sum())
        if nd1 > 0 or nd2 > 0:
            f = _gray_step(d, nd1, nd2, rs, f1m, f1, skmm, skm, fm, rho,
                           s, v, v3, v2, c, a, ng, ng1)
        if lu >= n - 1:
            break
        for i in range(ll, lu + 1):
            rs[ig[i]] -= 1
        fm = f
        f1m[:] = f1
        skmm[:] = skm
        ll = lu + 1
    for i in range(ng1):
        for j in range(i + 1):
            for k in range(ng):
                v[i, j] += c[i, k] * c[j, k] * v3[k]
                v[i, j] += c[i, k] * v2[j, k]
                v[i, j] += c[j, k] * v2[i, k]
    for i in range(ng1):
        for j in range(i):
            v[j, i] = v[i, j]
    return s, v


def _gray_step(d, nd1, nd2, rs, f1m, f1, skmm, skm, fm, rho, s, v, v3, v2, c, a, ng, ng1):
    """Process one failure time of Gray's test, updating the arrays in place.

    Args:
        d: The counts censored, failed from the cause and failed from another
            cause at this time, by group.
        nd1: The count failed from the cause at this time.
        nd2: The count failed from another cause at this time.
        rs: The risk set size of each group.
        f1m: The left-continuous CIF of each group.
        f1: The right-continuous CIF of each group.
        skmm: The left-continuous Kaplan-Meier estimate of each group.
        skm: The right-continuous Kaplan-Meier estimate of each group.
        fm: The left-continuous pooled CIF.
        rho: The power of the weight function.
        s: The scores.
        v: The variance matrix, lower triangle.
        v3: A variance work array.
        v2: A variance work array.
        c: A work matrix.
        a: A work matrix.
        ng: The count of groups.
        ng1: The count of groups less one.

    Returns:
        float: The right-continuous pooled CIF.
    """
    tr = 0.0
    tq = 0.0
    for i in range(ng):
        if rs[i] <= 0:
            continue
        td = float(d[1, i] + d[2, i])
        skm[i] = skmm[i] * (rs[i] - td) / rs[i]
        f1[i] = f1m[i] + (skmm[i] * d[1, i]) / rs[i]
        tr += rs[i] / skmm[i]
        tq += rs[i] * (1 - f1m[i]) / skmm[i]
    f = fm + nd1 / tr
    fb = (1 - fm) ** rho
    for i in range(ng):
        a[i, i:] = 0.0
        if rs[i] <= 0:
            continue
        t1 = rs[i] / skmm[i]
        a[i, i] = fb * t1 * (1 - t1 / tr)
        if a[i, i] != 0:
            c[i, i] += a[i, i] * nd1 / (tr * (1 - fm))
        for j in range(i + 1, ng):
            if rs[j] <= 0:
                continue
            a[i, j] = -fb * t1 * rs[j] / (skmm[j] * tr)
            if a[i, j] != 0:
                c[i, j] += a[i, j] * nd1 / (tr * (1 - fm))
    for i in range(1, ng):
        for j in range(i):
            a[i, j] = a[j, i]
            c[i, j] = c[j, i]
    for i in range(ng1):
        if rs[i] <= 0:
            continue
        s[i] += fb * (d[1, i] - nd1 * rs[i] * (1 - f1m[i]) / (skmm[i] * tq))
    if nd1 > 0:
        for k in range(ng):
            if rs[k] <= 0:
                continue
            t4 = 1.0
            if skm[k] > 0:
                t4 = 1 - (1 - f) / skm[k]
            t5 = 1.0
            if nd1 > 1:
                t5 = 1 - (nd1 - 1) / (tr * skmm[k] - 1)
            t3 = t5 * skmm[k] * nd1 / (tr * rs[k])
            v3[k] += t4 * t4 * t3
            for i in range(ng1):
                t1 = a[i, k] - t4 * c[i, k]
                v2[i, k] += t1 * t4 * t3
                for j in range(i + 1):
                    t2 = a[j, k] - t4 * c[j, k]
                    v[i, j] += t1 * t2 * t3
    if nd2 > 0:
        for k in range(ng):
            if skm[k] <= 0 or d[2, k] <= 0:
                continue
            t4 = (1 - f) / skm[k]
            t5 = 1.0
            if d[2, k] > 1:
                t5 = 1 - (d[2, k] - 1.0) / (rs[k] - 1.0)
            t6 = rs[k]
            t3 = t5 * ((skmm[k] ** 2) * d[2, k]) / (t6 ** 2)
            v3[k] += t4 * t4 * t3
            for i in range(ng1):
                t1 = t4 * c[i, k]
                v2[i, k] -= t1 * t4 * t3
                for j in range(i + 1):
                    t2 = t4 * c[j, k]
                    v[i, j] += t1 * t2 * t3
    return f


def gray_statistic(s, v):
    """The chi-square statistic s' inv(V) s, or NaN when V is singular.

    cmprsk checks the rank of V with a QR decomposition; the rank here comes from
    the singular values with the same relative tolerance of 1e-7.

    Args:
        s: The scores.
        v: The variance matrix.

    Returns:
        float: The statistic.
    """
    sv = np.linalg.svd(v, compute_uv=False)
    if len(sv) == 0 or sv[-1] <= 1e-7 * sv[0]:
        return np.nan
    return float(s @ np.linalg.solve(v, s))


raw = pd.read_csv("input_status.csv")
st = raw["status"].astype(str).str.replace(r"\.0$", "", regex=True).str.strip()
d = pd.DataFrame({"time": pd.to_numeric(raw["months"], errors="coerce"), "status": st,
                   "group": raw["arm"].astype(str).str.strip()}).dropna()
d = d[d["time"] >= 0].sort_values("time", kind="stable").reset_index(drop=True)
lv = sorted(d["group"].unique())
causes = sorted(d.loc[d["status"] != "0", "status"].unique())
fig, axes = plt.subplots(1, len(causes), figsize=(9.35, 5.60), squeeze=False)
for ax, cause in zip(axes[0], causes):
    ev, cs = (d["status"] != "0").to_numpy(), (d["status"] == cause).to_numpy()
    for g in lv:
        sel = (d["group"] == g).to_numpy()
        xs, fs = cif_curve(d["time"].to_numpy()[sel], ev[sel].astype(int), cs[sel].astype(int))
        ax.plot(xs, fs, label=g)
        print(cause, g, [float(fs[np.nonzero(xs <= t)[0][-1]]) if t <= xs[-1] else None
                         for t in [20.0, 30.0, 50.0]])
    if len(lv) > 1:
        code = np.where(~ev, 0, np.where(cs, 1, 2))
        s, v = gray_test(d["time"].to_numpy(float), code, d["group"].map(lv.index).to_numpy(), len(lv))
        stat = gray_statistic(s, v)
        print("Gray's test", cause, stat, stats.chi2.sf(stat, len(lv) - 1))
    ax.set_title("status = " + cause); ax.set_xlabel("months"); ax.set_ylabel("Cumulative incidence")
    ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figure_status_cif_py.pdf")
print("Reproduced figure.")
