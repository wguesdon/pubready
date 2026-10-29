#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container. The helper
# functions below group the subjects as ResourceSelection::hoslem.test does.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

def quantile7(x, probs):
    """Type-7 sample quantiles with the arithmetic of R's quantile.default.

    Args:
        x: The values.
        probs: The probabilities.

    Returns:
        numpy.ndarray: The quantiles.
    """
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    index = 1 + max(n - 1, 0) * np.asarray(probs, dtype=float)
    lo = np.floor(index).astype(int)
    hi = np.ceil(index).astype(int)
    qs = x[lo - 1].copy()
    sel = (index > lo) & (x[hi - 1] != qs)
    h = (index - lo)[sel]
    qs[sel] = (1 - h) * qs[sel] + h * x[hi[sel] - 1]
    return qs


def hl_groups(p, g):
    """The Hosmer-Lemeshow group of each subject, as hoslem.test forms it.

    Args:
        p: The predicted probabilities.
        g: The number of groups requested.

    Returns:
        numpy.ndarray: The group of each subject, from 1.
    """
    probs = np.minimum(0.0 + np.arange(g + 1) * (1 / g), 1.0)
    qq = np.unique(quantile7(p, probs))
    if len(qq) < 2:
        return np.ones(len(p), dtype=int)
    code = np.searchsorted(qq, p, side="left")
    code[code == 0] = 1
    return code


def group_table(p, y, g):
    """Per-group counts for the test and the figure.

    Args:
        p: The predicted probabilities.
        y: The outcome, 1 for an event.
        g: The number of groups requested.

    Returns:
        pandas.DataFrame: For each non-empty group: the size, the events, the
        expected events and non-events, the mean prediction, the observed
        proportion and its 95% Wilson interval.
    """
    grp = hl_groups(p, g)
    ids = np.unique(grp)
    n = np.array([np.sum(grp == i) for i in ids], dtype=float)
    o1 = np.array([np.sum(y[grp == i]) for i in ids], dtype=float)
    e1 = np.array([np.sum(p[grp == i]) for i in ids])
    e0 = np.array([np.sum(1 - p[grp == i]) for i in ids])
    obs = o1 / n
    z = stats.norm.ppf(0.975)
    centre = (obs + z ** 2 / (2 * n)) / (1 + z ** 2 / n)
    half = z * np.sqrt(obs * (1 - obs) / n + z ** 2 / (4 * n ** 2)) / (1 + z ** 2 / n)
    return pd.DataFrame({"group": ids, "n": n, "events": o1, "expected": e1,
                         "expected_non": e0, "mean_pred": e1 / n, "observed": obs,
                         "ci_low": centre - half, "ci_high": centre + half})


def hosmer_lemeshow(tab):
    """The Hosmer-Lemeshow statistic, its df and its p-value.

    Args:
        tab: The per-group table.

    Returns:
        dict: The statistic, df = groups - 2, the p-value and the group count.
    """
    stat = float(np.sum((tab["events"] - tab["expected"]) ** 2 / tab["expected"]) +
                 np.sum(((tab["n"] - tab["events"]) - tab["expected_non"]) ** 2 /
                        tab["expected_non"]))
    df = len(tab) - 2
    p = float(stats.chi2.sf(stat, df)) if df > 0 and np.isfinite(stat) else np.nan
    return {"statistic": stat, "df": df, "p_value": p, "groups": len(tab)}


raw = pd.read_csv("input_event_1y_risk_clinical_risk_extended.csv")
models = ['risk_clinical', 'risk_extended']
raw = raw.dropna(subset=["event_1y"] + models)
yb = (raw["event_1y"].astype(str).str.replace(r"\.0$", "", regex=True) == "yes").to_numpy().astype(int)
fig, ax = plt.subplots(figsize=(5.20, 8.56))
ax.plot([0, 1], [0, 1], ls="--", color="0.6")
for m in models:
    p = raw[m].to_numpy(float)
    lp = np.log(np.clip(p, 1e-8, 1 - 1e-8) / (1 - np.clip(p, 1e-8, 1 - 1e-8)))
    ones = np.ones((len(p), 1))
    print(m, "Brier", np.mean((p - yb) ** 2))
    print(sm.GLM(yb, ones, family=sm.families.Binomial(), offset=lp).fit().summary2().tables[1])
    print(sm.GLM(yb, np.column_stack([ones, lp]), family=sm.families.Binomial()).fit().summary2().tables[1])
    print("Hosmer-Lemeshow", hosmer_lemeshow(group_table(p, yb, 10)))
    s = lowess(yb, p, frac=2 / 3, it=0, delta=0.01 * (p.max() - p.min()))
    ax.plot(s[:, 0], s[:, 1], label=m)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax.set_xlabel("Predicted probability"); ax.set_ylabel("Observed proportion")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("figure_event_1y_risk_clinical_risk_extended_cal_py.pdf")
print("Reproduced figure.")
