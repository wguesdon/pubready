"""Small shared helpers (mirrors the R util functions)."""
import numpy as np
from scipy import stats


def cap_first(s):
    s = str(s)
    return s[:1].upper() + s[1:] if s else s


def p_stars(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "ns"
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def format_p(p, digits=2):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    if p < 0.001:
        return "< 0.001"
    return f"{p:.{digits}g}"


def shapiro_normal(x):
    """Return (p_value, is_normal) guarded for small/large n."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return (np.nan, None)
    if len(x) > 5000:
        return (np.nan, True)
    try:
        p = float(stats.shapiro(x).pvalue)
    except Exception:
        return (np.nan, None)
    return (p, p > 0.05)


def mann_whitney_r(a, b):
    """Mann-Whitney U test with the p-value rule of R's stats::wilcox.test.

    R computes the exact null distribution when both groups hold fewer than 50
    values and no value is tied, and otherwise the normal approximation with a
    continuity correction. SciPy takes the exact route only for very small
    samples, so the rule is spelled out here to keep the two engines on the same
    p-value.

    Args:
        a: Values of the first group.
        b: Values of the second group.

    Returns:
        tuple: The U statistic of the first group and the two-sided p-value.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    pooled = np.concatenate([a, b])
    exact = len(a) < 50 and len(b) < 50 and len(np.unique(pooled)) == len(pooled)
    if exact:
        res = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    else:
        res = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic",
                                 use_continuity=True)
    return float(res.statistic), float(res.pvalue)


def wilcoxon_signed_rank_r(d):
    """Wilcoxon signed-rank test with the statistic and p-value rule of R.

    R drops the zero differences, reports V, the sum of the ranks of the
    positive differences, and takes the exact route when fewer than 50 non-zero
    differences remain with no tie and no zero. SciPy reports the smaller rank
    sum and picks its own route, so both are restated here.

    Args:
        d: Paired differences, first condition minus second condition.

    Returns:
        tuple: The V statistic and the two-sided p-value.
    """
    d = np.asarray(d, dtype=float)
    nz = d[d != 0]
    if len(nz) == 0:
        return (np.nan, np.nan)
    rk = stats.rankdata(np.abs(nz))
    v = float(rk[nz > 0].sum())
    exact = len(nz) < 50 and len(nz) == len(d) and len(np.unique(np.abs(nz))) == len(nz)
    if exact:
        res = stats.wilcoxon(nz, alternative="two-sided", method="exact")
    else:
        res = stats.wilcoxon(nz, alternative="two-sided", method="approx", correction=True)
    return v, float(res.pvalue)


def rank_biserial_unpaired(u, n1, n2):
    """Rank-biserial correlation from the Mann-Whitney U of the first group.

    Args:
        u: The U statistic of the first group.
        n1: Size of the first group.
        n2: Size of the second group.

    Returns:
        float: The effect size, negative when the first group ranks lower.
    """
    return 2.0 * float(u) / (n1 * n2) - 1.0


def rank_biserial_paired(d):
    """Rank-biserial correlation for paired differences (Kerby's formula).

    Args:
        d: Paired differences, first condition minus second condition.

    Returns:
        float: The effect size, or NaN when every difference is zero.
    """
    d = np.asarray(d, dtype=float)
    nz = d[d != 0]
    if len(nz) == 0:
        return np.nan
    rk = stats.rankdata(np.abs(nz))
    return float((rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum())
