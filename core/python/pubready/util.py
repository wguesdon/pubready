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


def _cld_absorb(m):
    """Drop the redundant columns of a letter matrix (the absorb step).

    A column goes when it is a proper subset of another column, when it equals
    an earlier column, or when it is empty, as .cld_absorb() of the R engine.

    Args:
        m: A boolean matrix, one row per group and one column per letter.

    Returns:
        numpy.ndarray: The matrix without the redundant columns.
    """
    k = m.shape[1]
    keep = np.ones(k, dtype=bool)
    for a in range(k):
        if not m[:, a].any():
            keep[a] = False
            continue
        for b in range(k):
            if a == b or not m[m[:, a], b].all():
                continue
            if m[~m[:, a], b].any() or b < a:
                keep[a] = False
                break
    return m[:, keep]


def cld_letters(levels, group1, group2, p, alpha=0.05):
    """Compact letter display by the insert-and-absorb algorithm of Piepho (2004).

    This is the algorithm of multcompView::multcompLetters, with the same steps
    as cld_letters() of the R engine. Each pair that differs (p < alpha) splits
    every letter column that holds both groups into two, one without each group;
    then the redundant columns are absorbed. The pairs are taken in the order of
    the levels, the columns are sorted so that "a" starts at the first group, and
    the letters run a-z then A-Z. Groups that share a letter do not differ.

    Args:
        levels: The groups, in the order of the figure.
        group1: The first group of each pair.
        group2: The second group of each pair.
        p: The adjusted P value of each pair.
        alpha: The significance level.

    Returns:
        dict: The letters of each group, in the order of ``levels``.
    """
    import string
    levels = [str(g) for g in levels]
    pos = {g: i for i, g in enumerate(levels)}
    m = np.ones((len(levels), 1), dtype=bool)
    pairs = []
    for a, b, pv in zip(group1, group2, p):
        a, b = str(a), str(b)
        if pv is None or not np.isfinite(pv) or not pv < alpha or a not in pos or b not in pos:
            continue
        pairs.append((min(pos[a], pos[b]), max(pos[a], pos[b])))
    for i, j in sorted(pairs):
        both = np.where(m[i, :] & m[j, :])[0]
        if len(both) == 0:
            continue
        dup = m[:, both].copy()
        dup[i, :] = False
        m[j, both] = False
        m = _cld_absorb(np.hstack([m, dup]))
    keys = ["".join("0" if v else "1" for v in m[:, c]) for c in range(m.shape[1])]
    m = m[:, sorted(range(m.shape[1]), key=keys.__getitem__)]
    pool = list(string.ascii_lowercase + string.ascii_uppercase)
    if m.shape[1] > len(pool):
        pool += [f"z{i + 1}" for i in range(m.shape[1] - len(pool))]
    return {g: "".join(pool[c] for c in np.where(m[i, :])[0]) for i, g in enumerate(levels)}


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


def cohens_d_r(a, b, paired=False):
    """Signed Cohen's d with the formulas of rstatix::cohens_d.

    pingouin reports the magnitude only, so the sign that says which group is
    higher is restored here. Two independent groups use the pooled standard
    deviation; paired data use dz, the mean difference over its standard
    deviation.

    Args:
        a: Values of the first group.
        b: Values of the second group, paired with ``a`` when ``paired``.
        paired: Whether the two groups are paired.

    Returns:
        float: The effect size, negative when the first group is lower.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if paired:
        d = a - b
        return float(d.mean() / d.std(ddof=1))
    na, nb = len(a), len(b)
    pooled = np.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    return float((a.mean() - b.mean()) / pooled)


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


def p_phrase(p, digits=2):
    """Write a P value as a phrase for a methods sentence.

    Args:
        p: The P value.
        digits: The count of significant digits.

    Returns:
        str: "P < 0.001" below 0.001, else "P = " and the value. A table cell
        keeps the bare value of format_p.
    """
    text = format_p(p, digits)
    return f"P {text}" if text.startswith("<") else f"P = {text}"


def sig_text(x, digits=3):
    """Write a number with a fixed count of significant digits.

    Unlike "%.3g" it keeps the trailing zeros (3.70, not 3.7), and it drops the
    trailing point that "%#g" leaves on a whole number, as the R engine does.

    Args:
        x: The number.
        digits: The count of significant digits.

    Returns:
        str: The text, or "NA" for a missing value.
    """
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    text = f"{float(x):#.{digits}g}"
    return text[:-1] if text.endswith(".") else text


def df_text(d):
    """Write a degree of freedom for a methods sentence.

    Args:
        d: The degree of freedom.

    Returns:
        str: A whole number as an integer, and a fractional df (Welch,
        Greenhouse-Geisser) with two decimals.
    """
    d = float(d)
    return f"{int(round(d))}" if abs(d - round(d)) < 1e-9 else f"{d:.2f}"


def _level_label(v):
    """Format one value of a grouping column as level_labels() does.

    Args:
        v: One value.

    Returns:
        str or float: The label, or NaN for a missing value.
    """
    if v is None or (isinstance(v, (float, np.floating)) and np.isnan(v)):
        return np.nan
    if isinstance(v, (bool, np.bool_)):
        return "TRUE" if v else "FALSE"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if np.isfinite(v) and v == round(v) and abs(v) < 1e15:
            return str(int(v))
        return f"{v:.15g}"
    return str(v)


def level_labels(values):
    """Text labels for the values of a grouping column.

    A whole number prints with no decimal point and no exponent, so a dose of
    80 reads "80" in both engines, also when pandas reads the column as float
    because it holds a missing value. Other numbers print with 15 significant
    digits. The R engine has the same rule.

    Args:
        values: A pandas Series.

    Returns:
        pandas.Series: The labels as text, with each missing value kept as NaN.
    """
    return values.map(_level_label)


def _level_number(s):
    """Read one label as a number, as R as.numeric() does.

    Args:
        s: One label.

    Returns:
        float or None: The number, or None when the label is not a number.
    """
    if "_" in s:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return None if np.isnan(v) else v


def level_order(labels, levels=None, what="group"):
    """The order of the levels of a grouping column, one rule in both engines.

    With --levels the order is the listed order: a value of the data that the
    list does not name stops the run, and a listed value that no row holds is
    left out. Without --levels the levels sort as numbers when every label is a
    number, else in code-point order. The R engine has the same rule.

    Args:
        labels: The labels. A value that is not text is formatted with
            level_labels() first, and a missing value is left out.
        levels: The order that --levels gives, or None.
        what: The word for one level in the error message.

    Returns:
        list: The levels in order.

    Raises:
        ValueError: If a value of the data is missing from levels.
    """
    u = []
    for v in labels:
        lab = v if isinstance(v, str) else _level_label(v)
        if isinstance(lab, str) and lab not in u:
            u.append(lab)
    if levels:
        listed = []
        for s in levels:
            s = str(s).strip()
            if s not in listed:
                listed.append(s)
        miss = [v for v in u if v not in listed]
        if miss:
            raise ValueError("--levels does not list the {}(s) {} of the data".format(
                what, ", ".join(f"'{m}'" for m in miss)))
        return [v for v in listed if v in u]
    nums = [_level_number(v) for v in u]
    if u and all(n is not None for n in nums):
        return [v for _, v in sorted(zip(nums, u))]
    return sorted(u)
