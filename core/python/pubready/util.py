"""Small shared helpers (mirrors the R util functions)."""
import re

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
    """Shapiro-Wilk normality, guarded for small and large n.

    A sample with zero range (every value the same) is not normal and has no
    P, as shapiro_safe() of the R engine gives it.

    Args:
        x: The values.

    Returns:
        tuple: The P value and whether the values look normal (None below 3
        values).
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return (np.nan, None)
    if x.max() == x.min():
        return (np.nan, False)
    if len(x) > 5000:
        return (np.nan, True)
    try:
        p = float(stats.shapiro(x).pvalue)
    except Exception:
        return (np.nan, None)
    return (p, p > 0.05)


def is_zero_variance(v):
    """Whether a sample has zero variance: two or more finite values, all equal.

    Args:
        v: The values.

    Returns:
        bool: True for a sample with zero variance.
    """
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    return bool(len(v) >= 2 and v.max() == v.min())


def check_spread(v, recipe, what):
    """Stop when a variable cannot carry a test.

    The run stops when every value is the same, or when the variance is not a
    finite positive number, as for values near 1e300 (the squares overflow) or
    near 1e-300 (the squares underflow).

    Args:
        v: The values.
        recipe: The recipe name.
        what: The variable, for example "the outcome 'value'".

    Raises:
        ValueError: If the variable cannot carry a test, in the words of the R
            engine.
    """
    if is_zero_variance(v):
        raise ValueError(f"{recipe} cannot run because {what} has zero variance (all values are "
                         "equal)")
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) >= 2:
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            s = float(np.var(v, ddof=1))
        if not np.isfinite(s) or s <= 0:
            raise ValueError(f"{recipe} cannot run because the variance of {what} is not a finite "
                             f"positive number (the values run from {v.min():.3g} to "
                             f"{v.max():.3g}); rescale the values, for example to other units")


def zero_variance_groups(groups):
    """The names of the groups that have zero variance.

    Args:
        groups: A dict of the values of each group.

    Returns:
        list: The names, in the order of the dict.
    """
    return [str(g) for g, v in groups.items() if is_zero_variance(v)]


def check_test_method(recipe, method, allowed):
    """Stop on a --test value that a recipe does not know.

    Args:
        recipe: The recipe name.
        method: The --test value.
        allowed: The values that the recipe knows.

    Returns:
        str: The value in lower case.

    Raises:
        ValueError: If the recipe does not know the value. The message lists the
            known values, in the words of the R engine.
    """
    m = (method if method is not None else "auto").lower()
    if m and m not in allowed:
        raise ValueError(f"{recipe} does not know --test '{m}'; use one of: {', '.join(allowed)}")
    return m


def quote_join(names):
    """The quoted names of a list, joined as text, as quote_join() of the R engine.

    Args:
        names: The names.

    Returns:
        str: 'a'; 'a' and 'b'; 'a', 'b' and 'c'.
    """
    q = [f"'{n}'" for n in names]
    if len(q) <= 1:
        return "".join(q)
    return ", ".join(q[:-1]) + " and " + q[-1]


def zero_variance_message(test, subject, hint):
    """The message that stops a test that divides by a variance of zero.

    Args:
        test: The test, for example "Welch's one-way ANOVA".
        subject: The samples with their verb, for example "the group 'a' has".
        hint: What to run instead.

    Returns:
        str: The message, in the words of the R engine.
    """
    return f"{test} cannot run because {subject} zero variance (all values are equal); {hint}"


def zero_variance_subject(groups):
    """The subject of zero_variance_message() for a set of groups.

    Args:
        groups: The names of the groups.

    Returns:
        str: "the group 'a' has" or "the groups 'a' and 'b' have".
    """
    if len(groups) == 1:
        return f"the group {quote_join(groups)} has"
    return f"the groups {quote_join(groups)} have"


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


def rank_test_min_p(sizes, test):
    """The smallest P that a rank test can reach at the given group sizes.

    The same test function runs on data with no tie and complete separation.
    Independent groups take the ranks 1 to N in group order ("mann_whitney" for
    two groups, "kruskal" for the Kruskal-Wallis test); paired data take n
    distinct positive differences ("signed_rank"). The R engine has the same
    rule (rank_test_min_p).

    Args:
        sizes: The size of each group, or the count of pairs.
        test: "mann_whitney", "signed_rank" or "kruskal".

    Returns:
        float: The P value, or NaN when a group is empty.
    """
    sizes = [int(s) for s in sizes]
    if not sizes or min(sizes) < 1:
        return np.nan
    if test == "signed_rank":
        return wilcoxon_signed_rank_r(np.arange(1, sizes[0] + 1, dtype=float))[1]
    if len(sizes) < 2:
        return np.nan
    groups = np.split(np.arange(1, sum(sizes) + 1, dtype=float), np.cumsum(sizes)[:-1])
    if test == "mann_whitney":
        return mann_whitney_r(groups[0], groups[1])[1]
    return float(stats.kruskal(*groups).pvalue)


def rank_test_cannot_reach(min_p):
    """Whether the rank test cannot reach P < 0.05.

    Then --test auto keeps the parametric test after a failed normality check.

    Args:
        min_p: The smallest attainable P of the rank test.

    Returns:
        bool: True when the attainable P is above 0.05.
    """
    return bool(min_p is not None and np.isfinite(min_p) and min_p > 0.05)


def rank_floor_text(min_p, kept):
    """The sentence of a methods text on the attainable P of the rank test.

    Args:
        min_p: The smallest attainable P of the rank test.
        kept: Whether the parametric test was kept for this reason.

    Returns:
        str: The sentence, or "" when the rank test can reach P < 0.05. The R
        engine writes the same words.
    """
    if not rank_test_cannot_reach(min_p):
        return ""
    txt = ("The rank test cannot reach P < 0.05 at these group sizes (smallest attainable "
           f"P = {sig_text(min_p, 2)})")
    return txt + (", so the parametric test was kept." if kept else ".")


def cohens_d_r(a, b, paired=False):
    """Signed Cohen's d of two groups, first minus second.

    pingouin reports the magnitude only, so the sign that says which group is
    higher is restored here. Two independent groups use the pooled standard
    deviation; paired data use dz, the mean difference over its standard
    deviation. The R engine uses the same formulas (cohens_d_signed).

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


P_FLOOR_STUDENTIZED = 1e-12
"""Below this P the studentized range of Tukey and Games-Howell loses its precision.

R gives 0 or about 2.4e-13, Python 0, where the exact P can be 1e-17, so the P
is written as "< 1e-12", as p_floor_studentized of the R engine.
"""


def p_floor_column(p, floor=P_FLOOR_STUDENTIZED):
    """A column of P values as text when one of them is below the floor.

    Args:
        p: The P values.
        floor: The floor.

    Returns:
        list: The P values unchanged when none is below the floor, else text:
        "< 1e-12" below it and 15 significant digits for the others, as
        p_floor_column() of the R engine.
    """
    vals = [float(v) for v in p]
    if not any(np.isfinite(v) and v < floor for v in vals):
        return vals
    return [None if not np.isfinite(v) else (f"< {floor:g}" if v < floor else f"{v:.15g}")
            for v in vals]


def p_bracket_label(p, floored=False):
    """The label of a bracket in the p.format style.

    Args:
        p: The P value.
        floored: Whether the P comes from the studentized range.

    Returns:
        str: "p = 0.032", or "p < 1e-12" below the floor when floored.
    """
    if floored and np.isfinite(p) and p < P_FLOOR_STUDENTIZED:
        return f"p < {P_FLOOR_STUDENTIZED:g}"
    return f"p = {p:.2g}"


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


SCRIPT_LABEL_FN = '''def label(v):
    """Write a group value as the figure does: a whole number with no decimal point."""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)
'''
"""The label rule of level_labels() as source code, for an emitted script."""


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


#: One number in a time label: digits with an optional decimal part.
NUMBER = re.compile(r"[0-9]+(?:\.[0-9]+)?")


def one_number(s):
    """The one number that a time label holds, as .time_one_number() in R.

    A minus sign counts when it starts the label or follows a space ("Day -7").

    Args:
        s: The label.

    Returns:
        float or None: The number, or None when the label holds none or more
        than one.
    """
    found = list(NUMBER.finditer(s))
    if len(found) != 1:
        return None
    v = float(found[0].group())
    i = found[0].start()
    if i >= 1 and s[i - 1] == "-" and (i == 1 or s[i - 2] == " "):
        v = -v
    return v


def time_level_order(labels, levels=None, what="visit", option="--visit_levels"):
    """The order of a time column, one rule for sankey and timecourse.

    The listed order of --visit_levels (sankey) or --levels (timecourse) comes
    first; a value that the list does not name stops the run. Else the values
    sort as numbers when every value is a number; else by the one number that
    each value holds, as in "week 8" and "M12", with a tie in code-point order;
    else they keep the order of their first appearance. The R engine has the
    same rule (time_level_order).

    Args:
        labels: The time labels, in row order.
        levels: The listed order, or None.
        what: The word for one time value, for the messages.
        option: The option that lists the order, for the rule text.

    Returns:
        dict: The levels in order and the rule in words.

    Raises:
        ValueError: If a value of the data is missing from levels.
    """
    u = [v for v in dict.fromkeys(labels) if isinstance(v, str)]
    if levels:
        listed = [str(s).strip() for s in levels]
        miss = [v for v in u if v not in listed]
        if miss:
            raise ValueError(f"{option} does not list the {what}(s) "
                             + ", ".join(f"'{m}'" for m in miss) + " of the data")
        return {"levels": level_order(u, levels, what),
                "rule": f"the order that {option} gives"}
    if u and all(_level_number(v) is not None for v in u):
        return {"levels": level_order(u), "rule": "numeric order"}
    nums = [one_number(v) for v in u]
    if u and all(n is not None for n in nums):
        return {"levels": [v for _, v in sorted(zip(nums, u))],
                "rule": f"the order of the number that each {what} holds"}
    return {"levels": u,
            "rule": f"the order of first appearance, because the {what}s do not each hold "
                    "one number"}
