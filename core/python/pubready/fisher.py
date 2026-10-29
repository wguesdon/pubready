"""Fisher's exact test as stats::fisher.test computes it.

For a 2 x 2 table: the two-sided P value, the conditional maximum likelihood
estimate of the odds ratio and its exact 95% confidence interval. R finds the
estimate and the interval with uniroot(), whose default tolerance is about
1.2e-4, so this module carries a port of the zeroin routine of R to land on the
same numbers. For an r x c table: the two-sided P value from a full enumeration
of the tables with the observed margins.
"""
import functools
import math

import numpy as np
from scipy import stats

EPS = float(np.finfo(float).eps)
DBL_MAX = float(np.finfo(float).max)
REL_ERR = 1 + 1e-7
MAX_TABLES = 1_000_000


def zeroin(f, ax, bx, fa, fb, tol, maxit):
    """Find a root of f in [ax, bx] with the R_zeroin2 routine of R.

    A line-by-line port of src/library/stats/src/zeroin.c (Brent's method), so
    the iterates and the returned root are those of uniroot().

    Args:
        f: The function of one argument.
        ax: The left end of the interval.
        bx: The right end of the interval.
        fa: f(ax).
        fb: f(bx).
        tol: The tolerance on the root.
        maxit: The maximum count of iterations.

    Returns:
        float: The root.
    """
    a, b = ax, bx
    c, fc = a, fa
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    left = maxit + 1
    while left:
        left -= 1
        prev_step = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol_act = 2 * EPS * abs(b) + tol / 2
        new_step = (c - b) / 2
        if abs(new_step) <= tol_act or fb == 0.0:
            return b
        if abs(prev_step) >= tol_act and abs(fa) > abs(fb):
            cb = c - b
            if a == c:
                t1 = fb / fa
                p = cb * t1
                q = 1.0 - t1
            else:
                q = fa / fc
                t1 = fb / fc
                t2 = fb / fa
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0))
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0)
            if p > 0:
                q = -q
            else:
                p = -p
            if p < (0.75 * cb * q - abs(tol_act * q) / 2) and p < abs(prev_step * q / 2):
                new_step = p / q
        if abs(new_step) < tol_act:
            new_step = tol_act if new_step > 0 else -tol_act
        a, fa = b, fb
        b += new_step
        fb = f(b)
        if (fb > 0 and fc > 0) or (fb < 0 and fc < 0):
            c, fc = a, fa
    return b


def uniroot(f, lower, upper, tol=EPS ** 0.25, maxiter=1000):
    """The root of f in [lower, upper], as stats::uniroot() returns it.

    Args:
        f: The function of one argument.
        lower: The left end of the interval.
        upper: The right end of the interval.
        tol: The tolerance, by default the one of uniroot().
        maxiter: The maximum count of iterations.

    Returns:
        float: The root.
    """
    flo = min(max(f(lower), -DBL_MAX), DBL_MAX)
    fup = min(max(f(upper), -DBL_MAX), DBL_MAX)
    return zeroin(f, lower, upper, flo, fup, tol, maxiter)


def _context(table):
    """The conditional hypergeometric setting of a 2 x 2 table, as fisher.test sets it.

    Args:
        table: The 2 x 2 table of counts.

    Returns:
        dict: The cell x, the margins m, n and k, the support and the log density
        of the central hypergeometric distribution on it.
    """
    t = np.asarray(table, dtype=int)
    m, n, k, x = int(t[:, 0].sum()), int(t[:, 1].sum()), int(t[0, :].sum()), int(t[0, 0])
    lo, hi = max(0, k - n), min(k, m)
    support = np.arange(lo, hi + 1)
    logdc = stats.hypergeom.logpmf(support, m + n, m, k)
    return {"x": x, "m": m, "n": n, "k": k, "lo": lo, "hi": hi,
            "support": support, "logdc": logdc}


def dnhyper(ctx, ncp):
    """The noncentral hypergeometric density on the support.

    Args:
        ctx: The setting from _context().
        ncp: The odds ratio, finite and above zero.

    Returns:
        numpy.ndarray: The density.
    """
    d = ctx["logdc"] + math.log(ncp) * ctx["support"]
    d = np.exp(d - d.max())
    return d / d.sum()


def mnhyper(ctx, ncp):
    """The mean of the noncentral hypergeometric distribution.

    Args:
        ctx: The setting from _context().
        ncp: The odds ratio.

    Returns:
        float: The mean.
    """
    if ncp == 0:
        return float(ctx["lo"])
    if ncp == math.inf:
        return float(ctx["hi"])
    return float(np.sum(ctx["support"] * dnhyper(ctx, ncp)))


def pnhyper(ctx, q, ncp=1.0, upper_tail=False):
    """The tail probability of the noncentral hypergeometric distribution.

    Args:
        ctx: The setting from _context().
        q: The quantile.
        ncp: The odds ratio.
        upper_tail: True for P(X >= q), else P(X <= q).

    Returns:
        float: The probability.
    """
    m, n, k, x = ctx["m"], ctx["n"], ctx["k"], ctx["x"]
    if ncp == 1:
        if upper_tail:
            return float(stats.hypergeom.sf(x - 1, m + n, m, k))
        return float(stats.hypergeom.cdf(x, m + n, m, k))
    if ncp == 0:
        return float(q <= ctx["lo"]) if upper_tail else float(q >= ctx["lo"])
    if ncp == math.inf:
        return float(q <= ctx["hi"]) if upper_tail else float(q >= ctx["hi"])
    d = dnhyper(ctx, ncp)
    keep = ctx["support"] >= q if upper_tail else ctx["support"] <= q
    return float(d[keep].sum())


def _mean_gap(t, ctx, invert):
    """mnhyper(t) - x, or mnhyper(1 / t) - x when invert is True."""
    return mnhyper(ctx, 1 / t if invert else t) - ctx["x"]


def _tail_gap(t, ctx, alpha, upper_tail, invert):
    """pnhyper(x, t) - alpha, or with 1 / t when invert is True."""
    return pnhyper(ctx, ctx["x"], 1 / t if invert else t, upper_tail) - alpha


def _mle(ctx):
    """The conditional maximum likelihood estimate of the odds ratio."""
    x = ctx["x"]
    if x == ctx["lo"]:
        return 0.0
    if x == ctx["hi"]:
        return math.inf
    mu = mnhyper(ctx, 1.0)
    if mu > x:
        return uniroot(functools.partial(_mean_gap, ctx=ctx, invert=False), 0.0, 1.0)
    if mu < x:
        return 1 / uniroot(functools.partial(_mean_gap, ctx=ctx, invert=True), EPS, 1.0)
    return 1.0


def _ncp_upper(ctx, alpha):
    """The upper confidence limit of the odds ratio."""
    if ctx["x"] == ctx["hi"]:
        return math.inf
    p = pnhyper(ctx, ctx["x"], 1.0)
    if p < alpha:
        return uniroot(functools.partial(_tail_gap, ctx=ctx, alpha=alpha, upper_tail=False,
                                         invert=False), 0.0, 1.0)
    if p > alpha:
        return 1 / uniroot(functools.partial(_tail_gap, ctx=ctx, alpha=alpha, upper_tail=False,
                                             invert=True), EPS, 1.0)
    return 1.0


def _ncp_lower(ctx, alpha):
    """The lower confidence limit of the odds ratio."""
    if ctx["x"] == ctx["lo"]:
        return 0.0
    p = pnhyper(ctx, ctx["x"], 1.0, upper_tail=True)
    if p > alpha:
        return uniroot(functools.partial(_tail_gap, ctx=ctx, alpha=alpha, upper_tail=True,
                                         invert=False), 0.0, 1.0)
    if p < alpha:
        return 1 / uniroot(functools.partial(_tail_gap, ctx=ctx, alpha=alpha, upper_tail=True,
                                             invert=True), EPS, 1.0)
    return 1.0


def fisher_2x2(table, conf_level=0.95):
    """Fisher's exact test of a 2 x 2 table, two-sided, as stats::fisher.test.

    Args:
        table: The 2 x 2 table of counts, rows by columns.
        conf_level: The confidence level of the interval.

    Returns:
        dict: p_value, odds_ratio (the conditional MLE), ci_low and ci_high.
    """
    ctx = _context(table)
    d = dnhyper(ctx, 1.0)
    p = float(np.sum(d[d <= d[ctx["x"] - ctx["lo"]] * REL_ERR]))
    alpha = (1 - conf_level) / 2
    return {"p_value": p, "odds_ratio": _mle(ctx),
            "ci_low": _ncp_lower(ctx, alpha), "ci_high": _ncp_upper(ctx, alpha)}


def column_fills(total, caps):
    """Every way to put a column total into cells with upper bounds.

    Args:
        total: The column total.
        caps: The remaining row totals, one per cell.

    Yields:
        tuple: The cell counts of the column.
    """
    if len(caps) == 1:
        if total <= caps[0]:
            yield (total,)
        return
    rest = sum(caps[1:])
    for v in range(max(0, total - rest), min(total, caps[0]) + 1):
        for tail in column_fills(total - v, caps[1:]):
            yield (v,) + tail


def _sum_tables(cols, caps, acc, state):
    """Add up the probability of every completion of a partial table.

    Args:
        cols: The column totals that are not filled yet.
        caps: The remaining row totals.
        acc: The sum of log(count!) over the cells filled so far.
        state: The running sums: const, obs, p, count and limit.

    Raises:
        ValueError: If the table count passes the limit.
    """
    if len(cols) == 1:
        state["count"] += 1
        if state["count"] > state["limit"]:
            raise ValueError("Fisher's exact test of this table needs more than "
                             f"{state['limit']} tables; use --test chisq or the R engine")
        logp = state["const"] - acc - sum(math.lgamma(c + 1) for c in caps)
        if logp <= state["obs"] + math.log(REL_ERR):
            state["p"] += math.exp(logp)
        return
    for fill in column_fills(cols[0], caps):
        rest = tuple(c - f for c, f in zip(caps, fill))
        _sum_tables(cols[1:], rest, acc + sum(math.lgamma(f + 1) for f in fill), state)


def fisher_rxc_p(table, limit=MAX_TABLES):
    """The two-sided P value of Fisher's exact test of an r x c table.

    The sum of the probabilities of every table with the observed margins whose
    probability is at most that of the observed table, as the network algorithm
    of fisher.test computes it, here by full enumeration.

    Args:
        table: The table of counts.
        limit: The largest count of tables to enumerate.

    Returns:
        float: The P value.

    Raises:
        ValueError: If the table count passes the limit.
    """
    t = np.asarray(table, dtype=int)
    if t.shape[0] < t.shape[1]:
        t = t.T
    rows = tuple(int(v) for v in t.sum(axis=1))
    cols = tuple(int(v) for v in t.sum(axis=0))
    n = int(t.sum())
    const = (sum(math.lgamma(r + 1) for r in rows) + sum(math.lgamma(c + 1) for c in cols)
             - math.lgamma(n + 1))
    obs = const - sum(math.lgamma(int(v) + 1) for v in t.ravel())
    state = {"const": const, "obs": obs, "p": 0.0, "count": 0, "limit": limit}
    _sum_tables(cols, rows, 0.0, state)
    return min(1.0, max(0.0, state["p"]))
