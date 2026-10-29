"""The Fisher module: the port of zeroin, the 2 x 2 test of stats::fisher.test and
the r x c enumeration. The benchmark suite checks the numbers against R on public
data (mtcars_fisher, infert_fisher)."""
import math

import numpy as np
import pytest
from scipy import stats

from pubready.fisher import (_context, fisher_2x2, fisher_rxc_p, mnhyper, pnhyper,
                             uniroot)


def _square_gap(t):
    """t^2 - 2, whose positive root is the square root of 2."""
    return t * t - 2


def test_uniroot_lands_within_the_r_tolerance():
    root = uniroot(_square_gap, 0.0, 2.0)
    assert abs(root - math.sqrt(2)) < 1.3e-4


def test_the_2x2_p_value_matches_the_enumeration_and_scipy():
    table = np.array([[3, 1], [1, 3]])
    p = fisher_2x2(table)["p_value"]
    assert p == pytest.approx(fisher_rxc_p(table), rel=1e-12)
    assert p == pytest.approx(stats.fisher_exact(table)[1], rel=1e-9)


def test_the_estimate_solves_the_conditional_mean_and_the_ci_holds_alpha():
    table = np.array([[8, 4], [6, 14]])
    ft = fisher_2x2(table)
    ctx = _context(table)
    assert mnhyper(ctx, ft["odds_ratio"]) == pytest.approx(ctx["x"], abs=1e-3)
    assert pnhyper(ctx, ctx["x"], ft["ci_high"]) == pytest.approx(0.025, abs=1e-4)
    assert pnhyper(ctx, ctx["x"], ft["ci_low"], upper_tail=True) == pytest.approx(0.025, abs=1e-4)
    assert ft["ci_low"] < ft["odds_ratio"] < ft["ci_high"]


def test_a_zero_cell_gives_an_open_interval():
    ft = fisher_2x2(np.array([[0, 5], [6, 4]]))
    assert ft["odds_ratio"] == 0.0 and ft["ci_low"] == 0.0
    assert np.isfinite(ft["ci_high"])


def test_the_rxc_p_value_sums_to_one_for_a_proportional_table():
    assert fisher_rxc_p(np.array([[4, 8], [5, 10], [6, 12]])) == pytest.approx(1.0)


def test_the_rxc_enumeration_stops_above_its_limit():
    with pytest.raises(ValueError, match="use --test chisq or the R engine"):
        fisher_rxc_p(np.array([[10, 12, 9], [11, 8, 13], [9, 14, 10]]), limit=100)
