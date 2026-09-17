import numpy as np
import pytest
from scipy import stats

from pubready.util import (cap_first, format_p, mann_whitney_r, p_stars,
                           rank_biserial_paired, rank_biserial_unpaired,
                           shapiro_normal, wilcoxon_signed_rank_r)


def test_cap_first():
    assert cap_first("hello") == "Hello"
    assert cap_first("Hello") == "Hello"
    assert cap_first("") == ""


def test_p_stars_thresholds():
    assert p_stars(1e-5) == "****"
    assert p_stars(5e-4) == "***"
    assert p_stars(5e-3) == "**"
    assert p_stars(2e-2) == "*"
    assert p_stars(0.2) == "ns"
    assert p_stars(np.nan) == "ns"


def test_format_p():
    assert format_p(0.0005) == "< 0.001"
    assert format_p(0.0123) == "0.012"
    assert format_p(np.nan) == "NA"


def test_shapiro_normal_guards():
    p, is_normal = shapiro_normal([1.0, 2.0])
    assert np.isnan(p) and is_normal is None

    p, is_normal = shapiro_normal(np.arange(6000.0))
    assert np.isnan(p) and is_normal is True

    p, is_normal = shapiro_normal(np.random.default_rng(0).normal(size=200))
    assert is_normal is True


def test_mann_whitney_r_takes_the_exact_route_without_ties():
    a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    b = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    u, p = mann_whitney_r(a, b)
    assert u == 0.0
    # Two-sided exact p for complete separation is 2 / C(12, 6).
    assert p == pytest.approx(2 / 924)


def test_mann_whitney_r_falls_back_to_the_approximation_with_ties():
    a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    b = [6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
    _, p = mann_whitney_r(a, b)
    expected = stats.mannwhitneyu(a, b, alternative="two-sided", method="asymptotic",
                                  use_continuity=True).pvalue
    assert p == pytest.approx(expected)


def test_wilcoxon_signed_rank_r_reports_v_and_drops_zeros():
    v, _ = wilcoxon_signed_rank_r([1.0, 2.0, 3.0, -4.0])
    assert v == 6.0  # ranks 1, 2 and 3 carry the positive differences

    v_zero, p_zero = wilcoxon_signed_rank_r([0.0, 1.0, 2.0, 3.0, -4.0])
    expected = stats.wilcoxon([1.0, 2.0, 3.0, -4.0], alternative="two-sided",
                              method="approx", correction=True).pvalue
    assert v_zero == 6.0
    assert p_zero == pytest.approx(expected)


def test_rank_biserial_covers_both_designs():
    assert rank_biserial_unpaired(0, 6, 6) == -1.0
    assert rank_biserial_unpaired(36, 6, 6) == 1.0
    assert rank_biserial_unpaired(18, 6, 6) == 0.0
    assert rank_biserial_paired([1.0, 2.0, 3.0, -4.0]) == pytest.approx(0.2)
    assert np.isnan(rank_biserial_paired([0.0, 0.0]))
