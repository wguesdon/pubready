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


def test_p_stars_boundaries_match_the_r_engine():
    assert p_stars(1e-4) == "***"
    assert p_stars(0.05) == "ns"


def test_p_phrase_writes_p_below_0001_without_an_equals_sign():
    from pubready.util import p_phrase
    assert p_phrase(1e-5) == "P < 0.001"
    assert p_phrase(0.0091) == "P = 0.0091"
    assert p_phrase(0.22) == "P = 0.22"
    assert p_phrase(float("nan")) == "P = NA"


def test_sig_text_keeps_the_trailing_zeros_and_drops_a_trailing_point():
    from pubready.util import sig_text
    assert sig_text(3.7) == "3.70"
    assert sig_text(3) == "3.00"
    assert sig_text(213) == "213"
    assert sig_text(-22.04) == "-22.0"
    assert sig_text(0.0104, 2) == "0.010"
    assert sig_text(7.1e-5, 2) == "7.1e-05"
    assert sig_text(float("nan")) == "NA"


def test_integer_columns_drops_the_trailing_zero_of_whole_numbers():
    import pandas as pd
    from pubready.bundle import integer_columns
    frame = pd.DataFrame({"n": [16470.0, np.nan], "df": [2.0, 147.0], "p": [0.5, np.nan],
                          "label": ["a", "b"]})
    out = integer_columns(frame)
    assert out.to_csv(index=False).splitlines() == ["n,df,p,label", "16470,2,0.5,a", ",147,,b"]
    assert frame["n"].dtype.kind == "f"
