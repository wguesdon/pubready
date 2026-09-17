import numpy as np

from pubready.util import cap_first, format_p, p_stars, shapiro_normal


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
