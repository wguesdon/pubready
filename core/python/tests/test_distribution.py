"""Behaviour of the distribution recipe on small synthetic inputs. The numbers
themselves are checked against public benchmarks in test_benchmarks.py."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.distribution import bin_counts, bw_nrd0, fd_bins, q7, recipe_distribution
from pubready.spec import spec_from_args


def dist_spec(**kw):
    """Build a distribution spec as the plot entry does.

    Args:
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": "distribution", "data": "d.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def test_q7_matches_the_type_7_quantile():
    v = np.array([3.0, 1.0, 4.0, 1.5, 9.0, 2.6])
    for p in (0.1, 0.25, 0.5, 0.75, 0.9):
        assert q7(v, p) == pytest.approx(np.quantile(v, p), rel=1e-14)


def test_distribution_uses_the_freedman_diaconis_width():
    rng = np.random.default_rng(3)
    v = rng.normal(50, 10, 120)
    res = recipe_distribution(pd.DataFrame({"value": v}), dist_spec(y="value"))
    row = res["stats"].iloc[0]
    w = 2 * (q7(v, 0.75) - q7(v, 0.25)) * len(v) ** (-1 / 3)
    assert row["bin_width"] == pytest.approx(w, rel=1e-12)
    assert res["test_meta"]["bins"]["start"] == pytest.approx(np.floor(v.min() / w) * w)
    assert row["group"] == "all"
    assert row["shapiro_p"] == pytest.approx(stats.shapiro(v).pvalue, rel=1e-12)


def test_bw_nrd0_restates_silverman_and_its_fallbacks():
    v = np.array([1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 10.0])
    iqr = q7(v, 0.75) - q7(v, 0.25)
    lo = min(np.std(v, ddof=1), iqr / 1.34)
    assert bw_nrd0(v) == pytest.approx(0.9 * lo * len(v) ** -0.2)
    assert bw_nrd0(np.array([4.0, 4.0])) == pytest.approx(0.9 * 4.0 * 2 ** -0.2)


def test_bins_count_every_value_once_and_fall_back_to_sturges():
    v = np.array([0, 1, 1, 2, 3, 3, 3, 4, 5, 8], dtype=float)
    b = fd_bins(v)
    assert bin_counts(v, b).sum() == len(v)
    flat = np.array([5.0] * 20 + [1.0, 9.0])
    s = fd_bins(flat)
    assert s["rule"] == "Sturges"
    assert s["width"] == pytest.approx(8 / np.ceil(np.log2(22) + 1))


def test_distribution_reports_each_group_and_skips_shapiro_below_three():
    rng = np.random.default_rng(4)
    df = pd.DataFrame({"g": ["a"] * 30 + ["b"] * 2, "v": list(rng.normal(size=30)) + [0.1, 0.4]})
    res = recipe_distribution(df, dist_spec(x="g", y="v"))
    st = res["stats"]
    assert list(st["group"]) == ["a", "b"]
    assert list(st["n"]) == [30, 2]
    assert np.isnan(st["shapiro_p"].iloc[1])
    assert st["bin_width"].nunique() == 1
    assert list(res["qc"]["panels"]) == ["a", "b"]
