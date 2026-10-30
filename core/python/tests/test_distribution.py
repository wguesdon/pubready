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
    """Check that q7 matches the type 7 quantile."""
    v = np.array([3.0, 1.0, 4.0, 1.5, 9.0, 2.6])
    for p in (0.1, 0.25, 0.5, 0.75, 0.9):
        assert q7(v, p) == pytest.approx(np.quantile(v, p), rel=1e-14)


def test_distribution_uses_the_freedman_diaconis_width():
    """Check that distribution uses the freedman diaconis width."""
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
    """Check that bw nrd0 restates silverman and its fallbacks."""
    v = np.array([1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 10.0])
    iqr = q7(v, 0.75) - q7(v, 0.25)
    lo = min(np.std(v, ddof=1), iqr / 1.34)
    assert bw_nrd0(v) == pytest.approx(0.9 * lo * len(v) ** -0.2)
    assert bw_nrd0(np.array([4.0, 4.0])) == pytest.approx(0.9 * 4.0 * 2 ** -0.2)


def test_bins_count_every_value_once_and_fall_back_to_sturges():
    """Check that bins count every value once and fall back to sturges."""
    v = np.array([0, 1, 1, 2, 3, 3, 3, 4, 5, 8], dtype=float)
    b = fd_bins(v)
    assert bin_counts(v, b).sum() == len(v)
    flat = np.array([5.0] * 20 + [1.0, 9.0])
    s = fd_bins(flat)
    assert s["rule"] == "Sturges"
    assert s["width"] == pytest.approx(8 / np.ceil(np.log2(22) + 1))


def test_distribution_reports_each_group_and_skips_shapiro_below_three():
    """Check that distribution reports each group and skips shapiro below three."""
    rng = np.random.default_rng(4)
    df = pd.DataFrame({"g": ["a"] * 30 + ["b"] * 2, "v": list(rng.normal(size=30)) + [0.1, 0.4]})
    res = recipe_distribution(df, dist_spec(x="g", y="v"))
    st = res["stats"]
    assert list(st["group"]) == ["a", "b"]
    assert list(st["n"]) == [30, 2]
    assert np.isnan(st["shapiro_p"].iloc[1])
    assert st["bin_width"].nunique() == 1
    assert list(res["qc"]["panels"]) == ["a", "b"]


def test_numeric_group_labels_match_the_r_engine():
    """Check that numeric group labels match the r engine."""
    df = pd.DataFrame({"dose": [0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, np.nan],
                       "v": [1.0, 2.0, 3.0, 2.0, 3.0, 5.0, 4.0, 6.0, 7.0, 1.0]})
    res = recipe_distribution(df, dist_spec(x="dose", y="v"))
    assert list(res["stats"]["group"]) == ["0.5", "1", "2"]


def test_values_near_1e300_or_1e_minus_300_stop_the_run():
    """A variance that is not finite at double precision stops the run (F83)."""
    import pytest
    from pubready.recipes.correlation_heatmap import recipe_correlation_heatmap
    big = [-7.858e299, 4.285e299, 9.2e299, 1.1e299, -3.3e299]
    with pytest.raises(ValueError, match="the variance of column 'value' is not finite"):
        recipe_distribution(pd.DataFrame({"value": big}), dist_spec(y="value"))
    tiny = [1.2e-300, 3.4e-300, 2.2e-300, 5.1e-300, 0.7e-300]
    with pytest.raises(ValueError, match="not finite at double precision"):
        recipe_distribution(pd.DataFrame({"value": tiny}), dist_spec(y="value"))
    df = pd.DataFrame({"big": big, "b": [1, 3, 2, 5, 4], "c": [2, 1, 4, 3, 6]})
    spec = spec_from_args(SimpleNamespace(recipe="correlation_heatmap", data="v.csv", x=None, y=None))
    with pytest.raises(ValueError, match="correlation_heatmap: the variance of column 'big'"):
        recipe_correlation_heatmap(df, spec)


def test_a_group_with_no_spread_gets_no_density_curve():
    """A group with no spread gets no density and no bandwidth, as in R (F125)."""
    from pubready.recipes.distribution import _curves
    df = pd.DataFrame({"condition": ["ctrl"] * 6 + ["pair"] * 2,
                       "value": [11.4, 13.1, 12.2, 10.9, 12.8, 11.7, 12.0, 12.0]})
    res = recipe_distribution(df, dist_spec(y="value", x="condition"))
    st = res["stats"].set_index("group")
    assert np.isnan(st.loc["pair", "bandwidth"])
    assert np.isfinite(st.loc["ctrl", "bandwidth"])
    assert "The group(s) pair have no spread" in res["methods"]
    assert _curves(np.array([12.0, 12.0]), fd_bins(df["value"].to_numpy())) is None
