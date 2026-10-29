"""forest_meta: the DerSimonian-Laird truncation and the log-axis ticks. The R
suite (test-forest-meta.R) checks the same numbers."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.forest_meta import recipe_forest_meta
from pubready.spec import spec_from_args
from pubready.theme import log_ticks_ref


def fm_spec(**kw):
    """Build a forest_meta spec as the plot entry does.

    Args:
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": "forest_meta", "data": "d.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def homogeneous_studies():
    """Five studies that agree closely, so Cochran's Q is below its df.

    Returns:
        pandas.DataFrame: One row per study.
    """
    z = stats.norm.ppf(0.975)
    lr = np.log([0.80, 0.78, 0.82, 0.79, 0.81])
    se = np.array([0.20, 0.25, 0.30, 0.22, 0.28])
    return pd.DataFrame({"study": [f"s {i}" for i in range(1, 6)], "hr": np.exp(lr),
                         "lo": np.exp(lr - z * se), "hi": np.exp(lr + z * se)})


def test_forest_meta_truncates_tau2_and_i2_at_zero_when_q_is_below_its_df():
    res = recipe_forest_meta(homogeneous_studies(),
                             fm_spec(y="hr", ci_low="lo", ci_high="hi", label="study"))
    het = res["test_meta"]["heterogeneity"]
    assert het["q"] < het["df"]
    assert het["tau2"] == 0
    assert het["i2"] == 0
    assert res["test_meta"]["random"]["estimate"] == pytest.approx(
        res["test_meta"]["fixed"]["estimate"], rel=1e-10)
    st = res["stats"][res["stats"]["row_type"] == "study"]
    assert (st["weight_random_pct"] > 0).all()
    assert np.allclose(st["weight_random_pct"], st["weight_fixed_pct"], rtol=1e-10)


def test_log_ticks_ref_puts_a_tick_on_each_side_of_the_null_line():
    assert log_ticks_ref(0.4, 1.57, 1) == pytest.approx([0.5, 1, 2])
    assert log_ticks_ref(0.7, 0.95, 1) == pytest.approx([0.5, 1])
    assert log_ticks_ref(0.04, 18, 1) == pytest.approx([0.1, 1, 10])
    assert log_ticks_ref(0.15, 6, 1) == pytest.approx([0.2, 0.5, 1, 2, 5])
