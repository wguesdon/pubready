"""The v3 fixes of the comparison recipes.

core/r/tests/testthat/test-v3-compare.R holds the same cases for the R engine.
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.correlation import recipe_correlation
from pubready.recipes.feature_panel import recipe_feature_panel
from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.recipes.paired_compare import recipe_paired_compare
from pubready.recipes.two_group_compare import recipe_two_group_compare
from pubready.spec import spec_from_args
from pubready.util import shapiro_normal


def v3_spec(recipe, **kw):
    """Build a figure specification for a test.

    Args:
        recipe: The recipe name.
        **kw: The command-line options.

    Returns:
        dict: The specification.
    """
    opts = {"x": None, "y": None}
    opts.update(kw)
    return spec_from_args(SimpleNamespace(recipe=recipe, data="toy.csv", **opts))


# F43: a constant sample has one status, and a test that divides by a variance
# of zero stops with a message that names the samples.
def test_a_constant_sample_is_not_normal_and_has_no_shapiro_p():
    p, normal = shapiro_normal([5.0] * 6)
    assert np.isnan(p)
    assert normal is False


def test_a_forced_t_test_on_two_constant_groups_stops_and_names_them():
    df = pd.DataFrame({"group": ["a"] * 4 + ["b"] * 4, "value": [5.0] * 4 + [6.0] * 4})
    with pytest.raises(ValueError, match="the groups 'a' and 'b' have zero variance"):
        recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="value",
                                             test="welch_t"))


def test_welch_anova_on_a_constant_group_stops_and_names_it():
    df = pd.DataFrame({"group": ["het"] * 4 + ["ko"] * 4 + ["wt"] * 4,
                       "value": [5, 5, 5, 5, 1, 2, 3, 4, 2, 4, 6, 9]})
    with pytest.raises(ValueError, match="the group 'het' has zero variance"):
        recipe_multi_group_compare(df, v3_spec("multi_group_compare", x="group", y="value",
                                               test="welch"))


def test_a_paired_t_test_on_constant_differences_stops():
    df = pd.DataFrame({"subject": list(range(1, 6)) * 2, "condition": ["pre"] * 5 + ["post"] * 5,
                       "value": list(range(1, 6)) + [v + 3 for v in range(1, 6)]})
    with pytest.raises(ValueError, match="within-subject differences have zero variance"):
        recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value",
                                          id="subject", test="t"))


def test_correlation_stops_on_a_constant_column_and_names_it():
    df = pd.DataFrame({"a": [1, 2, 3, 4, 5.0], "b": [2.0] * 5})
    with pytest.raises(ValueError, match="the column 'b' is constant"):
        recipe_correlation(df, v3_spec("correlation", x="a", y="b"))


def test_feature_panel_reports_a_constant_feature_as_not_tested():
    df = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6,
                       "f1": [1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5],
                       "floor": [0.5] * 12})
    res = recipe_feature_panel(df, v3_spec("feature_panel", x="group"))
    row = res["stats"].set_index("feature").loc["floor"]
    assert row["note"] == "not tested: constant"
    assert np.isnan(row["p_value"])
    assert res["test_meta"]["n_constant"] == 1
    assert res["test_meta"]["assumptions"]["normality"]["n_checked"] == 2
    assert "1 feature(s) were constant" in res["methods"]
