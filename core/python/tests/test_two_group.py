from types import SimpleNamespace

import pandas as pd
import pytest

from pubready.recipes.two_group_compare import _resolve, recipe_two_group_compare
from pubready.spec import spec_from_args


def test_resolve_picks_test_from_assumptions():
    assert _resolve("auto", True, True, False)["label"] == "Student's two-sample t-test"
    assert _resolve("auto", True, False, False)["label"] == "Welch two-sample t-test"
    assert _resolve("auto", False, True, False)["label"] == "Mann-Whitney U test"
    assert _resolve("welch_t", False, True, False)["var_equal"] is False
    assert _resolve("t", True, True, True)["label"] == "paired t-test"
    assert _resolve("wilcoxon", False, True, True)["label"] == "Wilcoxon signed-rank test"
    with pytest.raises(ValueError, match="unknown test method"):
        _resolve("nonsense", True, True, False)


def test_recipe_two_group_runs_forced_welch_end_to_end():
    df = pd.DataFrame({
        "group": ["ctrl"] * 8 + ["treat"] * 8,
        "value": [10, 11, 12, 13, 14, 10.5, 11.5, 12.5,
                  20, 21, 22, 23, 24, 20.5, 21.5, 22.5],
    })
    a = SimpleNamespace(recipe="two_group_compare", data="toy.csv",
                        x="group", y="value", test="welch_t")
    spec = spec_from_args(a)

    res = recipe_two_group_compare(df, spec)

    row = res["stats"].iloc[0]
    assert row["test"] == "Welch two-sample t-test"
    assert row["n1"] == 8 and row["n2"] == 8
    assert row["p_value"] < 0.01
    assert {row["group1"], row["group2"]} == {"ctrl", "treat"}


def test_recipe_two_group_rejects_non_binary_group():
    df = pd.DataFrame({"group": ["a", "b", "c"] * 3, "value": range(9)})
    a = SimpleNamespace(recipe="two_group_compare", data="toy.csv",
                        x="group", y="value", test="auto")
    spec = spec_from_args(a)
    with pytest.raises(ValueError, match="exactly 2 groups"):
        recipe_two_group_compare(df, spec)


def test_the_points_sit_at_fixed_offsets_as_in_the_r_engine():
    """The points of a group sit at golden-ratio offsets by input row, so they
    keep their places in every run. The R suite (test-theme.R) expects the same
    offsets."""
    df = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3, "value": [1, 2, 3, 4, 5, 6.0]})
    spec = spec_from_args(SimpleNamespace(recipe="two_group_compare", data="toy.csv",
                                          x="group", y="value", test="welch_t"))
    res = recipe_two_group_compare(df, spec)
    groups = [c.get_offsets() for c in res["fig"].axes[0].collections if len(c.get_offsets()) == 3]
    assert list(groups[0][:, 0]) == pytest.approx([0.0283282, -0.0633437, 0.0849845], abs=1e-6)
    assert list(groups[1][:, 0]) == pytest.approx([1.0283282, 0.9366563, 1.0849845], abs=1e-6)


def test_two_group_methods_say_when_the_user_set_the_test():
    df = pd.DataFrame({"group": ["ctrl"] * 8 + ["treat"] * 8,
                       "value": [10, 11, 12, 13, 14, 10.5, 11.5, 12.5,
                                 20, 21, 22, 23, 24, 20.5, 21.5, 22.5]})
    base = dict(recipe="two_group_compare", data="toy.csv", x="group", y="value",
                paired=False, p_adjust="none")
    mw = recipe_two_group_compare(df, spec_from_args(SimpleNamespace(**base, test="wilcoxon")))
    assert "The test was set by the user." in mw["methods"]
    assert "normality assumption was not met" not in mw["methods"]
    auto = recipe_two_group_compare(df, spec_from_args(SimpleNamespace(**base, test="auto")))
    assert "set by the user" not in auto["methods"]
