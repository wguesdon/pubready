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
