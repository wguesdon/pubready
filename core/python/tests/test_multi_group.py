"""multi_group_compare: the bracket rule and the compact letter display. The R
suite (test-multi-group.R) uses the same inputs and expects the same letters, so
the two engines label the groups identically."""
import itertools
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.spec import spec_from_args
from pubready.util import cld_letters


def cld_pairs():
    """The ten pairs of five groups, six of them significant.

    Returns:
        tuple: The first groups, the second groups and the P values.
    """
    pairs = list(itertools.combinations("ABCDE", 2))
    sig = {"A-C", "A-D", "A-E", "B-D", "B-E", "C-E"}
    p = [0.001 if f"{a}-{b}" in sig else 0.4 for a, b in pairs]
    return [a for a, _ in pairs], [b for _, b in pairs], p


def test_cld_letters_follows_the_insert_and_absorb_algorithm():
    g1, g2, p = cld_pairs()
    assert list(cld_letters(list("ABCDE"), g1, g2, p).values()) == ["a", "ab", "bc", "cd", "d"]
    assert list(cld_letters(list("ABCDE"), g1, g2, [0.5] * 10).values()) == ["a"] * 5
    assert list(cld_letters(["x", "y", "z"], ["x", "x", "y"], ["y", "z", "z"],
                            [0.01, 0.01, 0.01]).values()) == ["a", "b", "c"]


def six_groups():
    """Six groups in three pairs of equal means, with the same fixed noise as in R.

    Returns:
        pandas.DataFrame: The tidy table.
    """
    g = np.repeat([f"G{i}" for i in range(1, 7)], 6)
    mu = np.repeat([0, 0, 5, 5, 10, 10], 6)
    return pd.DataFrame({"group": g, "value": mu + 0.3 * np.sin(np.arange(1, 37) * 1.3)})


def mg_spec():
    """Build a multi_group_compare spec with the one-way ANOVA forced.

    Returns:
        dict: The spec.
    """
    return spec_from_args(SimpleNamespace(recipe="multi_group_compare", data="d.csv", x="group",
                                          y="value", test="anova"))


def test_multi_group_compare_switches_to_letters_above_six_significant_pairs():
    res = recipe_multi_group_compare(six_groups(), mg_spec())
    assert res["test_meta"]["n_significant"] == 12
    assert res["test_meta"]["display"] == "compact letter display"
    assert res["test_meta"]["letters"] == {"G1": "a", "G2": "a", "G3": "b", "G4": "b",
                                           "G5": "c", "G6": "c"}
    assert "compact letter display" in res["methods"]
    assert res["width"] == pytest.approx(7.6)


def test_multi_group_compare_keeps_brackets_for_a_few_pairs():
    res = recipe_multi_group_compare(six_groups().iloc[:18], mg_spec())
    assert res["test_meta"]["display"] == "brackets"
    texts = [t.get_text() for t in res["fig"].axes[0].texts]
    assert texts.count("****") == 2
