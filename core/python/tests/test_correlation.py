from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.correlation import _resolve, recipe_correlation
from pubready.spec import spec_from_args


def test_resolve_picks_method_from_normality():
    assert _resolve("auto", True, True)["method"] == "pearson"
    assert _resolve("auto", True, False)["method"] == "spearman"
    assert _resolve("auto", False, True)["method"] == "spearman"
    assert _resolve("spearman", True, True)["method"] == "spearman"
    with pytest.raises(ValueError, match="unknown correlation method"):
        _resolve("nonsense", True, True)


def test_recipe_correlation_runs_end_to_end():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 40)
    y = 2 * x + rng.normal(0, 1, 40)
    df = pd.DataFrame({"x": x, "y": y})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="y", test="pearson")
    spec = spec_from_args(a)

    res = recipe_correlation(df, spec)

    row = res["stats"].iloc[0]
    assert row["method"] == "Pearson correlation"
    assert row["n"] == 40
    assert row["estimate"] > 0.9
    assert row["p_value"] < 1e-6


def test_recipe_correlation_needs_present_columns():
    df = pd.DataFrame({"x": [1, 2, 3]})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="missing", test="auto")
    spec = spec_from_args(a)
    with pytest.raises(ValueError, match="not found"):
        recipe_correlation(df, spec)


def test_the_correlation_subtitle_names_the_method_the_coefficient_and_p():
    from pubready.recipes.correlation import cor_label
    assert (cor_label("Pearson correlation", "r", 0.7123, 2.1e-10)
            == "Pearson correlation, r = 0.71, P = 2.1e-10")
    assert (cor_label("Spearman rank correlation", "rho", -0.456, 0.0312)
            == "Spearman rank correlation, rho = -0.46, P = 0.031")
    x = np.arange(1, 21, dtype=float)
    df = pd.DataFrame({"x": x, "y": x + np.sin(x)})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="y", test="kendall")
    res = recipe_correlation(df, spec_from_args(a))
    assert res["fig"].axes[0].get_title().startswith("Kendall rank correlation, tau = ")
