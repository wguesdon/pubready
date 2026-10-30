from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubready.recipes.correlation_heatmap import recipe_correlation_heatmap
from pubready.spec import spec_from_args


def test_recipe_correlation_heatmap_end_to_end():
    rng = np.random.default_rng(1)
    f = rng.normal(0, 1, 60)
    df = pd.DataFrame({
        "a": f + rng.normal(0, 0.2, 60),
        "b": f + rng.normal(0, 0.2, 60),
        "c": -f + rng.normal(0, 0.2, 60),
        "d": rng.normal(0, 1, 60),
    })
    a = SimpleNamespace(recipe="correlation_heatmap", data="vars.csv",
                        x=None, y=None, test="pearson")
    spec = spec_from_args(a)

    res = recipe_correlation_heatmap(df, spec)

    # 4 variables -> 6 unique pairs.
    assert len(res["stats"]) == 6
    assert set(res["stats"].columns) >= {"var1", "var2", "r", "p_value", "p_adj", "significant"}
    assert res["stats"]["r"].max() > 0.8    # a & b co-vary
    assert res["stats"]["r"].min() < -0.8   # c is anti-correlated
    assert res["test_meta"]["n_variables"] == 4


def ch_area_d():
    """Four variables on 12 rows with fixed values, as in the R suite.

    Returns:
        pandas.DataFrame: The table.
    """
    return pd.DataFrame({
        "sample": [f"S{i:02d}" for i in range(1, 13)],
        "a": [1.2, 2.3, 3.1, 4.8, 5.2, 6.9, 7.1, 8.4, 9.0, 10.2, 11.5, 12.1],
        "b": [2.0, 2.9, 3.8, 5.1, 5.0, 7.2, 6.8, 8.9, 9.4, 9.9, 12.0, 12.6],
        "c": [5.5, 3.2, 6.1, 2.2, 4.4, 1.9, 6.6, 3.3, 5.0, 2.7, 4.1, 3.9],
        "d": [9.1, 8.2, 7.7, 6.5, 6.6, 5.1, 4.9, 3.8, 3.1, 2.2, 1.7, 0.9],
    })


def ch_spec(**kw):
    """A correlation_heatmap spec.

    Args:
        **kw: Extra command-line options.

    Returns:
        dict: The spec.
    """
    return spec_from_args(SimpleNamespace(recipe="correlation_heatmap", data="v.csv",
                                          x=None, y=None, **kw))


def test_correlation_heatmap_keeps_a_mostly_numeric_column():
    """One text cell does not drop the column; the log names it, as in R (F42)."""
    df = ch_area_d()
    df["c"] = df["c"].astype(str)
    df.loc[3, "c"] = "<0.1"
    res = recipe_correlation_heatmap(df, ch_spec())
    assert res["test_meta"]["n_variables"] == 4
    steps = " | ".join(res["clean_steps"])
    assert "Coerced column 'c' to numeric; 1 value(s) became NA, for example '<0.1'." in steps
    assert "Left out the text column(s) 'sample'." in steps


def test_correlation_heatmap_stops_on_a_pair_with_too_few_shared_rows():
    """A pair with fewer than 3 shared rows stops the run, as in R (F79, F87)."""
    import pytest
    df = ch_area_d()
    nan = float("nan")
    df["e"] = [nan] * 6 + [1.0, 2.5, 1.7] + [nan] * 3
    df["f"] = [0.3, 0.9, 1.4, 2.2, 2.9, 3.1, nan, nan, 4.4, 5.0, 5.8, 6.1]
    with pytest.raises(ValueError, match=r"1 pair\(s\) of columns have no correlation: 'e' and "
                                         r"'f' \(1 shared complete row\(s\)\)"):
        recipe_correlation_heatmap(df, ch_spec())
