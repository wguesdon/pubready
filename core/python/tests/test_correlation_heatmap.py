from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubplot.recipes.correlation_heatmap import recipe_correlation_heatmap
from pubplot.spec import spec_from_args


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
