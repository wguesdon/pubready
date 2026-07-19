from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubplot.recipes.pca import recipe_pca
from pubplot.spec import spec_from_args


def test_pca_separates_two_groups():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 1, (10, 5))
    b = rng.normal(6, 1, (10, 5))
    df = pd.DataFrame(np.vstack([a, b]), columns=[f"f{i}" for i in range(5)])
    df["group"] = ["A"] * 10 + ["B"] * 10
    args = SimpleNamespace(recipe="pca", data="pca.csv", x=None, y=None, group="group")
    spec = spec_from_args(args)

    res = recipe_pca(df, spec)

    tm = res["test_meta"]
    assert tm["n_samples"] == 20
    assert tm["n_features"] == 5
    assert tm["permanova"]["p_value"] < 0.05
    assert res["stats"]["variance_explained_pct"].iloc[0] > 0
