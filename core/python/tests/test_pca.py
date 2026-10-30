from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubready.recipes.pca import recipe_pca
from pubready.spec import spec_from_args


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


def pca_id_data():
    """Twelve samples with a numeric id column and four noise features.

    Returns:
        pandas.DataFrame: The table.
    """
    rng = np.random.RandomState(2)
    df = pd.DataFrame(rng.normal(0, 1, (12, 4)), columns=["f1", "f2", "f3", "f4"])
    df.insert(0, "group", ["A"] * 4 + ["B"] * 4 + ["C"] * 4)
    df.insert(0, "sample_id", range(1, 13))
    return df


def test_pca_leaves_out_the_id_takes_features_and_names_them():
    """--id and --label are not features, --features sets them, the log names them (F47)."""
    import pytest
    df = pca_id_data()
    base = dict(recipe="pca", data="p.csv", x=None, y=None, group="group")
    plain = recipe_pca(df, spec_from_args(SimpleNamespace(**base)))
    assert plain["test_meta"]["n_features"] == 5
    assert "'sample_id', 'f1'" in plain["clean_steps"][0]
    res = recipe_pca(df, spec_from_args(SimpleNamespace(**base, id="sample_id")))
    assert res["test_meta"]["n_features"] == 4
    assert "4 numeric features: 'f1', 'f2', 'f3', 'f4'." in res["clean_steps"][0]
    assert any("Left out the --id and --label column(s) 'sample_id'" in s
               for s in res["clean_steps"])
    sel = recipe_pca(df, spec_from_args(SimpleNamespace(**base, features="f1,f3,f4")))
    assert sel["test_meta"]["n_features"] == 3
    with pytest.raises(ValueError, match="'f9'"):
        recipe_pca(df, spec_from_args(SimpleNamespace(**base, features="f1,f9")))
