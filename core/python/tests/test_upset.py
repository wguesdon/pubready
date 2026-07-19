from types import SimpleNamespace

import pandas as pd

from pubplot.recipes.upset import recipe_upset
from pubplot.spec import spec_from_args


def test_recipe_upset_summarizes_sets_and_intersections():
    df = pd.DataFrame({
        "gene": [f"g{i}" for i in range(6)],
        "A": [1, 1, 0, 1, 0, 0],
        "B": [1, 0, 1, 1, 0, 0],
        "C": [0, 0, 1, 1, 1, 0],
    })
    a = SimpleNamespace(recipe="upset", data="sets.csv", x=None, y=None)
    spec = spec_from_args(a)

    res = recipe_upset(df, spec)

    stats = res["stats"]
    sets = stats[stats["type"] == "set"]
    assert dict(zip(sets["sets"], sets["size"])) == {"A": 3, "B": 3, "C": 3}
    # g5 belongs to no set, so it is dropped: 5 elements remain.
    assert res["test_meta"]["n_elements"] == 5
    assert (stats["type"] == "intersection").sum() >= 1
