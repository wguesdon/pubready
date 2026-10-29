from types import SimpleNamespace

import pandas as pd

from pubready.recipes.upset import recipe_upset
from pubready.spec import spec_from_args


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


def test_recipe_upset_caps_the_drawn_intersections_at_top_n():
    # Ties sort by degree and then by set name, the order of the R engine
    # (test-upset.R has the same case).
    df = pd.DataFrame({
        "gene": [f"g{i}" for i in range(9)],
        "A": [1, 1, 1, 0, 0, 1, 0, 1, 0],
        "B": [0, 0, 0, 1, 1, 1, 0, 1, 0],
        "C": [0, 0, 0, 0, 0, 0, 1, 1, 1],
    })
    a = SimpleNamespace(recipe="upset", data="sets.csv", x=None, y=None, top_n=3)
    res = recipe_upset(df, spec_from_args(a))

    inter = res["stats"][res["stats"]["type"] == "intersection"]
    assert list(inter["sets"]) == ["A", "B", "C", "A & B", "A & B & C"]
    assert list(inter["size"]) == [3, 2, 2, 1, 1]
    assert res["test_meta"]["n_intersections"] == 5
    assert res["test_meta"]["n_intersections_shown"] == 3
    texts = [t.get_text() for t in res["fig"].texts]
    assert "The 3 largest of 5 intersections; 2 smaller ones are not shown" in texts
