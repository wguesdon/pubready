from types import SimpleNamespace

import pandas as pd

from pubplot.recipes.proportions import recipe_proportions
from pubplot.spec import spec_from_args


def test_proportions_chisq_2x2():
    df = pd.DataFrame({
        "arm": ["control"] * 20 + ["treated"] * 20,
        "response": (["no"] * 15 + ["yes"] * 5) + (["no"] * 5 + ["yes"] * 15),
    })
    a = SimpleNamespace(recipe="proportions", data="p.csv", x="arm", y="response", test="auto")
    spec = spec_from_args(a)

    res = recipe_proportions(df, spec)

    assert res["test_meta"]["name"] == "Pearson's chi-squared test"
    assert res["test_meta"]["p_value"] < 0.01
    assert res["test_meta"]["odds_ratio"] is not None
    assert len(res["stats"]) == 4  # 2 groups x 2 outcomes
