from types import SimpleNamespace

import pandas as pd

from pubready.recipes.proportions import recipe_proportions
from pubready.spec import spec_from_args


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


def test_the_first_outcome_sits_on_top_of_the_stack_as_in_the_legend():
    df = pd.DataFrame({
        "arm": ["control"] * 20 + ["treated"] * 20,
        "response": (["no"] * 15 + ["yes"] * 5) + (["no"] * 5 + ["yes"] * 15),
    })
    a = SimpleNamespace(recipe="proportions", data="p.csv", x="arm", y="response", test="auto")
    res = recipe_proportions(df, spec_from_args(a))
    ax = res["fig"].axes[0]
    leg = [t.get_text() for t in ax.get_legend().get_texts()]
    assert leg == ["no", "yes"]
    # The bar of "no" in the first group starts at 0.25 and ends at the top.
    bars = {c.get_label(): c for c in ax.containers}
    no_bar = bars["no"].patches[0]
    assert no_bar.get_y() == 0.25
    assert no_bar.get_y() + no_bar.get_height() == 1.0
