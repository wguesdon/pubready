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


def test_proportions_writes_the_omnibus_test_and_names_the_odds_ratio():
    df = pd.DataFrame({"arm": ["control"] * 20 + ["treated"] * 20,
                       "response": ["no"] * 15 + ["yes"] * 5 + ["no"] * 5 + ["yes"] * 15})
    res = recipe_proportions(df, spec_from_args(SimpleNamespace(
        recipe="proportions", data="p.csv", x="arm", y="response", test="chisq")))
    assert res["stats"]["omnibus_statistic"].unique()[0] == res["test_meta"]["statistic"]
    assert res["stats"]["omnibus_p"].unique()[0] == res["test_meta"]["p_value"]
    assert "(chi-square(1) = " in res["methods"]
    assert "The odds ratio of no for control versus treated was 9.00." in res["methods"]


def test_multi_group_compare_writes_the_omnibus_test_on_each_row():
    from pubready.recipes.multi_group_compare import recipe_multi_group_compare
    df = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6 + ["c"] * 6,
                       "value": list(range(1, 7)) + [v + 2 for v in range(1, 7)]
                       + [v + 5 for v in range(1, 7)]})
    res = recipe_multi_group_compare(df, spec_from_args(SimpleNamespace(
        recipe="multi_group_compare", data="d.csv", x="group", y="value", test="anova")))
    assert list(res["stats"]["omnibus_p"].unique()) == [res["test_meta"]["omnibus"]["p_value"]]
    assert res["stats"]["omnibus_df2"].unique()[0] == 15


def test_proportions_draws_the_groups_in_the_order_of_levels():
    import pytest
    df = pd.DataFrame({"arm": ["b"] * 10 + ["c"] * 10 + ["a"] * 10,
                       "response": ["no", "yes"] * 15})
    base = dict(recipe="proportions", data="p.csv", x="arm", y="response", test="chisq")
    res = recipe_proportions(df, spec_from_args(SimpleNamespace(**base, levels="a,b,c")))
    assert list(dict.fromkeys(res["stats"]["group"])) == ["a", "b", "c"]
    plain = recipe_proportions(df, spec_from_args(SimpleNamespace(**base, levels=None)))
    assert list(dict.fromkeys(plain["stats"]["group"])) == ["b", "c", "a"]
    with pytest.raises(ValueError, match="'c'"):
        recipe_proportions(df, spec_from_args(SimpleNamespace(**base, levels="a,b")))
