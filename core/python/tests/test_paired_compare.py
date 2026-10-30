from types import SimpleNamespace

import pandas as pd

from pubready.recipes.paired_compare import recipe_paired_compare
from pubready.spec import spec_from_args


def test_paired_compare_paired_t_end_to_end():
    pre = [10, 11, 12, 13, 14, 15, 16, 17]
    post = [pre[i] - 5 + 0.3 * (-1) ** i for i in range(len(pre))]
    df = pd.DataFrame({
        "subject": [f"s{i}" for i in range(8)] * 2,
        "condition": ["pre"] * 8 + ["post"] * 8,
        "value": pre + post,
    })
    a = SimpleNamespace(recipe="paired_compare", data="p.csv", x="condition", y="value",
                        id="subject", test="t", levels="pre,post")
    spec = spec_from_args(a)

    res = recipe_paired_compare(df, spec)

    row = res["stats"].iloc[0]
    assert row["n_pairs"] == 8
    assert row["test"] == "paired t-test"
    assert row["condition1"] == "pre" and row["condition2"] == "post"
    assert row["p_value"] < 0.001


def test_paired_compare_needs_id():
    df = pd.DataFrame({"condition": ["pre", "post"], "value": [1.0, 2.0]})
    a = SimpleNamespace(recipe="paired_compare", data="p.csv", x="condition", y="value",
                        id=None, test="auto")
    spec = spec_from_args(a)
    try:
        recipe_paired_compare(df, spec)
        assert False, "expected ValueError"
    except ValueError as e:
        assert "--id" in str(e)


def test_paired_compare_labels_p_below_1e4_with_four_stars():
    pre = [10, 11, 12, 13, 14, 15, 16, 17]
    post = [pre[i] - 5 + 0.3 * (-1) ** (i + 1) for i in range(len(pre))]
    df = pd.DataFrame({"subject": [f"s{i}" for i in range(8)] * 2,
                       "condition": ["pre"] * 8 + ["post"] * 8, "value": pre + post})
    a = SimpleNamespace(recipe="paired_compare", data="p.csv", x="condition", y="value",
                        id="subject", test="t")
    res = recipe_paired_compare(df, spec_from_args(a))
    assert res["stats"]["p_value"].iloc[0] < 1e-4
    assert "****" in [t.get_text() for t in res["fig"].axes[0].texts]


def test_paired_compare_methods_name_the_source_of_the_test_choice():
    pre = [10, 11, 12, 13, 14, 15, 16, 17]
    post = [pre[i] - 5 + 0.3 * (-1) ** i for i in range(len(pre))]
    df = pd.DataFrame({"subject": [f"s{i}" for i in range(8)] * 2,
                       "condition": ["pre"] * 8 + ["post"] * 8,
                       "value": pre + post})
    base = dict(recipe="paired_compare", data="p.csv", x="condition", y="value", id="subject")
    forced = recipe_paired_compare(df, spec_from_args(SimpleNamespace(**base, test="t")))
    assert "The test was set by the user." in forced["methods"]
    assert "chosen from the Shapiro-Wilk" not in forced["methods"]
    auto = recipe_paired_compare(df, spec_from_args(SimpleNamespace(**base, test="auto")))
    assert "chosen from the Shapiro-Wilk normality" in auto["methods"]
