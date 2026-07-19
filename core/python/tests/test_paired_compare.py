from types import SimpleNamespace

import pandas as pd

from pubplot.recipes.paired_compare import recipe_paired_compare
from pubplot.spec import spec_from_args


def test_paired_compare_paired_t_end_to_end():
    pre = [10, 11, 12, 13, 14, 15, 16, 17]
    post = [pre[i] - 5 + 0.3 * (-1) ** i for i in range(len(pre))]
    df = pd.DataFrame({
        "subject": [f"s{i}" for i in range(8)] * 2,
        "condition": ["pre"] * 8 + ["post"] * 8,
        "value": pre + post,
    })
    a = SimpleNamespace(recipe="paired_compare", data="p.csv", x="condition", y="value",
                        id="subject", test="t")
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
