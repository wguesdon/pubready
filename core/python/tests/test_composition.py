"""Behaviour of the composition recipe on small synthetic inputs. The numbers
themselves are checked against public benchmarks in test_benchmarks.py."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.composition import bh_adjust, natural_key, recipe_composition
from pubready.spec import spec_from_args
from pubready.util import mann_whitney_r


def comp_spec(**kw):
    """Build a composition spec as the plot entry does.

    Args:
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": "composition", "data": "d.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def comp_data():
    """Ten samples in two conditions, three cell types, more B cells when treated.

    Returns:
        pandas.DataFrame: The long counts.
    """
    rng = np.random.default_rng(7)
    rows = []
    for i in range(1, 11):
        cond = "control" if i <= 5 else "treated"
        for cell in ("T", "B", "NK"):
            lam = 80 if (cell == "B" and cond == "treated") else 40
            rows.append({"sample": f"S{i}", "condition": cond, "cell": cell, "n": rng.poisson(lam)})
    return pd.DataFrame(rows)


def test_composition_compares_two_groups_with_the_r_wilcox_rule():
    """Check that composition compares two groups with the r wilcox rule."""
    res = recipe_composition(comp_data(), comp_spec(id="sample", x="cell", y="n", group="condition"))
    long = res["df_used"]
    assert long.groupby("sample")["proportion"].sum().to_numpy() == pytest.approx(np.ones(10))
    b = long[long["category"] == "B"]
    u, p = mann_whitney_r(b.loc[b["group"] == "control", "proportion"],
                          b.loc[b["group"] == "treated", "proportion"])
    row = res["stats"][res["stats"]["category"] == "B"].iloc[0]
    assert row["p_value"] == pytest.approx(p, rel=1e-12)
    assert row["statistic"] == u
    assert set(res["stats"]["test"]) == {"Mann-Whitney U test"}
    assert "compositional" in res["methods"]


def test_bh_adjust_matches_r_and_skips_missing_values():
    """Check that bh adjust matches r and skips missing values."""
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.2])
    # p.adjust(c(0.01, 0.04, NA, 0.03, 0.2), "BH") in R
    expected = [0.04, 0.0533333333333333, np.nan, 0.0533333333333333, 0.2]
    assert bh_adjust(p) == pytest.approx(expected, nan_ok=True)


def test_samples_sort_by_group_then_natural_name():
    """Check that samples sort by group then natural name."""
    df = comp_data()
    df.loc[df["sample"] == "S1", "sample"] = "S10b"
    res = recipe_composition(df, comp_spec(id="sample", x="cell", y="n", group="condition"))
    ids = list(pd.unique(res["df_used"]["sample"]))
    assert ids[:4] == ["S2", "S3", "S4", "S5"]
    assert sorted(["S10", "S2", "S1"], key=natural_key) == ["S1", "S2", "S10"]


def test_rows_are_counted_without_y_and_three_groups_use_kruskal():
    """Check that rows are counted without y and three groups use kruskal."""
    rng = np.random.default_rng(8)
    df = pd.DataFrame({"sample": np.repeat([f"s{i:02d}" for i in range(1, 10)], 30),
                       "cell": rng.choice(["x", "y", "z"], 270),
                       "grp": np.repeat(["a", "b", "c"], 90)})
    res = recipe_composition(df, comp_spec(id="sample", x="cell", group="grp"))
    assert set(res["stats"]["test"]) == {"Kruskal-Wallis test"}
    assert res["df_used"]["count"].sum() == 270
    assert any("counts as one observation" in s for s in res["clean_steps"])


def test_without_group_is_descriptive_and_a_sample_in_two_groups_is_refused():
    """Check that without group is descriptive and a sample in two groups is refused."""
    df = comp_data()
    res = recipe_composition(df, comp_spec(id="sample", x="cell", y="n"))
    assert res["stats"]["p_value"].isna().all()
    assert res["test_meta"]["name"] == "descriptive (no test)"
    df.loc[0, "condition"] = "treated"
    with pytest.raises(ValueError, match="one group"):
        recipe_composition(df, comp_spec(id="sample", x="cell", y="n", group="condition"))


def test_each_composition_strip_names_the_group_column():
    res = recipe_composition(comp_data(), comp_spec(id="sample", x="cell", y="n", group="condition"))
    titles = [ax.get_title() for ax in res["fig"].axes if ax.get_title()]
    assert titles and all(t.startswith("condition: ") for t in titles)


def test_composition_drops_a_category_with_a_total_count_of_zero():
    """A category with no count is dropped and logged, as in R (F88)."""
    from pubready.recipes.composition import recipe_composition
    df = comp_data()
    extra = pd.DataFrame({"sample": [f"S{i}" for i in range(1, 11)],
                          "condition": ["control"] * 5 + ["treated"] * 5,
                          "cell": "Platelet", "n": 0})
    res = recipe_composition(pd.concat([df, extra], ignore_index=True),
                             comp_spec(id="sample", x="cell", y="n", group="condition"))
    assert "Platelet" not in set(res["stats"]["category"])
    assert res["test_meta"]["n_categories"] == 3
    assert "Dropped 1 category(ies) with a total count of zero: 'Platelet'." in res["clean_steps"]
    assert "adjusted across the 3 categories" in res["methods"]
