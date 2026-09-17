"""The spider recipe, Python engine. Mirrors core/r/tests/testthat/test-spider_response.R."""
from types import SimpleNamespace

import pandas as pd
import pytest

from pubready.recipes.spider_response import recipe_spider_response
from pubready.spec import spec_from_args


def sample_spider():
    """Build three patients whose paths land in the three categories.

    Returns:
        The table, one row per measurement.
    """
    return pd.DataFrame({
        "patient": ["p1"] * 3 + ["p2"] * 3 + ["p3"] * 3,
        "arm": ["treated"] * 3 + ["control"] * 3 + ["treated"] * 3,
        "week": [0, 8, 16] * 3,
        "size": [100, 80, 50,      # p1 falls to -50%: partial response
                 100, 110, 130,    # p2 rises to +30%: progressive disease
                 100, 95, 105],    # p3 stays inside both lines: stable
    })


def spider_spec(**extra):
    """Build the specification a spider call arrives with.

    Args:
        **extra: The fields to override.

    Returns:
        The specification.
    """
    fields = {"recipe": "spider_response", "data": "spider.csv", "x": "week",
              "y": "size", "id": "patient"}
    fields.update(extra)
    return spec_from_args(SimpleNamespace(**fields))


def test_the_change_is_measured_from_the_first_visit():
    res = recipe_spider_response(sample_spider(), spider_spec())

    first = res["df_used"][res["df_used"]["week"] == 0]
    assert set(first["change"]) == {0.0}
    assert res["test_meta"]["n_patients"] == 3
    assert res["test_meta"]["n_measurements"] == 9


def test_the_two_thresholds_put_each_patient_in_a_category():
    res = recipe_spider_response(sample_spider(), spider_spec())

    got = dict(zip(res["stats"]["patient"], res["stats"]["category_from_change"]))
    assert got["p1"] == "partial response"
    assert got["p2"] == "progressive disease"
    assert got["p3"] == "stable disease"
    assert res["test_meta"]["partial_response"] == 1
    assert res["test_meta"]["progressive_disease"] == 1
    assert res["test_meta"]["stable_disease"] == 1


def test_a_threshold_the_reader_moves_changes_the_category():
    # p3 reaches +5% and -5%. A pair of thresholds inside that range makes it a
    # partial response, which is what says the parameters are read.
    res = recipe_spider_response(sample_spider(),
                                 spider_spec(pd_threshold=4, pr_threshold=-4))

    got = dict(zip(res["stats"]["patient"], res["stats"]["category_from_change"]))
    assert got["p3"] == "partial response"
    assert res["test_meta"]["pd_threshold"] == 4
    assert res["test_meta"]["pr_threshold"] == -4


def test_a_patient_with_a_zero_baseline_is_dropped_and_recorded():
    df = pd.concat([sample_spider(), pd.DataFrame({
        "patient": ["p4"] * 3, "arm": ["control"] * 3,
        "week": [0, 8, 16], "size": [0, 10, 20]})], ignore_index=True)

    res = recipe_spider_response(df, spider_spec())

    assert res["test_meta"]["n_patients"] == 3
    assert res["test_meta"]["n_dropped_patients"] == 1
    assert "p4" not in set(res["stats"]["patient"])
    assert any("non-zero baseline" in step for step in res["clean_steps"])


def test_the_recipe_names_the_column_it_cannot_find():
    with pytest.raises(ValueError, match="subject"):
        recipe_spider_response(sample_spider(), spider_spec(id="subject"))
    with pytest.raises(ValueError, match="--id"):
        recipe_spider_response(sample_spider(), spider_spec(id=None))


def test_the_group_column_reaches_the_table():
    res = recipe_spider_response(sample_spider(), spider_spec(group="arm"))

    assert "group" in res["stats"].columns
    assert sorted(set(res["stats"]["group"])) == ["control", "treated"]


def test_a_patient_who_only_grows_carries_its_smallest_growth():
    """The baseline reads zero for everyone, so it cannot be the best change."""
    df = pd.DataFrame({
        "patient": ["p1"] * 3,
        "week": [0, 8, 16],
        "size": [100, 105, 112],
    })

    res = recipe_spider_response(df, spider_spec())

    row = res["stats"].iloc[0]
    assert row["best_change_pct"] == pytest.approx(5.0)
    assert row["worst_change_pct"] == pytest.approx(12.0)
    assert row["category_from_change"] == "stable disease"


def test_a_patient_with_no_visit_after_the_baseline_is_dropped():
    """A response needs a second measurement, and RECIST calls that patient
    not evaluable."""
    df = pd.concat([sample_spider(), pd.DataFrame({
        "patient": ["p5"], "arm": ["treated"], "week": [0], "size": [100.0]})],
        ignore_index=True)

    res = recipe_spider_response(df, spider_spec())

    assert res["test_meta"]["n_without_follow_up"] == 1
    assert "p5" not in set(res["stats"]["patient"])
    assert any("after the baseline" in step for step in res["clean_steps"])
