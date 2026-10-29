"""Behaviour of the six recipes added in v1.0.0 on small synthetic inputs. The
numbers themselves are checked against public benchmarks in test_benchmarks.py."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.bland_altman import recipe_bland_altman
from pubready.recipes.dose_response import fpl, recipe_dose_response
from pubready.recipes.forest_meta import recipe_forest_meta
from pubready.recipes.roc_curve import recipe_roc_curve
from pubready.recipes.swimmer import recipe_swimmer
from pubready.recipes.timecourse import recipe_timecourse
from pubready.recipes.waterfall import recipe_waterfall
from pubready.spec import spec_from_args


def new_spec(recipe, **kw):
    """Build a figure spec as the plot entry does.

    Args:
        recipe: The recipe name.
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": recipe, "data": "d.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def test_dose_response_reports_an_ic50_and_drops_the_zero_dose():
    rng = np.random.default_rng(1)
    dose = np.repeat([0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10], 3)
    resp = fpl(np.log10(np.maximum(dose, 1e-6)), 5, 100, np.log10(0.3), -1.2) + rng.normal(0, 1, 24)
    res = recipe_dose_response(pd.DataFrame({"dose": dose, "response": resp}),
                               new_spec("dose_response", x="dose", y="response"))
    row = res["stats"].iloc[0]
    assert row["midpoint"] == "IC50"
    assert row["hill"] < 0
    assert row["ec50"] == pytest.approx(0.3, rel=0.1)
    assert any("zero or below" in s for s in res["clean_steps"])
    assert row["n"] == 21


def test_dose_response_needs_four_doses_per_curve():
    df = pd.DataFrame({"dose": [1, 2, 3, 1, 2, 3], "response": [1, 2, 3, 1, 2, 3]})
    with pytest.raises(ValueError, match="at least 4 distinct doses"):
        recipe_dose_response(df, new_spec("dose_response", x="dose", y="response"))


def test_roc_curve_reads_a_marker_whose_low_values_mark_the_case():
    df = pd.DataFrame({"status": ["case"] * 20 + ["control"] * 20,
                       "marker": list(range(1, 21)) + list(range(11, 31))})
    res = recipe_roc_curve(df, new_spec("roc_curve", y="status", x="marker", positive="case"))
    assert res["stats"]["direction"].iloc[0] == ">"
    assert res["stats"]["auc"].iloc[0] > 0.5
    assert res["test_meta"]["case_level"] == "case"


def test_roc_curve_takes_the_second_sorted_level_as_the_case():
    df = pd.DataFrame({"y": [0] * 15 + [1] * 15, "m": list(range(1, 16)) + list(range(6, 21))})
    res = recipe_roc_curve(df, new_spec("roc_curve", y="y", x="m"))
    assert res["test_meta"]["case_level"] == "1"
    assert res["stats"]["direction"].iloc[0] == "<"


def test_roc_curve_needs_two_outcome_levels():
    df = pd.DataFrame({"y": ["a"] * 5 + ["b"] * 5 + ["c"] * 5, "m": range(15)})
    with pytest.raises(ValueError, match="exactly 2 outcome levels"):
        recipe_roc_curve(df, new_spec("roc_curve", y="y", x="m"))


def test_roc_auc_is_the_mann_whitney_probability():
    rng = np.random.default_rng(5)
    a, b = rng.normal(0, 1, 30), rng.normal(0.8, 1, 25)
    df = pd.DataFrame({"y": ["ctrl"] * 30 + ["sick"] * 25, "m": np.concatenate([a, b])})
    res = recipe_roc_curve(df, new_spec("roc_curve", y="y", x="m"))
    u = stats.mannwhitneyu(b, a).statistic
    assert res["stats"]["auc"].iloc[0] == pytest.approx(u / (30 * 25), rel=1e-12)


def test_bland_altman_recovers_a_constant_offset():
    rng = np.random.default_rng(2)
    a = rng.uniform(10, 100, 30)
    df = pd.DataFrame({"a": a + 2 + rng.normal(0, 0.5, 30), "b": a})
    res = recipe_bland_altman(df, new_spec("bland_altman", x="a", y="b"))
    s = res["stats"].iloc[0]
    assert s["bias"] == pytest.approx(2, abs=0.2)
    assert s["loa_low"] < s["bias"] < s["loa_high"]
    assert list(res["qc"]["panels"]) == ["a - b"]


def test_bland_altman_labels_sit_in_the_right_margin():
    # The bias and limit labels sit outside the panel, so they cannot cover a point.
    rng = np.random.default_rng(3)
    a = rng.uniform(10, 100, 40)
    df = pd.DataFrame({"a": a + 2 + rng.normal(0, 1, 40), "b": a})
    res = recipe_bland_altman(df, new_spec("bland_altman", x="a", y="b"))
    fig = res["fig"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax = fig.axes[0]
    labels = [t for t in ax.texts if "SD" in t.get_text() or "Bias" in t.get_text()]
    assert len(labels) == 3
    right = ax.get_window_extent(renderer).x1
    assert all(t.get_window_extent(renderer).x0 > right for t in labels)
    assert all(t.get_window_extent(renderer).x1 <= fig.bbox.x1 for t in labels)


def test_forest_meta_pools_a_difference_measure_by_inverse_variance():
    df = pd.DataFrame({"study": ["s1", "s2", "s3"], "md": [1.0, 2.0, 1.5], "se": [0.5, 0.25, 1.0]})
    res = recipe_forest_meta(df, new_spec("forest_meta", y="md", se="se", label="study"))
    w = 1 / df["se"] ** 2
    fe = res["stats"].set_index("label").loc["Fixed effect"]
    assert res["test_meta"]["measure"] == "difference"
    assert fe["estimate"] == pytest.approx(float((w * df["md"]).sum() / w.sum()), rel=1e-10)
    studies = res["stats"][res["stats"]["row_type"] == "study"]
    assert studies["weight_fixed_pct"].sum() == pytest.approx(100)


def test_forest_meta_picks_the_ratio_scale_for_log_symmetric_cis():
    z = stats.norm.ppf(0.975)
    lr, se = np.log([0.6, 0.8, 0.7, 0.5]), np.array([0.2, 0.3, 0.25, 0.4])
    df = pd.DataFrame({"or": np.exp(lr), "lo": np.exp(lr - z * se), "hi": np.exp(lr + z * se)})
    res = recipe_forest_meta(df, new_spec("forest_meta", y="or", ci_low="lo", ci_high="hi"))
    assert res["test_meta"]["measure"] == "ratio"
    studies = res["stats"][res["stats"]["row_type"] == "study"]
    assert np.allclose(studies["se"], se, rtol=1e-10)


def test_timecourse_without_a_group_runs_a_one_way_rm_anova():
    rng = np.random.default_rng(3)
    df = pd.DataFrame([(f"s{i}", t) for t in (0, 1, 2, 3) for i in range(1, 9)],
                      columns=["id", "time"])
    df["value"] = 10 + 2 * df["time"] + rng.normal(0, 1, len(df))
    res = recipe_timecourse(df, new_spec("timecourse", x="time", y="value", id="id"))
    assert res["test_meta"]["name"] == "one-way repeated-measures ANOVA"
    assert len(res["stats"]) == 1
    assert res["stats"]["p_value"].iloc[0] < 0.001


def test_timecourse_drops_an_incomplete_subject_and_honours_gg():
    rng = np.random.default_rng(4)
    df = pd.DataFrame([(f"s{i}", t) for t in (1, 2, 3) for i in range(1, 11)],
                      columns=["id", "time"])
    df["arm"] = np.where(df["id"].isin([f"s{i}" for i in range(1, 6)]), "a", "b")
    df["value"] = rng.normal(0, 1, len(df))
    df = df[~((df["id"] == "s1") & (df["time"] == 3))]
    res = recipe_timecourse(df, new_spec("timecourse", x="time", y="value", id="id",
                                         group="arm", test="gg"))
    assert res["test_meta"]["n_subjects"] == 9
    assert res["test_meta"]["correction_applied"] is True
    assert res["stats"].set_index("effect").loc["time", "correction"] == "Greenhouse-Geisser"


def test_timecourse_refuses_two_values_for_one_subject_and_time():
    df = pd.DataFrame({"id": ["a", "a", "a", "b", "b"], "time": [1, 1, 2, 1, 2], "value": range(5)})
    with pytest.raises(ValueError, match="one value per subject and time"):
        recipe_timecourse(df, new_spec("timecourse", x="time", y="value", id="id"))


def test_timecourse_refuses_a_subject_in_two_groups():
    df = pd.DataFrame({"id": ["a", "a", "b", "b"], "time": [1, 2, 1, 2], "value": range(4),
                       "arm": ["x", "y", "x", "x"]})
    with pytest.raises(ValueError, match="one group only"):
        recipe_timecourse(df, new_spec("timecourse", x="time", y="value", id="id", group="arm"))


def test_waterfall_reads_one_change_per_patient_and_sorts_the_bars():
    df = pd.DataFrame({"patient": ["p1", "p2", "p3", "p4"], "best": [-45, 25, -10, 5]})
    res = recipe_waterfall(df, new_spec("waterfall", y="best", id="patient"))
    assert list(res["stats"]["patient"]) == ["p2", "p4", "p3", "p1"]
    assert list(res["stats"]["category_from_change"]) == [
        "progressive disease", "stable disease", "stable disease", "partial response"]


def test_waterfall_without_time_refuses_two_rows_for_one_patient():
    df = pd.DataFrame({"patient": ["p1", "p1"], "best": [-45, 25]})
    with pytest.raises(ValueError, match="one row per patient"):
        recipe_waterfall(df, new_spec("waterfall", y="best", id="patient"))


def test_swimmer_sorts_the_courses_and_counts_ongoing_and_markers():
    df = pd.DataFrame({"pt": ["a", "b", "c"], "months": [3, 10, 6], "ended": [1, 0, 1],
                       "response": [np.nan, 2, 4]})
    res = recipe_swimmer(df, new_spec("swimmer", id="pt", time="months", event="ended",
                                      markers="response"))
    assert list(res["stats"]["patient"]) == ["b", "c", "a"]
    assert res["test_meta"]["n_ongoing"] == 1
    assert res["test_meta"]["markers"]["response"] == 2


def test_swimmer_refuses_a_status_not_coded_one_or_zero():
    df = pd.DataFrame({"pt": ["a", "b"], "months": [3, 10], "ended": [1, 2]})
    with pytest.raises(ValueError, match="1 = the course ended"):
        recipe_swimmer(df, new_spec("swimmer", id="pt", time="months", event="ended"))
