"""Advice routing, scientist decisions, conflicts and deterministic output."""

import os
from pathlib import Path

import pandas as pd
import pytest

from pubready.advise import advise, collect_facts, load_tree, render_result

BENCH = Path(os.environ.get("PUBREADY_BENCHMARKS", "/opt/pubready/benchmarks"))
CASES = [
    ("sleep", dict(x="group", y="extra"), ["aim=compare", "design=paired"], "paired_compare", 0),
    ("sleep", dict(x="group", y="extra"), ["aim=compare", "design=independent"], "two_group_compare", 0),
    ("sleep", dict(x="group", y="extra"), ["aim=compare"], None, 3),
    ("plant_growth", dict(x="group", y="weight"), ["aim=compare", "design=independent"], "multi_group_compare", 0),
    ("insect_sprays", dict(x="spray", y="count", test="kruskal"), ["aim=compare", "design=independent"], "multi_group_compare", 0),
    ("insect_sprays", dict(x="spray", y="count", test="wilcoxon"), ["aim=compare", "design=independent"], "multi_group_compare", 2),
    ("tooth_growth", dict(x="dose", y="len", fill="supp"), ["aim=compare", "design=independent"], "factorial_anova", 0),
    ("chick_weight", dict(x="Time", y="weight", id="Chick", group="Diet"), ["aim=compare", "design=repeated"], "timecourse", 0),
    ("chick_weight", dict(x="Time", y="weight", fill="Diet", engine="python"), ["aim=compare", "design=independent"], "factorial_anova", 0),
    ("orthodont", dict(x="age", y="distance", id="Subject", group="Sex"), ["aim=compare", "design=repeated"], "timecourse", 0),
    ("anorexia_ft", dict(x="condition", y="weight_lbs"), ["aim=compare", "design=paired"], "paired_compare", 0),
    ("chickwts", dict(x="feed", y="weight", ref="casein"), ["aim=compare", "design=independent", "presentation=estimation"], "estimation_plot", 0),
    ("chickwts", dict(x="feed", y="weight"), ["aim=compare", "design=independent", "presentation=estimation"], "estimation_plot", 3),
    ("iris", dict(x="Species"), ["aim=compare", "design=independent"], "feature_panel", 0),
    ("titanic", dict(x="Class", y="Survived"), ["aim=compare", "category_design=independent"], "proportions", 0),
    ("pbcseq_edema", dict(id="id", x="year", y="edema"), ["aim=compare", "outcome_type=category", "category_design=repeated"], "sankey", 0),
    ("mtcars", dict(x="wt", y="mpg"), ["aim=association"], "correlation", 0),
    ("pefr", dict(x="wright_1", y="mini_1"), ["aim=agreement"], "bland_altman", 0),
    ("asah", dict(x="s100b,ndka", y="outcome", positive="Poor"), ["aim=discrimination"], "roc_curve", 0),
    ("asah", dict(x="s100b,ndka", y="outcome", positive="Bad"), ["aim=discrimination"], "roc_curve", 2),
    ("pima_te", dict(x="predicted", y="type", positive="Yes"), ["aim=calibration"], "calibration", 0),
    ("bcg", dict(y="risk_ratio", ci_low="ci_low", ci_high="ci_high", label="trial", measure="ratio"), ["aim=meta_analysis"], "forest_meta", 0),
    ("lung", dict(time="time", event="event", x="sex"), ["aim=time_to_event", "tte_aim=compare_groups"], "survival_km", 0),
    ("veteran", dict(time="time", event="status", covariates="trt,celltype,karno,age"), ["aim=time_to_event", "tte_aim=adjusted_hazards"], "cox_forest", 0),
    ("mgus2", dict(time="etime", event="event", x="sex"), ["aim=time_to_event"], None, 3),
    ("mgus2", dict(time="etime", event="event", x="sex", cause="1"), ["aim=time_to_event", "competing_events=present"], "cumulative_incidence", 0),
    ("mgus2", dict(time="etime", event="event", x="sex"), ["aim=time_to_event", "competing_events=absent"], None, 4),
    ("dnase", dict(x="conc", y="density", group="Run"), ["aim=curve", "curve_kind=dose_response"], "dose_response", 0),
    ("puromycin", dict(x="conc", y="rate", group="state"), ["aim=curve", "curve_kind=michaelis_menten"], "enzyme_kinetics", 0),
    ("dnase_standard_curve", dict(x="conc", y="density", id="sample"), ["aim=curve"], None, 3),
]


@pytest.mark.parametrize("dataset,flags,answers,recipe,code", CASES)
def test_benchmark_routes(dataset, flags, answers, recipe, code):
    """Check the approved benchmark routes and identical repeat output."""
    options = dict(data=str(BENCH / "data" / f"{dataset}.csv"), **flags)
    result = advise(options, answers)
    assert result["exit_code"] == code, result
    assert result["recipe"] == recipe
    assert render_result(result) == render_result(advise(options, answers))
    if code:
        assert result["command"] is None


@pytest.mark.parametrize("design", ["stratified", "covariate_adjusted", "rate_of_change", "control_only"])
def test_unsupported_designs(design):
    """Never send an unsupported design to a pooled comparison."""
    result = advise(dict(data=str(BENCH / "data/sleep.csv"), x="group", y="extra"),
                    ["aim=compare", f"design={design}"])
    assert result["exit_code"] == 5
    assert result["command"] is None


def test_ambiguous_facet():
    """Ask before treating a facet as a feature or a third factor."""
    df = pd.DataFrame({"x": ["a", "b"] * 10, "y": range(20),
                       "fill": ["u"] * 10 + ["v"] * 10, "facet": ["f", "g"] * 10})
    flags = dict(x="x", y="y", facet="facet")
    answers = ["aim=compare", "design=independent"]
    missing = advise(flags, answers, frame=df)
    assert missing["open_questions"][0]["id"] == "n_outcomes"
    panel = advise(flags, answers + ["n_outcomes=panel"], frame=df)
    assert panel["recipe"] == "feature_panel"
    factorial = advise(dict(flags, fill="fill"), answers + ["n_outcomes=one"], frame=df)
    assert factorial["recipe"] == "factorial_anova"


def test_explicit_answers_and_unused_flags(tmp_path):
    """Reject duplicate conflicts, invalid keys and unused recipe flags."""
    flags = dict(data=str(BENCH / "data/sleep.csv"), x="group", y="extra")
    assert advise(flags, ["aim=compare", "aim=curve"])["exit_code"] == 2
    assert advise(flags, ["aim=unknown"])["exit_code"] == 2
    assert advise(flags, ["unknown=x"])["exit_code"] == 2
    path = tmp_path / "answers.yaml"
    path.write_text("aim: compare\naim: curve\n")
    assert advise(flags, filename=path)["exit_code"] == 2
    path.write_text("aim: compare\ndesign: independent\n")
    assert advise(flags, filename=path)["exit_code"] == 0
    bad = advise(dict(flags, group="ID"), ["aim=compare", "design=independent"])
    assert bad["exit_code"] == 2


def test_numeric_and_time_orders():
    """Keep the shared numeric order and the paired time rule."""
    df = pd.DataFrame({"group": ["week 16", "week 8"] * 10,
                       "value": range(20), "id": [i for i in range(10) for _ in range(2)]})
    result = advise(dict(x="group", y="value", id="id"), ["aim=compare", "design=paired"], frame=df)
    assert result["facts"]["x.plot_order"] == ["week 8", "week 16"]
    df["group"] = ["16", "8"] * 10
    facts = collect_facts(df, dict(x="group", y="value"))
    assert facts["x.levels"] == ["8", "16"]


def test_case_control_cause_and_forced_art():
    """Require domain choices and enforce the Python ART restriction."""
    roc = advise(dict(data=str(BENCH / "data/asah.csv"), x="s100b", y="outcome"),
                 ["aim=discrimination"])
    assert roc["exit_code"] == 3
    assert roc["open_questions"][0]["id"] == "positive"
    art = advise(dict(data=str(BENCH / "data/tooth_growth.csv"), x="dose", y="len", fill="supp",
                      engine="python", test="art"), ["aim=compare", "design=independent"])
    assert art["exit_code"] == 4


@pytest.mark.parametrize("recipe", list(load_tree()[0]["questions"]["describe_kind"]["answers"]))
def test_every_descriptive_route(recipe):
    """Reach every descriptive leaf, without requiring a statistical fit."""
    result = advise(dict(data=str(BENCH / "data/iris.csv")), ["aim=describe", f"describe_kind={recipe}"])
    assert result["recipe"] == recipe
    assert result["exit_code"] in (0, 3)


def test_paired_flag_cannot_emit_row_order_pairing():
    """An explicit paired flag uses the subject-id route and rejects a conflict."""
    flags = dict(data=str(BENCH / "data/sleep.csv"), x="group", y="extra", paired=True)
    result = advise(flags, ["aim=compare"])
    assert result["recipe"] == "paired_compare"
    assert "--paired" not in result["command"]
    assert advise(flags, ["aim=compare", "design=independent"])["exit_code"] == 2


def test_binary_event_normalization_and_prediction_range():
    """Normalize logical events and reject impossible predicted probabilities."""
    df = pd.DataFrame({"time": range(8), "event": [True, False] * 4,
                       "group": ["a", "b"] * 4, "prediction": [1.2] * 8})
    result = advise(dict(time="time", event="event", x="group"),
                    ["aim=time_to_event", "tte_aim=compare_groups"], frame=df)
    assert result["exit_code"] == 0
    bad = advise(dict(x="prediction", y="event", positive="TRUE"), ["aim=calibration"], frame=df)
    assert bad["exit_code"] == 4
