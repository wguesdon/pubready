import numpy as np
import pandas as pd
import pytest

from pubready.recipes.cox_forest import forest_log_ticks, forest_stars, recipe_cox_forest


@pytest.mark.parametrize("p, stars", [
    (1e-9, "***"), (5e-4, "***"), (5e-3, "**"), (0.02, "*"), (0.05, ""), (0.3, ""),
])
def test_forest_stars_stop_at_three(p, stars):
    # The forest plot grades the stars as survminer::ggforest did, in both engines.
    assert forest_stars(p) == stars


@pytest.mark.parametrize("lo, hi, ticks", [
    # The 1-2-5 sequence, extended to 2 so that the CIs above 1 are covered.
    (0.34, 1.94, [0.5, 1, 2]),
    (0.44, 3.54, [0.5, 1, 2]),
    # A narrow range takes finer steps rather than a nearly empty axis.
    (0.8, 1.3, [0.7, 1, 1.5]),
    (0.99, 1.02, [0.8, 1, 1.25]),
    # A range over many decades keeps the powers of ten only.
    (0.001, 500, [0.001, 0.01, 0.1, 1, 10, 100]),
])
def test_forest_log_ticks_cover_both_sides_of_one(lo, hi, ticks):
    # The R engine has the same cases in test-cox-forest.R.
    assert forest_log_ticks(lo, hi) == pytest.approx(ticks)


def test_cox_forest_shows_a_reference_row_per_factor():
    rng = np.random.RandomState(1)
    n = 120
    df = pd.DataFrame({
        "time": rng.exponential(10, n), "event": rng.binomial(1, 0.7, n),
        "arm": rng.choice(["control", "treated"], n), "age": rng.normal(60, 8, n),
        "stage": rng.choice(["I", "II", "III"], n),
    })
    spec = {"data": {"time": "time", "event": "event", "covariates": ["arm", "age", "stage"]},
            "appearance": {}}
    res = recipe_cox_forest(df, spec)
    assert list(res["stats"]["term"]) == ["armtreated", "age", "stageII", "stageIII"]
    texts = [t.get_text() for t in res["fig"].texts]
    for label in ("arm: control", "arm: treated", "age", "stage: I", "stage: III"):
        assert label in texts
    table = [t.get_text() for t in res["fig"].axes[1].texts]
    assert table.count("reference") == 2
    assert "HR (95% CI)" in table and "P" in table


def test_the_cox_subtitle_warns_when_the_proportional_hazards_check_fails():
    from pubready.recipes.cox_forest import cox_subtitle
    assert cox_subtitle(137, 128, 1e-5) == "137 subjects, 128 events; PH assumption violated (P < 0.001)"
    assert cox_subtitle(227, 164, 0.22) == "227 subjects, 164 events"
    assert cox_subtitle(227, 164, float("nan")) == "227 subjects, 164 events"
