"""survival_km: the number at risk counts time >= t at breaks that both engines
compute in the same way. The R suite (test-survival-km.R) expects the same
numbers."""
from types import SimpleNamespace

import pandas as pd
import pytest

from pubready.recipes.survival_km import km_break_step, km_breaks, recipe_survival_km
from pubready.spec import spec_from_args


def km_data():
    """Two arms with events that fall on the breaks.

    Returns:
        pandas.DataFrame: One row per subject.
    """
    return pd.DataFrame({
        "arm": ["A"] * 9 + ["B"] * 9,
        "months": [2, 5, 5, 8, 10, 10, 14, 18, 21.8, 3, 6, 10, 12, 15, 17, 19, 20, 20],
        "event": [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1],
    })


def test_the_number_at_risk_counts_the_subjects_at_or_after_each_break():
    spec = spec_from_args(SimpleNamespace(recipe="survival_km", data="d.csv", x="arm", y=None,
                                          test="auto", time="months", event="event"))
    res = recipe_survival_km(km_data(), spec)
    ar = res["test_meta"]["at_risk"]
    assert ar["times"] == pytest.approx([0, 5, 10, 15, 20])
    assert ar["counts"]["A"] == [9, 8, 5, 2, 1]
    assert ar["counts"]["B"] == [9, 8, 7, 5, 2]
    # The x axis reaches the last time, so the last event is on the figure.
    assert res["fig"].axes[0].get_xlim()[1] >= 21.8


def test_the_x_breaks_follow_the_1_2_5_step_that_gives_at_most_six_intervals():
    assert km_break_step(21.8) == 5
    assert km_break_step(60) == 10
    assert km_break_step(1) == pytest.approx(0.2)
    assert km_breaks(21.8, 5) == pytest.approx([0, 5, 10, 15, 20])
