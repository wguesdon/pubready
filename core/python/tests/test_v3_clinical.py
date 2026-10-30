"""The v3 fixes of the clinical recipes. The R suite (test-v3-clinical.R) runs
the same cases and expects the same results."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.cox_forest import collinear_columns, recipe_cox_forest
from pubready.recipes.cumulative_incidence import recipe_cumulative_incidence
from pubready.recipes.survival_km import (event_codes, km_break_step, km_median,
                                          recipe_survival_km)
from pubready.recipes.swimmer import recipe_swimmer
from pubready.spec import spec_from_args


def v3c_spec(recipe, **kw):
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


def km_space():
    """Two arms whose column names hold a space.

    Returns:
        pandas.DataFrame: One row per subject.
    """
    return pd.DataFrame({
        "follow up months": [2, 5, 5, 8, 10, 10, 14, 18, 21.8, 3, 6, 10, 12, 15, 17, 19, 20, 20],
        "death event": [1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1],
        "treatment arm": ["A"] * 9 + ["B"] * 9,
    })


def test_survival_km_reads_columns_whose_names_hold_a_space():
    """F44: the R engine stopped on such names; both engines draw them."""
    res = recipe_survival_km(km_space(), v3c_spec("survival_km", time="follow up months",
                                                  event="death event", x="treatment arm"))
    assert list(res["stats"]["group"]) == ["A", "B"]
    assert list(res["stats"]["events"]) == [7, 6]
    plain = km_space().set_axis(["t", "e", "g"], axis=1)
    ref = recipe_survival_km(plain, v3c_spec("survival_km", time="t", event="e", x="g"))
    assert res["test_meta"]["log_rank"]["chisq"] == pytest.approx(
        ref["test_meta"]["log_rank"]["chisq"], rel=1e-12)
    assert "groups of treatment arm" in res["methods"]


def test_an_event_column_holds_0_1_or_a_logical_else_the_message_lists_its_values():
    """F56: a text or a decimal event code stops the run and names the values."""
    assert list(event_codes(pd.Series([True, False, np.nan], dtype=object), "e").fillna(-1)) == [1, 0, -1]
    assert list(event_codes(pd.Series(["1", " 0", np.nan, ""]), "e").fillna(-1)) == [1, 0, -1, -1]
    assert list(event_codes(pd.Series([1.0, 0.0, 1.0, np.nan]), "e").fillna(-1)) == [1, 0, 1, -1]
    with pytest.raises(ValueError, match="event column 'status' must be coded 1 = event, "
                                         "0 = censored, but it holds 'alive', 'dead'"):
        event_codes(pd.Series(["dead", "alive", "dead"]), "status")
    with pytest.raises(ValueError, match="but it holds '1.5', '2', '10'"):
        event_codes(pd.Series([1, 1.5, 0, 2, 10]), "e")
    with pytest.raises(ValueError, match="'6', '7' and 2 more"):
        event_codes(pd.Series([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), "e")
    d = pd.DataFrame({"time": [3, 5, 8, 9], "event": ["dead", "alive", "dead", "alive"],
                      "arm": ["a", "a", "b", "b"]})
    with pytest.raises(ValueError, match="holds 'alive', 'dead'"):
        recipe_survival_km(d, v3c_spec("survival_km", time="time", event="event", x="arm"))
    with pytest.raises(ValueError, match="holds 'alive', 'dead'"):
        recipe_cox_forest(d, v3c_spec("cox_forest", time="time", event="event", covariates="arm"))
    s = pd.DataFrame({"patient": ["P1", "P2"], "months": [4, 7], "ended": ["yes", "no"]})
    with pytest.raises(ValueError, match="must be coded 1 = the course ended, 0 = ongoing, "
                                         "but it holds 'no', 'yes'"):
        recipe_swimmer(s, v3c_spec("swimmer", id="patient", time="months", event="ended"))


def test_the_median_survival_takes_the_midpoint_of_a_step_at_exactly_one_half():
    """F67: lifelines took the start of the step; survfit takes the midpoint."""
    d = pd.DataFrame({"time": list(range(1, 11)) + list(range(2, 21, 2)) + [1, 2, 2.5, 3],
                      "event": [1] * 20 + [1, 1, 0, 1],
                      "arm": ["a"] * 10 + ["b"] * 10 + ["c"] * 4})
    res = recipe_survival_km(d, v3c_spec("survival_km", time="time", event="event", x="arm"))
    # The medians that survival::survfit reports (see test-v3-clinical.R).
    assert list(res["stats"]["median_survival"]) == pytest.approx([5.5, 11, 2.5])
    assert np.isnan(km_median(np.array([0.0, 1, 2, 3]), np.array([1.0, 0.8, 0.6, 0.55])))
    assert km_median(np.array([0.0, 1, 2]), np.array([1.0, 0.75, 0.5])) == 2


def test_survival_km_names_one_group_and_draws_times_above_5e6():
    """F92: one group gave chi-square(0) in Python; times above 5e6 stopped both engines."""
    d = pd.DataFrame({"time": [3, 5, 8, 9], "event": [1, 0, 1, 1], "arm": ["control"] * 4})
    with pytest.raises(ValueError, match=r"column 'arm' holds 1 group \('control', n = 4\)"):
        recipe_survival_km(d, v3c_spec("survival_km", time="time", event="event", x="arm"))
    assert km_break_step(1.038e8) == 2e7
    assert km_break_step(5.5e12) == 1e12
    d2 = pd.DataFrame({"time": [v * 1e7 for v in (1, 3, 5, 8, 9, 10)],
                       "event": [1, 0, 1, 1, 1, 0], "arm": ["a"] * 3 + ["b"] * 3})
    res = recipe_survival_km(d2, v3c_spec("survival_km", time="time", event="event", x="arm"))
    assert res["test_meta"]["at_risk"]["times"] == pytest.approx([v * 1e7 for v in (0, 2, 4, 6, 8, 10)])


def test_survival_cox_and_cumulative_incidence_share_one_level_order():
    """F64: numeric groups sort as numbers and text levels in code-point order."""
    d = pd.DataFrame({"time": [5, 8, 12, 3, 9, 14, 7, 10, 15, 4, 6, 11],
                      "event": [1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0],
                      "dose": [80] * 4 + [160] * 4 + [235] * 4,
                      "inherit": ["X-linked", "autosomal"] * 6,
                      "status": [1, 2, 0, 1, 1, 0, 2, 1, 0, 1, 2, 1]})
    km = recipe_survival_km(d, v3c_spec("survival_km", time="time", event="event", x="dose"))
    assert list(km["stats"]["group"]) == ["80", "160", "235"]
    cx = recipe_cox_forest(d, v3c_spec("cox_forest", time="time", event="event",
                                       covariates="inherit"))
    assert list(cx["stats"]["term"]) == ["inheritautosomal"]
    ci = recipe_cumulative_incidence(d, v3c_spec("cumulative_incidence", time="time",
                                                 event="status", x="dose"))
    assert ci["test_meta"]["groups"] == ["80", "160", "235"]
    assert list(pd.unique(ci["stats"]["group"])) == ["80", "160", "235"]


def cox_data():
    """Sixteen subjects with an arm and an age.

    Returns:
        pandas.DataFrame: One row per subject.
    """
    return pd.DataFrame({
        "time": [5, 8, 12, 3, 9, 14, 7, 10, 15, 4, 6, 11, 13, 2, 16, 9],
        "event": [1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1],
        "arm": ["control", "treated"] * 8,
        "age": [61, 55, 70, 48, 66, 59, 52, 73, 45, 68, 57, 63, 50, 71, 58, 64]})


def test_a_text_cell_in_a_numeric_cox_covariate_becomes_nan_and_is_logged():
    """F42: one "45 y" made age a factor of many levels."""
    d = cox_data()
    d["age"] = d["age"].astype(str)
    d.loc[2, "age"] = "45 y"
    res = recipe_cox_forest(d, v3c_spec("cox_forest", time="time", event="event",
                                        covariates="arm,age"))
    assert list(res["stats"]["term"]) == ["armtreated", "age"]
    assert res["test_meta"]["n"] == 15
    assert res["clean_steps"][0] == ("Coerced column 'age' to numeric; 1 value(s) became NA, "
                                     "for example '45 y'.")
    assert res["clean_steps"][1].startswith("Dropped 1 row")


def cox_run(data, covs):
    """Run cox_forest on the time and event columns of cox_data.

    Args:
        data: The table.
        covs: The covariates, comma separated.

    Returns:
        dict: The bundle payload.
    """
    return recipe_cox_forest(data, v3c_spec("cox_forest", time="time", event="event",
                                            covariates=covs))


def test_cox_forest_names_a_constant_or_collinear_covariate_and_no_event():
    """F91: these inputs stopped R, gave NA rows or a traceback in Python."""
    d = cox_data()
    d["site"] = "A"
    d["dose"] = 10
    d["age_copy"] = 2 * d["age"] + 1
    d["arm2"] = d["arm"]
    with pytest.raises(ValueError, match=r"covariate 'site' holds one value \('A'\) in the 16 "
                                         r"subject\(s\) kept"):
        cox_run(d, "arm,age,site")
    with pytest.raises(ValueError, match=r"covariate 'dose' holds one value \('10'\)"):
        cox_run(d, "arm,dose")
    with pytest.raises(ValueError, match="covariate 'age_copy' is collinear with the covariates "
                                         "before it"):
        cox_run(d, "arm,age,age_copy")
    with pytest.raises(ValueError, match="covariate 'arm2' is collinear"):
        cox_run(d, "arm,age,arm2")
    d["event"] = 0
    with pytest.raises(ValueError, match=r"column 'event' holds no 1 among the 16 subject\(s\) kept"):
        cox_run(d, "arm,age")
    x = np.array([[1, 2], [2, 1], [3, 4], [5, 4]], dtype=float)
    assert collinear_columns(np.column_stack([x, x[:, 0] + 2 * x[:, 1]])) == [2]
    assert collinear_columns(np.column_stack([x[:, 0], 3 - x[:, 0], x[:, 1]])) == [1]
    assert collinear_columns(x) == []
