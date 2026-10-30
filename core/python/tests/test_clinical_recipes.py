"""Behaviour of the cumulative_incidence, calibration and sankey recipes on small
synthetic inputs. The reference values of the tied-time and tied-decile cases
come from the R recipes, which call cmprsk::cuminc and match
ResourceSelection::hoslem.test (see test-clinical-recipes.R)."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.calibration import quantile7, recipe_calibration
from pubready.recipes.cumulative_incidence import (cif_curve, cif_value_at,
                                                   recipe_cumulative_incidence)
from pubready.recipes.sankey import recipe_sankey
from pubready.spec import spec_from_args


def clin_spec(recipe, **kw):
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


def cif_tiny():
    """Three groups of twelve subjects with tied times and two causes.

    Returns:
        pandas.DataFrame: The data.
    """
    return pd.DataFrame({
        "time": [2, 3, 3, 5, 6, 6, 8, 9, 11, 12, 12, 14, 1, 3, 4, 4, 6, 7, 9, 9, 10, 13, 15, 15,
                 2, 2, 5, 7, 8, 8, 10, 11, 12, 14, 14, 16],
        "status": [1, 0, 2, 1, 1, 0, 2, 1, 0, 1, 2, 0, 2, 1, 1, 2, 0, 1, 1, 0, 2, 1, 0, 1,
                   1, 1, 0, 2, 1, 2, 0, 1, 1, 2, 0, 0],
        "group": ["A"] * 12 + ["B"] * 12 + ["C"] * 12})


def test_cumulative_incidence_matches_cmprsk_on_tied_times():
    """Gray's test and the CIF match cmprsk on tied times with three groups."""
    res = recipe_cumulative_incidence(cif_tiny(), clin_spec(
        "cumulative_incidence", time="time", event="status", x="group"))
    st = res["stats"]
    r1, r2 = st[st["cause"] == "1"], st[st["cause"] == "2"]
    # Gray's test and the CIF at 7.5 and 12.5 as cmprsk::cuminc reports them.
    assert r1["gray_statistic"].iloc[0] == pytest.approx(0.09669301500815727, rel=1e-9)
    assert r1["gray_p"].iloc[0] == pytest.approx(0.9528035762814548, rel=1e-9)
    assert r2["gray_statistic"].iloc[0] == pytest.approx(0.053571669025497974, rel=1e-9)
    assert list(r1["cif_7.5"]) == pytest.approx([0.2685185185185185, 0.26190476190476186,
                                                 0.16666666666666666], rel=1e-12)
    assert list(r2["cif_12.5"]) == pytest.approx([0.3353909465020576, 0.2857142857142857,
                                                  0.1851851851851852], rel=1e-12)
    assert list(st.columns[4:7]) == ["cif_5", "cif_7.5", "cif_12.5"]
    assert r1["gray_df"].iloc[0] == 2


def test_cumulative_incidence_puts_the_named_cause_first():
    """The cause that --cause names comes first, and an unknown code fails."""
    res = recipe_cumulative_incidence(cif_tiny(), clin_spec(
        "cumulative_incidence", time="time", event="status", x="group", cause="2"))
    assert res["stats"]["cause"].iloc[0] == "2"
    assert res["test_meta"]["causes"] == ["2", "1"]
    with pytest.raises(ValueError, match="is not an event code"):
        recipe_cumulative_incidence(cif_tiny(), clin_spec(
            "cumulative_incidence", time="time", event="status", cause="7"))


def test_cumulative_incidence_of_one_cause_is_one_minus_kaplan_meier():
    """With one cause and one group, the CIF equals 1 - Kaplan-Meier."""
    t = np.array([1, 2, 2, 3, 5, 6, 7, 9], dtype=float)
    e = np.array([1, 0, 1, 1, 0, 1, 0, 1])
    xs, fs = cif_curve(t, e, e)
    km = (7 / 8) * (6 / 7) * (4 / 5) * (2 / 3)
    assert cif_value_at(6, xs, fs) == pytest.approx(1 - km, rel=1e-12)
    assert np.isnan(cif_value_at(10, xs, fs))
    res = recipe_cumulative_incidence(pd.DataFrame({"time": t, "status": e}), clin_spec(
        "cumulative_incidence", time="time", event="status"))
    assert np.isnan(res["stats"]["gray_p"].iloc[0])


def test_quantile7_follows_r_type_7():
    """The type-7 quantiles follow the arithmetic of R."""
    x = np.arange(1, 11, dtype=float)
    probs = np.arange(11) * 0.1
    assert quantile7(x, probs) == pytest.approx(1 + 0.9 * np.arange(11), rel=1e-12)
    assert quantile7([0.1, 0.1, 0.1, 0.5], [0.5]) == pytest.approx([0.1])


def test_calibration_groups_tied_predictions_as_hoslem_test():
    """Tied deciles merge, and an empty group is left out, as in the R engine."""
    p = np.repeat([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [3, 4, 3, 4, 3, 4, 2, 3, 4])
    y = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
    res = recipe_calibration(pd.DataFrame({"y": y, "p": p}), clin_spec("calibration", y="y", x="p"))
    row = res["stats"].iloc[0]
    assert row["hl_statistic"] == pytest.approx(3.812047479810638, rel=1e-10)
    assert row["hl_p"] == pytest.approx(0.7020938638256268, rel=1e-9)
    assert row["hl_groups"] == 8 and row["hl_df"] == 6
    assert row["slope"] == pytest.approx(0.6753078755364041, rel=1e-8)
    assert row["intercept"] == pytest.approx(0.38448260942044593, rel=1e-8)
    assert row["intercept_ci_low"] == pytest.approx(-0.4570281394993177, rel=1e-8)
    assert row["case_level"] == "1"


def test_calibration_recovers_the_ideal_slope_for_true_risks():
    """True risks give a slope near 1 and an intercept near 0."""
    rng = np.random.default_rng(7)
    p = 1 / (1 + np.exp(-rng.normal(-0.5, 1.2, 4000)))
    y = rng.binomial(1, p)
    d = pd.DataFrame({"outcome": np.where(y == 1, "event", "none"), "risk": p})
    row = recipe_calibration(d, clin_spec("calibration", y="outcome", x="risk",
                                          positive="event"))["stats"].iloc[0]
    assert row["slope_ci_low"] < 1 < row["slope_ci_high"]
    assert row["intercept_ci_low"] < 0 < row["intercept_ci_high"]


def test_calibration_rejects_values_outside_the_unit_interval():
    """A prediction outside [0, 1] is an error."""
    d = pd.DataFrame({"y": [0, 1, 0, 1], "p": [0.2, 1.3, 0.1, 0.8]})
    with pytest.raises(ValueError, match=r"outside \[0, 1\]"):
        recipe_calibration(d, clin_spec("calibration", y="y", x="p"))


def test_sankey_counts_the_transitions_of_patients_present_at_both_visits():
    """Only patients present at both visits make a transition."""
    d = pd.DataFrame({"pid": [1, 1, 1, 2, 2, 3, 3, 3, 4, 4],
                      "visit": [8, 16, 24, 8, 24, 8, 16, 24, 8, 16],
                      "state": ["SD", "PR", "PR", "SD", "PD", "PR", "PR", "CR", "SD", "SD"]})
    res = recipe_sankey(d, clin_spec("sankey", id="pid", x="visit", y="state"))
    assert res["test_meta"]["visits"] == ["8", "16", "24"]
    assert res["test_meta"]["states"] == ["SD", "PR", "PD", "CR"]
    st = res["stats"]
    assert st.loc[st["from_visit"] == "8", "n"].sum() == 3
    assert res["test_meta"]["no_test_reason"] == "more than two visits"


def test_sankey_runs_bowker_and_drops_the_empty_pairs():
    """Bowker's test leaves out a pair of cells that holds no patient."""
    frm = ["a"] * 10 + ["b"] * 8 + ["c"] * 6
    to = ["a"] * 5 + ["b"] * 4 + ["c"] + ["b"] * 6 + ["a"] * 2 + ["c"] * 6
    d = pd.DataFrame({"id": list(range(24)) * 2, "visit": ["pre"] * 24 + ["post"] * 24,
                      "state": frm + to})
    meta = recipe_sankey(d, clin_spec("sankey", id="id", x="visit", y="state"))["test_meta"]
    assert meta["visits"] == ["post", "pre"]
    assert meta["name"] == "Bowker's test of symmetry"
    assert meta["empty_pairs"] == 1 and meta["df"] == 2
    stat = (2 - 4) ** 2 / 6 + (0 - 1) ** 2 / 1
    assert meta["statistic"] == pytest.approx(stat, rel=1e-12)
    assert meta["p_value"] == pytest.approx(stats.chi2.sf(stat, 2), rel=1e-12)


def test_sankey_uses_mcnemar_for_two_states():
    """Two states give McNemar's test without continuity correction."""
    frm = ["no"] * 12 + ["yes"] * 8
    to = ["no"] * 7 + ["yes"] * 5 + ["no"] * 2 + ["yes"] * 6
    d = pd.DataFrame({"id": list(range(20)) * 2, "visit": [0] * 20 + [1] * 20, "state": frm + to})
    meta = recipe_sankey(d, clin_spec("sankey", id="id", x="visit", y="state"))["test_meta"]
    assert meta["name"] == "McNemar's test"
    assert meta["statistic"] == pytest.approx((5 - 2) ** 2 / 7, rel=1e-12)


def test_sankey_needs_one_row_per_patient_and_visit():
    """Two rows for one patient at one visit are an error."""
    d = pd.DataFrame({"id": [1, 1, 2], "visit": [1, 1, 2], "state": ["a", "b", "a"]})
    with pytest.raises(ValueError, match="one row per patient and visit"):
        recipe_sankey(d, clin_spec("sankey", id="id", x="visit", y="state"))


def sk_visits():
    """Four patients over three visits, the data of the R suite.

    Returns:
        pandas.DataFrame: The data.
    """
    return pd.DataFrame({"pid": [1, 1, 1, 2, 2, 3, 3, 3, 4, 4],
                         "visit": [8, 16, 24, 8, 24, 8, 16, 24, 8, 16],
                         "state": ["SD", "PR", "PR", "SD", "PD", "PR", "PR", "CR", "SD", "SD"]})


def test_sankey_levels_set_the_state_order_and_leave_out_an_unused_state():
    """--levels orders the states, as in the R engine."""
    res = recipe_sankey(sk_visits(), clin_spec("sankey", id="pid", x="visit", y="state",
                                               levels="CR, PR,SD,PD,NE"))
    assert res["test_meta"]["states"] == ["CR", "PR", "SD", "PD"]
    assert "states in the order set by the analyst (CR, PR, SD, PD)" in res["methods"]


def test_sankey_levels_stop_when_a_state_of_the_data_is_not_listed():
    """A state that --levels leaves out is an error with the R message."""
    with pytest.raises(ValueError) as err:
        recipe_sankey(sk_visits(), clin_spec("sankey", id="pid", x="visit", y="state",
                                             levels="CR,PR,SD"))
    assert str(err.value) == "--levels does not list the state(s) 'PD' of the data"


def test_calibration_reports_the_slope_as_not_estimable_for_constant_predictions():
    """Constant predictions: slope and HL test not estimable, as in R (F78)."""
    from scipy.special import logit
    y = np.array([1, 0, 0] * 20)
    res = recipe_calibration(pd.DataFrame({"y": y, "p": 0.3}), clin_spec("calibration", y="y", x="p"))
    st = res["stats"]
    assert np.isnan(st["slope"].iloc[0]) and np.isnan(st["hl_p"].iloc[0])
    assert res["test_meta"]["models"]["p"]["slope_estimable"] is False
    assert "slope and Hosmer-Lemeshow test not estimable" in res["methods"]
    assert st["intercept"].iloc[0] == pytest.approx(logit(1 / 3) - logit(0.3), rel=1e-6)
