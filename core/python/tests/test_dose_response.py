"""dose_response: the "not determined" rule. The R suite (test-dose-response.R)
runs the same data and expects the same result."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.dose_response import fpl, is_flat, recipe_dose_response, status_reason
from pubready.spec import spec_from_args


def drc_spec(**kw):
    """Build a dose_response spec as the plot entry does.

    Args:
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": "dose_response", "data": "d.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def drc_three_curves():
    """Two falling curves and one flat curve, with the same fixed noise as in R.

    Returns:
        pandas.DataFrame: The tidy table.
    """
    dose = np.repeat([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30], 3)
    noise = 2 * np.sin(np.arange(1, len(dose) + 1) * 1.7)
    lx = np.log10(dose)
    return pd.DataFrame({
        "compound": np.repeat(["active A", "active B", "inactive"], len(dose)),
        "dose": np.tile(dose, 3),
        "response": np.concatenate([fpl(lx, 5, 100, np.log10(0.3), -1) + noise,
                                    fpl(lx, 5, 100, np.log10(3), -1) + noise,
                                    100 + noise]),
    })


def test_dose_response_reports_a_flat_curve_as_not_determined():
    res = recipe_dose_response(drc_three_curves(),
                               drc_spec(x="dose", y="response", group="compound"))
    st = res["stats"]
    assert list(st["fit_status"]) == ["fitted", "fitted", "not determined"]
    assert np.isnan(st["ec50"].iloc[2])
    assert st[["ec50_ci_low", "ec50_ci_high", "hill"]].iloc[2].isna().all()
    assert st["status_reason"].iloc[2] == "the mean responses vary by less than 10% of the data range"
    assert st["ec50"].iloc[:2].to_numpy() == pytest.approx([0.3, 3], rel=0.1)
    assert res["test_meta"]["n_determined"] == 2
    assert res["test_meta"]["shared_ec50_test"]["df1"] == 1
    assert "not determined" in res["methods"]
    assert "Not determined, left out of the test: inactive" in res["fig"].axes[0].get_title()


def test_the_not_determined_rule_names_the_same_reason_as_the_r_engine():
    lx = np.array([-2.0, -1.0, 0.0, 1.0])
    assert status_reason(False, 0, 100, 0, lx, 0.1, 10, 100) == "the fit did not converge"
    assert (status_reason(True, 0, 5, 0, lx, 0.1, 10, 100)
            == "the response span is below 10% of the data range")
    assert (status_reason(True, 0, 100, 1.5, lx, 0.1, 10, 100)
            == "the EC50 lies outside the tested doses")
    assert (status_reason(True, 0, 100, 0, lx, 0, np.inf, 100)
            == "the 95% CI of the EC50 is not finite")
    assert status_reason(True, 0, 100, 0, lx, 0.1, 10, 100) is None


def test_the_check_before_the_fit_finds_a_flat_curve_as_the_r_engine_does():
    lx = np.repeat([-2.0, -1.0, 0.0, 1.0], 2)
    assert is_flat(lx, np.array([100, 102, 99, 101, 103, 101, 100, 100.0]), 100)
    assert not is_flat(lx, np.array([100, 102, 90, 92, 50, 52, 10, 12.0]), 100)
    # A single curve is its own reference, so it is never flat.
    y = np.array([100, 101, 100, 101, 101, 102, 100, 101.0])
    assert not is_flat(lx, y, 1.5)


def test_the_legend_key_of_a_fitted_curve_shows_a_point_on_a_line():
    res = recipe_dose_response(drc_three_curves(),
                               drc_spec(x="dose", y="response", group="compound"))
    legend = res["fig"].axes[0].get_legend()
    keys = legend.legend_handles
    assert [k.get_marker() for k in keys[:2]] == ["o", "o"]
    assert [k.get_linestyle() for k in keys[:2]] == ["-", "-"]
    assert legend.get_texts()[2].get_text() == "inactive (not determined)"
