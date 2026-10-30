"""The v3 fixes of the clinical recipes. The R suite (test-v3-clinical.R) runs
the same cases and expects the same results."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.cox_forest import recipe_cox_forest
from pubready.recipes.survival_km import event_codes, recipe_survival_km
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
