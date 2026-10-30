"""The v3 fixes of the clinical recipes. The R suite (test-v3-clinical.R) runs
the same cases and expects the same results."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.survival_km import recipe_survival_km
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
