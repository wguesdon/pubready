import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import stats

from pubplot.qc import qc_moment_stats, qc_normality
from pubplot.recipes.two_group_compare import recipe_two_group_compare
from pubplot.spec import spec_from_args


def test_qc_moment_stats_reports_moments_and_guards_small_n():
    small = qc_moment_stats([1, 2])
    assert small["n"] == 2
    assert math.isnan(small["skewness"])

    q = stats.norm.ppf((np.arange(1, 201) - 0.5) / 200)   # exact normal quantiles
    norm = qc_moment_stats(q)
    assert abs(norm["skewness"]) < 0.1                    # symmetric
    assert not math.isnan(norm["kurtosis"])
    assert norm["shapiro_p"] > 0.05


def test_qc_normality_one_panel_per_named_vector():
    rng = np.random.RandomState(0)
    qc = {"quantity": "group values",
          "panels": {"ctrl": rng.normal(size=30), "treat": rng.normal(size=30)}}
    fig, panels = qc_normality(qc)
    assert len(panels) == 2
    assert {p["panel"] for p in panels} == {"ctrl", "treat"}
    assert all("shapiro_p" in p and "skewness" in p and "kurtosis" in p for p in panels)


def test_two_group_exposes_checked_group_values_for_qc():
    df = pd.DataFrame({"group": ["a"] * 10 + ["b"] * 10,
                       "value": list(range(10)) + list(range(3, 13))})
    a = SimpleNamespace(recipe="two_group_compare", data="toy.csv",
                        x="group", y="value", test="auto")
    res = recipe_two_group_compare(df, spec_from_args(a))
    assert res["qc"]["quantity"] == "group values"
    assert set(res["qc"]["panels"].keys()) == {"a", "b"}
