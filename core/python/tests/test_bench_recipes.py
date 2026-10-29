"""Behaviour of the bench-assay recipes (estimation_plot, enzyme_kinetics and
standard_curve) on small fixed inputs. The numbers on public data are checked in
test_benchmarks.py."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.enzyme_kinetics import recipe_enzyme_kinetics
from pubready.recipes.estimation_plot import (BLOCK, block_multipliers, park_miller_uniforms,
                                              recipe_estimation_plot, wrap)
from pubready.recipes.standard_curve import recipe_standard_curve
from pubready.spec import spec_from_args


def bench_recipe_spec(recipe, **kw):
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


def fixed_unpaired():
    """The fixed two-group dataset of the R test.

    Returns:
        pandas.DataFrame: The data.
    """
    return pd.DataFrame({"group": ["ctrl"] * 7 + ["trt"] * 6,
                         "value": [10.2, 11.5, 9.8, 12.1, 10.9, 11.3, 9.4,
                                   13.4, 12.8, 14.1, 12.2, 13.9, 11.7]})


def fixed_paired():
    """The fixed paired dataset of the R test.

    Returns:
        pandas.DataFrame: The data.
    """
    return pd.DataFrame({"id": [f"s{i}" for i in range(1, 6)] * 2,
                         "cond": ["pre"] * 5 + ["post"] * 5,
                         "value": [5.1, 4.8, 6.2, 5.5, 5.9, 5.9, 5.1, 7.0, 6.1, 6.8]})


def fpl_signal(conc, bottom, top, ec50, hill):
    """The 4PL of dose_response on the concentration scale.

    Args:
        conc: The concentrations.
        bottom: The lower plateau.
        top: The upper plateau.
        ec50: The midpoint concentration.
        hill: The Hill slope.

    Returns:
        numpy.ndarray: The signals.
    """
    conc = np.asarray(conc, dtype=float)
    return bottom + (top - bottom) / (1 + 10 ** ((np.log10(ec50) - np.log10(conc)) * hill))


def test_park_miller_reaches_the_published_state():
    # Park, Miller and Stockmeyer (1993): seed 1, multiplier 48271, state 399268537.
    u, state = park_miller_uniforms(10000, 1, block_multipliers(BLOCK))
    assert state == 399268537
    assert np.all((u > 0) & (u < 1))


def test_estimation_plot_matches_the_r_bca_interval_unpaired():
    df = fixed_unpaired()
    res = recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x="group", y="value"))
    row = res["stats"].iloc[0]
    # The R engine computes these two values from the same Park-Miller resamples.
    assert row["ci_low"] == pytest.approx(1.3148939612822748, rel=1e-9)
    assert row["ci_high"] == pytest.approx(3.2523809523809533, rel=1e-9)
    t = stats.ttest_ind(df["value"][7:], df["value"][:7], equal_var=False)
    assert row["statistic"] == pytest.approx(t.statistic, rel=1e-12)
    assert row["p_value"] == pytest.approx(t.pvalue, rel=1e-12)
    assert res["test_meta"]["reference"] == "ctrl"


def test_estimation_plot_matches_the_r_bca_interval_paired():
    res = recipe_estimation_plot(fixed_paired(), bench_recipe_spec(
        "estimation_plot", x="cond", y="value", id="id", paired=True, ref="pre"))
    row = res["stats"].iloc[0]
    assert row["ci_low"] == pytest.approx(0.4, rel=1e-9)
    assert row["ci_high"] == pytest.approx(0.8199999999999996, rel=1e-9)
    assert row["test"] == "paired t-test"
    assert row["n_ref"] == 5


def test_estimation_plot_puts_the_reference_first():
    df = pd.DataFrame({"g": ["b"] * 4 + ["c"] * 4 + ["a"] * 4,
                       "v": list(range(1, 5)) + list(range(5, 9)) + list(range(3, 7))})
    res = recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x="g", y="v", ref="c"))
    assert list(res["stats"]["reference"]) == ["c", "c"]
    assert list(res["stats"]["group"]) == ["a", "b"]
    with pytest.raises(ValueError, match="is not a level"):
        recipe_estimation_plot(df, bench_recipe_spec("estimation_plot", x="g", y="v", ref="z"))


def test_estimation_plot_paired_leaves_out_an_incomplete_subject():
    df = fixed_paired().iloc[:-1]
    res = recipe_estimation_plot(df, bench_recipe_spec(
        "estimation_plot", x="cond", y="value", id="id", paired=True, ref="pre"))
    assert res["stats"]["n_group"].iloc[0] == 4
    assert any("left out 1 subject" in s for s in res["clean_steps"])


def test_wrap_breaks_after_spaces_and_underscores():
    assert wrap("Heterozygous knockout vehicle") == "Heterozygous\nknockout\nvehicle"
    assert wrap("stimulated_plus_inhibitor") == "stimulated_\nplus_\ninhibitor"


def test_enzyme_kinetics_recovers_vmax_and_km():
    s = np.repeat([0.5, 1, 2, 4, 8, 16, 32], 2)
    df = pd.DataFrame({"s": s, "v": 10 * s / (2 + s) * np.tile([1.01, 0.99], 7)})
    res = recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x="s", y="v"))
    r = res["stats"].iloc[0]
    assert r["vmax"] == pytest.approx(10, rel=0.01)
    assert r["km"] == pytest.approx(2, rel=0.02)
    tq = stats.t.ppf(0.975, len(df) - 2)
    assert r["km_ci_high"] - r["km"] == pytest.approx(tq * r["km_se"], rel=1e-10)
    assert res["test_meta"]["shared_km_test"] is None


def test_enzyme_kinetics_tests_shared_parameters():
    s = np.tile([0.5, 1, 2, 4, 8, 16, 32], 2).astype(float)
    v = np.concatenate([10 * s[:7] / (2 + s[:7]), 10 * s[7:] / (6 + s[7:])])
    v = v * np.resize([1.01, 0.99, 1.0], 14)
    df = pd.DataFrame({"s": s, "grp": ["none"] * 7 + ["inhibitor"] * 7, "v": v})
    res = recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x="s", y="v", group="grp"))
    assert list(res["stats"]["group"]) == ["none", "inhibitor"]
    assert res["test_meta"]["shared_km_test"]["p_value"] < 0.001
    assert res["test_meta"]["shared_vmax_test"]["p_value"] > 0.05
    assert res["test_meta"]["shared_km_test"]["df1"] == 1


def test_enzyme_kinetics_needs_three_concentrations():
    df = pd.DataFrame({"s": [1, 2, 1, 2], "v": [1, 2, 1.1, 2.1]})
    with pytest.raises(ValueError, match="at least 3 distinct positive substrate"):
        recipe_enzyme_kinetics(df, bench_recipe_spec("enzyme_kinetics", x="s", y="v"))


def test_standard_curve_4pl_flags_the_out_of_range_unknowns():
    conc = np.repeat([0, 1, 3, 10, 30, 100, 300], 2).astype(float)
    sig = np.concatenate([[0.05, 0.05], fpl_signal(conc[2:], 0.05, 2.5, 25, 1.1)]) * np.tile([1.01, 0.99], 7)
    unk = np.concatenate([[0.02], fpl_signal([2, 20, 200], 0.05, 2.5, 25, 1.1), [2.6]])
    df = pd.DataFrame({"id": [f"std{i}" for i in range(1, 15)] + [f"u{i}" for i in range(1, 6)],
                       "conc": np.concatenate([conc, [np.nan] * 5]), "od": np.concatenate([sig, unk])})
    res = recipe_standard_curve(df, bench_recipe_spec("standard_curve", x="conc", y="od", id="id"))
    st = res["stats"]
    assert res["test_meta"]["model"] == "four-parameter logistic"
    assert list(st["flag"]) == ["below range", "in range", "in range", "in range", "above range"]
    assert st["concentration"].iloc[[0, 4]].isna().all()
    assert st["concentration"].iloc[1:4].to_numpy() == pytest.approx([2, 20, 200], rel=0.03)
    rec = res["test_meta"]["recovery"]
    assert rec["n_levels"] == 6
    # Near the plateaus a 1% error in signal is a large error in concentration, so
    # only the levels on the steep part of the curve are held to a tight recovery.
    assert all(abs(lev["recovery_pct"] - 100) < 3 for lev in rec["levels"][2:4])
    assert any("concentration of zero out of the four-parameter" in s for s in res["clean_steps"])


def test_standard_curve_falls_back_to_a_line_and_names_rows():
    df = pd.DataFrame({"conc": [0, 0, 10, 10, 20, 20, 40, 40, np.nan, np.nan],
                       "signal": [0.02, 0.03, 0.21, 0.23, 0.41, 0.44, 0.8, 0.83, 0.33, 0.95]})
    res = recipe_standard_curve(df, bench_recipe_spec("standard_curve", x="conc", y="signal"))
    st = res["stats"]
    assert res["test_meta"]["model"] == "linear"
    assert list(st["id"]) == ["row 9", "row 10"]
    slope, intercept = np.polyfit(df["conc"][:8], df["signal"][:8], 1)
    assert st["concentration"].iloc[0] == pytest.approx((0.33 - intercept) / slope, rel=1e-12)
    assert st["flag"].iloc[1] == "above range"
    with pytest.raises(ValueError, match="needs at least 4 distinct positive"):
        recipe_standard_curve(df, bench_recipe_spec("standard_curve", x="conc", y="signal", test="4pl"))


def test_standard_curve_reads_a_falling_curve():
    conc = np.array([1, 3, 10, 30, 100, 300], dtype=float)
    od = fpl_signal(conc, 0.1, 2, 20, -1) * np.array([1.01, 0.99, 1.02, 0.98, 1.01, 0.99])
    df = pd.DataFrame({"conc": np.concatenate([conc, [np.nan, np.nan]]),
                       "od": np.concatenate([od, [2.05, 0.05]])})
    res = recipe_standard_curve(df, bench_recipe_spec("standard_curve", x="conc", y="od"))
    # A competitive assay: a high signal means a low concentration.
    assert list(res["stats"]["flag"]) == ["below range", "above range"]
