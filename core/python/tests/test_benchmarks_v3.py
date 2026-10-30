"""Benchmarks of the v3 numbers audit: 32 cases on public datasets, each checked
against a reference value that an independent call to the canonical R function
computed (benchmarks/reference_values.csv, sections 40 to 62 of
tools/make_benchmarks/make_benchmark_data.R). ``figkit test`` mounts the
benchmarks folder; the tests skip when it is absent. The checks mirror
core/r/tests/testthat/test-benchmarks-v3.R quantity for quantity."""
import itertools
import os
import re
from types import SimpleNamespace

import pandas as pd
import pytest

from pubready.io import read_tidy
from pubready.recipes.bland_altman import recipe_bland_altman
from pubready.recipes.composition import recipe_composition
from pubready.recipes.correlation import recipe_correlation
from pubready.recipes.cox_forest import recipe_cox_forest
from pubready.recipes.dose_response import recipe_dose_response
from pubready.recipes.enzyme_kinetics import recipe_enzyme_kinetics
from pubready.recipes.estimation_plot import recipe_estimation_plot
from pubready.recipes.factorial_anova import recipe_factorial_anova
from pubready.recipes.feature_panel import recipe_feature_panel
from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.recipes.paired_compare import recipe_paired_compare
from pubready.recipes.pca import recipe_pca
from pubready.recipes.roc_curve import recipe_roc_curve
from pubready.recipes.standard_curve import recipe_standard_curve
from pubready.recipes.survival_km import recipe_survival_km
from pubready.recipes.timecourse import recipe_timecourse
from pubready.recipes.two_group_compare import recipe_two_group_compare
from pubready.spec import spec_from_args

BENCH = os.environ.get("PUBREADY_BENCHMARKS", "")
pytestmark = pytest.mark.skipif(not BENCH or not os.path.isdir(BENCH),
                                reason="benchmarks not mounted")

ROUTE_TO_R = "--engine r"
"""The routing text of a Python stop for a method that only the R engine has."""


def bench_data(name):
    """Read one benchmark dataset.

    Args:
        name: The dataset name, without the extension.

    Returns:
        pandas.DataFrame: The data, read as the pipeline reads it.
    """
    return read_tidy(os.path.join(BENCH, "data", f"{name}.csv"))


def bench_spec(recipe, **kw):
    """Build a figure spec as the plot entry does.

    Args:
        recipe: The recipe name.
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": recipe, "data": "bench.csv", "x": None, "y": None, "test": "auto"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def bench_refs():
    """Read the reference values.

    Returns:
        pandas.DataFrame: One row per case and quantity.
    """
    return pd.read_csv(os.path.join(BENCH, "reference_values.csv"))


def expect_ref(actual, case, quantity):
    """Assert that a value matches its reference within the stated tolerance.

    Args:
        actual: The value the recipe computed.
        case: The benchmark case.
        quantity: The quantity name in reference_values.csv.
    """
    refs = bench_refs()
    r = refs[(refs["case"] == case) & (refs["quantity"] == quantity)]
    assert len(r) == 1, f"{case} {quantity}: no single reference"
    expected = float(r["value"].iloc[0])
    assert float(actual) == pytest.approx(expected, rel=float(r["rel_tol"].iloc[0])), \
        f"{case} {quantity}: {actual} against {expected}"


def pair_row(stats, a, b, c1="group1", c2="group2"):
    """Select the row of one unordered pair.

    Args:
        stats: The stats table.
        a: One member of the pair.
        b: The other member.
        c1: The first pair column.
        c2: The second pair column.

    Returns:
        pandas.Series: The row.
    """
    sel = ((stats[c1] == a) & (stats[c2] == b)) | ((stats[c1] == b) & (stats[c2] == a))
    return stats[sel].iloc[0]


def dropped_rows(steps, reason):
    """The count that a "Dropped <n> row(s) with <reason>" cleaning step names.

    Args:
        steps: The cleaning steps of a recipe result.
        reason: The text that identifies the step.

    Returns:
        float: The count, or NaN when no single step holds the reason.
    """
    hit = [s for s in steps if reason in s]
    if len(hit) != 1:
        return float("nan")
    return float(re.match(r"^Dropped (\d+) row", hit[0]).group(1))


# --- two_group_compare -------------------------------------------------------

def test_cats_student_t_and_pooled_cohens_d():
    res = recipe_two_group_compare(bench_data("cats"), bench_spec(
        "two_group_compare", x="Sex", y="Hwt", test="student_t"))
    assert res["test_meta"]["name"] == "Student's two-sample t-test"
    s = res["stats"].iloc[0]
    assert s["group1"] == "F"
    expect_ref(s["statistic"], "cats_student", "statistic")
    expect_ref(s["df"], "cats_student", "df")
    expect_ref(s["p_value"], "cats_student", "p_value")
    expect_ref(s["effect_size"], "cats_student", "cohens_d")


def test_cats_mann_whitney_with_continuity_correction():
    res = recipe_two_group_compare(bench_data("cats"), bench_spec(
        "two_group_compare", x="Sex", y="Hwt", test="wilcoxon"))
    assert res["test_meta"]["name"] == "Mann-Whitney U test"
    s = res["stats"].iloc[0]
    expect_ref(s["statistic"], "cats_mann_whitney", "statistic")
    expect_ref(s["p_value"], "cats_mann_whitney", "p_value")
    expect_ref(s["effect_size"], "cats_mann_whitney", "rank_biserial")


def test_plant_growth_exact_mann_whitney():
    res = recipe_two_group_compare(bench_data("plant_growth_ctrl_trt2"), bench_spec(
        "two_group_compare", x="group", y="weight", test="wilcoxon"))
    s = res["stats"].iloc[0]
    expect_ref(s["statistic"], "plant_growth_mw_exact", "statistic")
    expect_ref(s["p_value"], "plant_growth_mw_exact", "p_value")
    expect_ref(s["effect_size"], "plant_growth_mw_exact", "rank_biserial")


def test_sleep_two_group_paired_t_and_dz():
    res = recipe_two_group_compare(bench_data("sleep"), bench_spec(
        "two_group_compare", x="group", y="extra", paired=True, test="t"))
    assert res["test_meta"]["name"] == "paired t-test"
    s = res["stats"].iloc[0]
    assert s["effect_size_name"] == "Cohen's dz"
    expect_ref(s["statistic"], "sleep_two_group_paired", "statistic")
    expect_ref(s["df"], "sleep_two_group_paired", "df")
    expect_ref(s["p_value"], "sleep_two_group_paired", "p_value")
    expect_ref(s["effect_size"], "sleep_two_group_paired", "cohens_dz")


def test_sleep_two_group_paired_auto_picks_signed_rank():
    res = recipe_two_group_compare(bench_data("sleep"), bench_spec(
        "two_group_compare", x="group", y="extra", paired=True))
    s = res["stats"].iloc[0]
    expect_ref(s["shapiro_p_diff"], "sleep_two_group_paired_auto", "shapiro_p_diff")
    expect_ref(res["test_meta"]["assumptions"]["normality_of_differences"]["p"],
               "sleep_two_group_paired_auto", "shapiro_p_diff")
    assert res["test_meta"]["name"] == "Wilcoxon signed-rank test"
    expect_ref(s["statistic"], "sleep_two_group_paired_auto", "statistic")
    expect_ref(s["p_value"], "sleep_two_group_paired_auto", "p_value")


def test_nitrofen_two_group_sorts_doses_as_numbers():
    res = recipe_two_group_compare(bench_data("nitrofen_80_160"), bench_spec(
        "two_group_compare", x="conc", y="total", test="welch_t"))
    s = res["stats"].iloc[0]
    assert [s["group1"], s["group2"]] == ["80", "160"]
    expect_ref(s["statistic"], "nitrofen_two_group", "statistic")
    expect_ref(s["p_value"], "nitrofen_two_group", "p_value")


def test_cars93_origin_usa_is_group_one():
    res = recipe_two_group_compare(bench_data("cars93"), bench_spec(
        "two_group_compare", x="Origin", y="MPG.highway", test="wilcoxon"))
    s = res["stats"].iloc[0]
    assert [s["group1"], s["group2"]] == ["USA", "non-USA"]
    expect_ref(s["statistic"], "cars93_origin", "statistic")
    expect_ref(s["p_value"], "cars93_origin", "p_value")
    expect_ref(s["effect_size"], "cars93_origin", "rank_biserial")


# --- feature_panel -----------------------------------------------------------

def test_cats_feature_panel_welch_pooled_d_and_bh():
    res = recipe_feature_panel(bench_data("cats"), bench_spec(
        "feature_panel", x="Sex", test="welch_t"))
    st = res["stats"].set_index("feature")
    for f in ("Bwt", "Hwt"):
        expect_ref(st.loc[f, "statistic"], "cats_feature_panel", f"statistic_{f}")
        expect_ref(st.loc[f, "df1"], "cats_feature_panel", f"df_{f}")
        expect_ref(st.loc[f, "p_value"], "cats_feature_panel", f"p_{f}")
        expect_ref(st.loc[f, "effect_size"], "cats_feature_panel", f"cohens_d_{f}")
        expect_ref(st.loc[f, "p_adj"], "cats_feature_panel", f"p_adj_{f}")


# --- multi_group_compare -----------------------------------------------------

def test_plant_growth_every_tukey_pair_and_eta_squared():
    res = recipe_multi_group_compare(bench_data("plant_growth"), bench_spec(
        "multi_group_compare", x="group", y="weight", test="anova"))
    for a, b in (("ctrl", "trt1"), ("ctrl", "trt2"), ("trt1", "trt2")):
        r = pair_row(res["stats"], a, b)
        expect_ref(r["estimate"], "plant_growth_tukey", f"estimate_{b}_{a}")
        expect_ref(r["p_adj"], "plant_growth_tukey", f"p_{b}_{a}")
    expect_ref(res["test_meta"]["omnibus"]["effect_size"]["value"], "plant_growth_tukey",
               "eta_squared")


def test_chickwts_every_games_howell_pair():
    res = recipe_multi_group_compare(bench_data("chickwts"), bench_spec(
        "multi_group_compare", x="feed", y="weight", test="welch"))
    assert res["test_meta"]["post_hoc"] == "the Games-Howell test"
    feeds = ("casein", "horsebean", "linseed", "meatmeal", "soybean", "sunflower")
    for a, b in itertools.combinations(feeds, 2):
        r = pair_row(res["stats"], a, b)
        expect_ref(r["estimate"], "chickwts_games_howell", f"estimate_{b}_{a}")
        expect_ref(r["p_adj"], "chickwts_games_howell", f"p_{b}_{a}")


INSECT_PAIRS = tuple(itertools.combinations("ABCDEF", 2))


def test_insect_sprays_dunn_bh_and_eta_squared_h():
    res = recipe_multi_group_compare(bench_data("insect_sprays"), bench_spec(
        "multi_group_compare", x="spray", y="count", test="kruskal"))
    for a, b in INSECT_PAIRS:
        expect_ref(pair_row(res["stats"], a, b)["p_adj"], "insect_sprays_dunn", f"p_adj_{a}_{b}")
    expect_ref(res["test_meta"]["omnibus"]["effect_size"]["value"], "insect_sprays_dunn",
               "eta_squared_h")


@pytest.mark.xfail(strict=True, reason="insect_sprays_dunn: the Dunn table has no z column (N-11)")
def test_insect_sprays_dunn_z():
    res = recipe_multi_group_compare(bench_data("insect_sprays"), bench_spec(
        "multi_group_compare", x="spray", y="count", test="kruskal"))
    st = res["stats"]
    for a, b in INSECT_PAIRS:
        r = st[(st["group1"] == a) & (st["group2"] == b)].iloc[0]
        expect_ref(r["statistic"], "insect_sprays_dunn", f"z_{a}_{b}")


def test_nitrofen_multi_sorts_doses_and_matches_anova_and_tukey():
    res = recipe_multi_group_compare(bench_data("nitrofen"), bench_spec(
        "multi_group_compare", x="conc", y="total", test="anova"))
    assert list(res["stats"]["group1"].unique()) == ["0", "80", "160", "235"]
    om = res["test_meta"]["omnibus"]
    expect_ref(om["statistic"], "nitrofen_multi", "statistic")
    expect_ref(om["p_value"], "nitrofen_multi", "p_value")
    r = pair_row(res["stats"], "80", "160")
    expect_ref(r["estimate"], "nitrofen_multi", "estimate_160_80")
    expect_ref(r["p_adj"], "nitrofen_multi", "p_160_80")
    r = pair_row(res["stats"], "0", "80")
    expect_ref(r["estimate"], "nitrofen_multi", "estimate_80_0")
    expect_ref(r["p_adj"], "nitrofen_multi", "p_80_0")


# --- factorial_anova ---------------------------------------------------------

QUINE_TERMS = ("Eth", "Sex", "Lrn", "Eth:Sex", "Eth:Lrn", "Sex:Lrn", "Eth:Sex:Lrn")


def quine_spec(test):
    """The spec of the quine three-way design.

    Args:
        test: The --test value.

    Returns:
        dict: The spec.
    """
    return bench_spec("factorial_anova", x="Eth", fill="Sex", facet="Lrn", y="Days", test=test)


def test_quine_threeway_type_two_anova():
    res = recipe_factorial_anova(bench_data("quine"), quine_spec("anova"))
    st = res["stats"].set_index("effect")
    for term in QUINE_TERMS:
        expect_ref(st.loc[term, "statistic"], "quine_threeway_anova", f"F_{term}")
        expect_ref(st.loc[term, "p_value"], "quine_threeway_anova", f"p_{term}")
        expect_ref(st.loc[term, "effect_size"], "quine_threeway_anova", f"pes_{term}")


def test_quine_threeway_art_routes_to_r():
    """quine_threeway_art is R only: the Python engine stops and routes to R."""
    with pytest.raises(ValueError, match=ROUTE_TO_R):
        recipe_factorial_anova(bench_data("quine"), quine_spec("art"))


def test_quine_threeway_auto_residual_shapiro_and_route_to_r():
    """The residual P matches; under auto the ART path stops and routes to R."""
    res = recipe_factorial_anova(bench_data("quine"), quine_spec("anova"))
    expect_ref(res["test_meta"]["residual_normality"]["p"], "quine_threeway_auto",
               "residual_shapiro_p")
    expect_ref(res["stats"]["residual_shapiro_p"].iloc[0], "quine_threeway_auto",
               "residual_shapiro_p")
    with pytest.raises(ValueError, match=ROUTE_TO_R):
        recipe_factorial_anova(bench_data("quine"), quine_spec("auto"))


def test_higgins_art_routes_to_r():
    """higgins_art is R only: the Python engine stops and routes to R."""
    with pytest.raises(ValueError, match=ROUTE_TO_R):
        recipe_factorial_anova(bench_data("higgins1990_table1"), bench_spec(
            "factorial_anova", x="Row", fill="Column", y="Response", test="art"))


# --- estimation_plot and paired_compare --------------------------------------

def test_shoes_paired_estimation_and_hedges_g():
    res = recipe_estimation_plot(bench_data("shoes"), bench_spec(
        "estimation_plot", x="material", y="wear", paired=True, id="boy", ref="A"))
    st = res["stats"]
    assert (st["group"] == "B").sum() == 1
    r = st[st["group"] == "B"].iloc[0]
    expect_ref(r["mean_diff"], "shoes_paired_estimation", "mean_diff")
    expect_ref(r["statistic"], "shoes_paired_estimation", "statistic")
    expect_ref(r["df"], "shoes_paired_estimation", "df")
    expect_ref(r["p_value"], "shoes_paired_estimation", "p_value")
    expect_ref(r["hedges_g"], "shoes_paired_estimation", "hedges_g")


def test_sleep_paired_compare_effect_is_condition_two_minus_one():
    res = recipe_paired_compare(bench_data("sleep"), bench_spec(
        "paired_compare", x="group", y="extra", id="ID", test="t"))
    s = res["stats"].iloc[0]
    expect_ref(s["effect_size"], "sleep_paired_effect", "cohens_dz")
    expect_ref(s["mean_diff"], "sleep_paired_effect", "mean_diff")


# --- timecourse --------------------------------------------------------------

def test_oxboys_one_way_repeated_measures():
    res = recipe_timecourse(bench_data("oxboys"), bench_spec(
        "timecourse", x="Occasion", y="height", id="Subject"))
    st = res["stats"]
    assert (st["effect"] == "Occasion").sum() == 1
    r = st[st["effect"] == "Occasion"].iloc[0]
    expect_ref(r["F"], "oxboys_rm", "F_time")
    expect_ref(r["p_value"], "oxboys_rm", "p_time_uncorrected")
    expect_ref(r["partial_eta2"], "oxboys_rm", "pes_time")
    expect_ref(r["gg_epsilon"], "oxboys_rm", "gg_epsilon")
    expect_ref(r["p_gg"], "oxboys_rm", "p_time_gg")
    expect_ref(res["test_meta"]["sphericity"]["p_value"], "oxboys_rm", "mauchly_p")


# --- survival ----------------------------------------------------------------

def veteran_trt_km():
    """Run survival_km on veteran by treatment arm.

    Returns:
        dict: The recipe result.
    """
    return recipe_survival_km(bench_data("veteran"), bench_spec(
        "survival_km", time="time", event="status", x="trt"))


def test_veteran_trt_medians_and_log_rank():
    res = veteran_trt_km()
    st = res["stats"].set_index("group")
    for arm in ("standard", "test"):
        expect_ref(st.loc[arm, "median_survival"], "veteran_trt_km", f"median_{arm}")
    expect_ref(res["test_meta"]["log_rank"]["chisq"], "veteran_trt_km", "log_rank_chisq")
    expect_ref(res["test_meta"]["log_rank"]["p_value"], "veteran_trt_km", "log_rank_p")


@pytest.mark.xfail(strict=True, reason="veteran_trt_km: the Python survival_km table has no CI "
                                       "of the median (an Engine limits entry)")
def test_veteran_trt_median_confidence_intervals():
    st = veteran_trt_km()["stats"].set_index("group")
    for arm in ("standard", "test"):
        expect_ref(st.loc[arm, "ci_low"], "veteran_trt_km", f"median_ci_low_{arm}")
        expect_ref(st.loc[arm, "ci_high"], "veteran_trt_km", f"median_ci_high_{arm}")


def test_colon_cox_three_level_rx():
    res = recipe_cox_forest(bench_data("colon_death"), bench_spec(
        "cox_forest", time="time", event="status", covariates="rx,sex,age"))
    st = res["stats"].set_index("term")
    keys = {"rxLev+5FU": "rxLev5FU", "rxObs": "rxObs", "sex": "sex", "age": "age"}
    for term, key in keys.items():
        expect_ref(st.loc[term, "HR"], "colon_cox", f"HR_{key}")
        expect_ref(st.loc[term, "ci_low"], "colon_cox", f"ci_low_{key}")
        expect_ref(st.loc[term, "ci_high"], "colon_cox", f"ci_high_{key}")
        expect_ref(st.loc[term, "p_value"], "colon_cox", f"p_{key}")
    expect_ref(res["test_meta"]["proportional_hazards"]["p_value"], "colon_cox", "ph_global_p")


def test_cgd_cox_x_linked_is_the_reference():
    res = recipe_cox_forest(bench_data("cgd_first"), bench_spec(
        "cox_forest", time="time", event="status", covariates="treat,inherit,age"))
    st = res["stats"].set_index("term")
    keys = {"treatrIFN-g": "treatrIFNg", "inheritautosomal": "inheritautosomal", "age": "age"}
    for term, key in keys.items():
        expect_ref(st.loc[term, "HR"], "cgd_cox", f"HR_{key}")
        expect_ref(st.loc[term, "p_value"], "cgd_cox", f"p_{key}")
    expect_ref(res["test_meta"]["proportional_hazards"]["p_value"], "cgd_cox", "ph_global_p")


# --- diagnostic and agreement ------------------------------------------------

def test_birthwt_roc_auc_delong_ci_and_youden():
    res = recipe_roc_curve(bench_data("birthwt"), bench_spec(
        "roc_curve", y="low", x="lwt", positive="1"))
    r = res["stats"].set_index("marker").loc["lwt"]
    assert r["direction"] == ">"
    expect_ref(r["auc"], "birthwt_roc", "auc")
    expect_ref(r["auc_ci_low"], "birthwt_roc", "auc_ci_low")
    expect_ref(r["auc_ci_high"], "birthwt_roc", "auc_ci_high")
    expect_ref(r["youden_threshold"], "birthwt_roc", "youden_threshold")
    expect_ref(r["sensitivity"], "birthwt_roc", "youden_sensitivity")
    expect_ref(r["specificity"], "birthwt_roc", "youden_specificity")


def test_davis_bland_altman_limits_and_proportional_bias():
    res = recipe_bland_altman(bench_data("davis"), bench_spec(
        "bland_altman", x="weight", y="repwt"))
    s = res["stats"].iloc[0]
    for q in ("n", "bias", "bias_ci_low", "bias_ci_high", "sd_diff", "loa_low", "loa_high",
              "loa_low_ci_low", "loa_high_ci_high", "prop_bias_slope", "prop_bias_p"):
        expect_ref(s[q], "davis_bland_altman", q)


# --- curve fits --------------------------------------------------------------

def test_nitrendipene_falling_curves_and_shared_ic50():
    res = recipe_dose_response(bench_data("nitrendipene"), bench_spec(
        "dose_response", x="NIF", y="activity", group="Tissue"))
    expect_ref(dropped_rows(res["clean_steps"], "with a dose of zero or below"),
               "nitrendipene_ic50", "n_dropped_zero_dose")
    cols = {"bottom": "bottom", "top": "top", "ic50": "midpoint", "hill": "hill",
            "ic50_ci_low": "midpoint_ci_low", "ic50_ci_high": "midpoint_ci_high", "rss": "rss"}
    st = res["stats"].set_index("group")
    for t in ("tissue1", "tissue2", "tissue3", "tissue4"):
        for q, col in cols.items():
            expect_ref(st.loc[t, col], "nitrendipene_ic50", f"{q}_{t}")
    sh = res["test_meta"]["shared_ec50_test"]
    expect_ref(sh["statistic"], "nitrendipene_ic50", "shared_F")
    expect_ref(sh["p_value"], "nitrendipene_ic50", "shared_p")


def test_dnase_ec50_standard_error_and_ci():
    res = recipe_dose_response(bench_data("dnase_run1"), bench_spec(
        "dose_response", x="conc", y="density"))
    s = res["stats"].iloc[0]
    expect_ref(s["se_log_midpoint"], "dnase_ec50_ci", "se_log_ec50")
    expect_ref(s["midpoint_ci_low"], "dnase_ec50_ci", "ec50_ci_low")
    expect_ref(s["midpoint_ci_high"], "dnase_ec50_ci", "ec50_ci_high")


def test_formaldehyde_linear_standard_curve():
    res = recipe_standard_curve(bench_data("formaldehyde_standard_curve"), bench_spec(
        "standard_curve", x="carb", y="optden", id="sample", test="linear"))
    tm = res["test_meta"]
    expect_ref(tm["parameters"]["intercept"], "formaldehyde_linear", "intercept")
    expect_ref(tm["parameters"]["slope"], "formaldehyde_linear", "slope")
    expect_ref(tm["r_squared"], "formaldehyde_linear", "r_squared")
    expect_ref(tm["n_in_range"], "formaldehyde_linear", "n_in_range")
    expect_ref(tm["n_below"], "formaldehyde_linear", "n_below")
    refs = bench_refs()
    conc_q = refs.loc[(refs["case"] == "formaldehyde_linear")
                      & refs["quantity"].str.startswith("conc_"), "quantity"]
    assert len(conc_q) == 5
    st = res["stats"].set_index("id")
    for q in conc_q:
        expect_ref(st.loc[q[len("conc_"):], "concentration"], "formaldehyde_linear", q)


def test_puromycin_vmax_and_km_confidence_intervals():
    res = recipe_enzyme_kinetics(bench_data("puromycin"), bench_spec(
        "enzyme_kinetics", x="conc", y="rate", group="state"))
    st = res["stats"].set_index("group")
    for s in ("treated", "untreated"):
        for q in ("vmax_ci_low", "vmax_ci_high", "km_ci_low", "km_ci_high"):
            expect_ref(st.loc[s, q], "puromycin_ci", f"{q}_{s}")


# --- composition, PCA and correlation ----------------------------------------

def test_housing_composition_kruskal_by_influence():
    res = recipe_composition(bench_data("housing"), bench_spec(
        "composition", id="sample", x="Sat", y="Freq", group="Infl"))
    st = res["stats"].drop_duplicates("category").set_index("category")
    for s in ("Low", "Medium", "High"):
        expect_ref(st.loc[s, "statistic"], "housing_kruskal", f"H_{s}")
        expect_ref(st.loc[s, "p_value"], "housing_kruskal", f"p_{s}")
        expect_ref(st.loc[s, "p_adj"], "housing_kruskal", f"p_adj_{s}")


def test_airway_and_iris_permanova_pseudo_f():
    ap = recipe_pca(bench_data("airway_pca"), bench_spec("pca", group="dex"))
    perm = ap["test_meta"]["permanova"]
    assert perm["method"] == "exact"
    expect_ref(perm["statistic_F"], "airway_pca_permanova", "F")
    expect_ref(perm["p_value"], "airway_pca_permanova", "p_exact")
    ir = recipe_pca(bench_data("iris"), bench_spec("pca", group="Species"))
    expect_ref(ir["test_meta"]["permanova"]["statistic_F"], "iris_pca_permanova", "F")


def test_mtcars_pearson_confidence_interval():
    res = recipe_correlation(bench_data("mtcars"), bench_spec(
        "correlation", x="wt", y="mpg", test="pearson"))
    s = res["stats"].iloc[0]
    expect_ref(s["conf_low"], "mtcars_pearson_ci", "conf_low")
    expect_ref(s["conf_high"], "mtcars_pearson_ci", "conf_high")
