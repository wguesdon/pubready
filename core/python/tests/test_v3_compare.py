"""The v3 fixes of the comparison recipes.

core/r/tests/testthat/test-v3-compare.R holds the same cases for the R engine.
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.correlation import recipe_correlation
from pubready.recipes.feature_panel import recipe_feature_panel
from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.recipes.paired_compare import recipe_paired_compare
from pubready.recipes.two_group_compare import recipe_two_group_compare
from pubready.spec import spec_from_args
from pubready.util import rank_test_min_p, shapiro_normal


def v3_spec(recipe, **kw):
    """Build a figure specification for a test.

    Args:
        recipe: The recipe name.
        **kw: The command-line options.

    Returns:
        dict: The specification.
    """
    opts = {"x": None, "y": None}
    opts.update(kw)
    return spec_from_args(SimpleNamespace(recipe=recipe, data="toy.csv", **opts))


# F43: a constant sample has one status, and a test that divides by a variance
# of zero stops with a message that names the samples.
def test_a_constant_sample_is_not_normal_and_has_no_shapiro_p():
    p, normal = shapiro_normal([5.0] * 6)
    assert np.isnan(p)
    assert normal is False


def test_a_forced_t_test_on_two_constant_groups_stops_and_names_them():
    df = pd.DataFrame({"group": ["a"] * 4 + ["b"] * 4, "value": [5.0] * 4 + [6.0] * 4})
    with pytest.raises(ValueError, match="the groups 'a' and 'b' have zero variance"):
        recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="value",
                                             test="welch_t"))


def test_welch_anova_on_a_constant_group_stops_and_names_it():
    df = pd.DataFrame({"group": ["het"] * 4 + ["ko"] * 4 + ["wt"] * 4,
                       "value": [5, 5, 5, 5, 1, 2, 3, 4, 2, 4, 6, 9]})
    with pytest.raises(ValueError, match="the group 'het' has zero variance"):
        recipe_multi_group_compare(df, v3_spec("multi_group_compare", x="group", y="value",
                                               test="welch"))


def test_a_paired_t_test_on_constant_differences_stops():
    df = pd.DataFrame({"subject": list(range(1, 6)) * 2, "condition": ["pre"] * 5 + ["post"] * 5,
                       "value": list(range(1, 6)) + [v + 3 for v in range(1, 6)]})
    with pytest.raises(ValueError, match="within-subject differences have zero variance"):
        recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value",
                                          id="subject", test="t"))


def test_correlation_stops_on_a_constant_column_and_names_it():
    df = pd.DataFrame({"a": [1, 2, 3, 4, 5.0], "b": [2.0] * 5})
    with pytest.raises(ValueError, match="the column 'b' is constant"):
        recipe_correlation(df, v3_spec("correlation", x="a", y="b"))


def test_feature_panel_reports_a_constant_feature_as_not_tested():
    df = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6,
                       "f1": [1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8.5],
                       "floor": [0.5] * 12})
    res = recipe_feature_panel(df, v3_spec("feature_panel", x="group"))
    row = res["stats"].set_index("feature").loc["floor"]
    assert row["note"] == "not tested: constant"
    assert np.isnan(row["p_value"])
    assert res["test_meta"]["n_constant"] == 1
    assert res["test_meta"]["assumptions"]["normality"]["n_checked"] == 2
    assert "1 feature(s) were constant" in res["methods"]


# F61: --test auto keeps the parametric test when the rank test cannot reach
# P < 0.05 at the group sizes, and the methods text gives the attainable P.
def test_rank_test_min_p_gives_the_smallest_attainable_p_of_each_rank_test():
    assert rank_test_min_p([3, 3], "mann_whitney") == pytest.approx(0.1)
    assert rank_test_min_p([5], "signed_rank") == pytest.approx(0.0625)
    assert rank_test_min_p([3, 3, 3], "kruskal") < 0.05
    assert np.isnan(rank_test_min_p([3, 0], "mann_whitney"))


def test_auto_keeps_welch_t_when_the_rank_test_cannot_reach_005():
    df = pd.DataFrame({"condition": ["control"] * 3 + ["treated"] * 3,
                       "ct": [5, 5, 6.2, 7.1, 7.4, 7.6]})
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="condition", y="ct"))
    assert res["stats"]["test"].iloc[0] == "Welch two-sample t-test"
    assert res["test_meta"]["rank_test_attainable_p"] == pytest.approx(0.1)
    assert "smallest attainable P = 0.10), so the parametric test was kept" in res["methods"]
    forced = recipe_two_group_compare(df, v3_spec("two_group_compare", x="condition", y="ct",
                                                  test="wilcoxon"))
    assert forced["stats"]["test"].iloc[0] == "Mann-Whitney U test"
    assert "smallest attainable P = 0.10)." in forced["methods"]


def test_paired_compare_keeps_the_paired_t_test_at_5_pairs():
    pre = [1, 2, 3, 4, 5]
    post = [v + d for v, d in zip(pre, [1, 1, 1, 1.2, 5])]
    df = pd.DataFrame({"subject": list(range(1, 6)) * 2, "condition": ["pre"] * 5 + ["post"] * 5,
                       "value": pre + post})
    res = recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value",
                                            id="subject"))
    assert res["stats"]["test"].iloc[0] == "paired t-test"
    assert res["test_meta"]["rank_test_attainable_p"] == pytest.approx(0.0625)
    assert "so the parametric test was kept" in res["methods"]


def test_multi_group_compare_keeps_the_anova_at_2_values_per_group():
    df = pd.DataFrame({"group": ["a", "a", "b", "b", "c", "c"],
                       "value": [1, 2, 5, 6, 9, 10.5]})
    res = recipe_multi_group_compare(df, v3_spec("multi_group_compare", x="group", y="value"))
    assert res["test_meta"]["name"] in ("one-way ANOVA", "Welch's one-way ANOVA")
    assert res["test_meta"]["rank_test_attainable_p"] > 0.05
    assert "so the parametric test was kept" in res["methods"]


def test_feature_panel_keeps_the_parametric_family_at_3_against_3():
    df = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3,
                       "f1": [5, 5, 6.2, 7.1, 7.4, 7.6], "f2": [1, 2, 3, 2, 3, 4.5]})
    res = recipe_feature_panel(df, v3_spec("feature_panel", x="group"))
    assert res["test_meta"]["test"] == "welch_t"
    assert res["test_meta"]["kept_parametric"] is True
    assert "so the parametric test was kept" in res["methods"]


# F66: the paired path of two_group_compare names the effect Cohen's dz and
# writes no variance sentence.
def test_the_paired_t_test_of_two_group_compare_reports_cohens_dz():
    a = [1, 2, 3, 4, 5, 6]
    b = [2.1, 2.9, 4.4, 4.8, 6.3, 7.2]
    df = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6, "value": a + b})
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="value",
                                               paired=True, test="t"))
    row = res["stats"].iloc[0]
    assert row["effect_size_name"] == "Cohen's dz"
    d = np.array(a) - np.array(b)
    assert row["effect_size"] == pytest.approx(d.mean() / d.std(ddof=1), rel=1e-9)
    assert "variances" not in res["methods"]
    assert "Cohen's dz = -4.82 (a minus b)" in res["methods"]


# F36: under --paired, --test auto checks the normality of the within-pair
# differences, as paired_compare does.
def test_two_group_compare_paired_chooses_the_test_from_the_differences():
    sleep_df = pd.DataFrame({"group": ["1"] * 10 + ["2"] * 10,
                             "extra": [0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0,
                                       1.9, 0.8, 1.1, 0.1, -0.1, 4.4, 5.5, 1.6, 4.6, 3.4]})
    res = recipe_two_group_compare(sleep_df, v3_spec("two_group_compare", x="group", y="extra",
                                                     paired=True))
    row = res["stats"].iloc[0]
    assert row["test"] == "Wilcoxon signed-rank test"
    assert row["shapiro_p_diff"] == pytest.approx(0.03334161, abs=1e-6)
    assert res["qc"]["quantity"] == "within-pair differences"
    assert "normality of the within-pair differences (1 minus 2)" in res["methods"]
    with pytest.raises(ValueError, match="need the same size"):
        recipe_two_group_compare(sleep_df.iloc[1:], v3_spec("two_group_compare", x="group",
                                                           y="extra", paired=True))


# F65: Cohen's d of two_group_compare uses the pooled SD in both engines.
def test_two_group_compare_reports_cohens_d_over_the_pooled_sd():
    a = [7, 7.4, 9.5, 7.2]
    b = [9.6, 11.2, 10.1, 12.5, 8.8, 13, 11.4]
    df = pd.DataFrame({"group": ["F"] * 4 + ["M"] * 7, "hwt": a + b})
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="hwt",
                                               test="student_t"))
    sp = np.sqrt((3 * np.var(a, ddof=1) + 6 * np.var(b, ddof=1)) / 9)
    assert res["stats"]["effect_size"].iloc[0] == pytest.approx((np.mean(a) - np.mean(b)) / sp,
                                                                rel=1e-12)
    assert "Cohen's d is the mean difference over the pooled standard deviation." in res["methods"]


# F38, F64: one level rule in both engines, with --levels in each comparison
# recipe.
def test_numeric_groups_sort_as_numbers_in_the_comparison_recipes():
    df = pd.DataFrame({"conc": [160.0] * 5 + [80.0] * 5,
                       "total": [20, 22, 19, 25, 23, 30, 31, 28, 33, 29]})
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="conc", y="total",
                                               test="welch_t"))
    row = res["stats"].iloc[0]
    assert (row["group1"], row["group2"]) == ("80", "160")
    assert row["statistic"] > 0
    mg = pd.DataFrame({"dose": np.repeat([0, 160, 235, 310, 80], 4),
                       "y": [1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6]})
    res = recipe_multi_group_compare(mg, v3_spec("multi_group_compare", x="dose", y="y",
                                                 test="anova"))
    assert list(res["stats"]["group1"][:4]) == ["0", "0", "0", "0"]
    assert list(res["stats"]["group2"][:4]) == ["80", "160", "235", "310"]


def test_text_groups_sort_by_code_point_and_levels_sets_the_order():
    df = pd.DataFrame({"group": ["control"] * 5 + ["Drug"] * 5,
                       "v": [1, 2, 3, 4, 5, 3, 4, 5, 6, 7.5]})
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="v",
                                               test="wilcoxon"))
    assert res["stats"]["group1"].iloc[0] == "Drug"
    res = recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="v",
                                               test="wilcoxon", levels="control,Drug"))
    assert res["stats"]["group1"].iloc[0] == "control"


def test_paired_compare_takes_the_level_rule_not_the_row_order():
    df = pd.DataFrame({"subject": list(range(1, 7)) * 2, "condition": ["pre"] * 6 + ["post"] * 6,
                       "value": [1, 2, 3, 4, 5, 6, 2.1, 2.9, 4.4, 4.8, 6.3, 7.2]})
    res = recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value",
                                            id="subject", test="t"))
    assert res["stats"]["condition1"].iloc[0] == "post"
    res = recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value",
                                            id="subject", test="t", levels="pre,post"))
    assert res["stats"]["condition1"].iloc[0] == "pre"
    assert res["stats"]["contrast"].iloc[0] == "post - pre"


# F63: each recipe rejects a --test value that it does not know and lists the
# values that it knows.
def test_an_unknown_test_value_stops_the_comparison_recipes():
    from pubready.recipes.factorial_anova import recipe_factorial_anova
    from pubready.recipes.proportions import recipe_proportions
    pr = pd.DataFrame({"arm": ["a"] * 10 + ["b"] * 10, "response": ["yes", "no"] * 10})
    with pytest.raises(ValueError, match="proportions does not know --test 'bogus'; use one of: "
                                         "auto, chisq"):
        recipe_proportions(pr, v3_spec("proportions", x="arm", y="response", test="bogus"))
    fa = pd.DataFrame({"g": ["wt"] * 8 + ["ko"] * 8, "t": ["v", "d"] * 8, "y": range(1, 17)})
    with pytest.raises(ValueError, match="factorial_anova does not know --test 'bogus'; use one "
                                         "of: auto, anova, aov, art"):
        recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y", test="bogus"))
    with pytest.raises(ValueError, match="multi_group_compare does not know --test 'welch_t'"):
        recipe_multi_group_compare(fa, v3_spec("multi_group_compare", x="g", y="y", test="welch_t"))


# F55: a subject with two rows in one condition stops paired_compare.
def test_paired_compare_stops_on_a_subject_with_two_rows_in_one_condition():
    df = pd.DataFrame({"subject": [f"S{i:02d}" for i in range(1, 5)] * 2 + ["S02"],
                       "condition": ["pre"] * 4 + ["post"] * 4 + ["pre"],
                       "value": [1, 2, 3, 4, 2, 3, 5, 6, 250]})
    with pytest.raises(ValueError, match="'subject' repeats the subject 'S02' within a condition"):
        recipe_paired_compare(df, v3_spec("paired_compare", x="condition", y="value", id="subject"))


# F90: a constant outcome stops with a message that names it.
def test_a_constant_outcome_stops_the_comparison_recipes_and_factorial_anova():
    from pubready.recipes.factorial_anova import recipe_factorial_anova
    df = pd.DataFrame({"group": ["a"] * 6 + ["b"] * 6, "value": [50.0] * 12})
    with pytest.raises(ValueError, match="the outcome 'value' has zero variance"):
        recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="value"))
    fa = pd.DataFrame({"g": ["wt"] * 8 + ["ko"] * 8, "t": ["v", "d"] * 8, "y": [3.0] * 16})
    with pytest.raises(ValueError, match="factorial_anova cannot run because the outcome 'y' has "
                                         "zero variance"):
        recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y"))


# F83: a variance that is not finite stops with a clear message.
def test_values_near_1e300_or_1e_minus_300_stop_with_a_clear_message():
    df = pd.DataFrame({"group": ["a"] * 4 + ["b"] * 4,
                       "value": [1.1e300, 1.0e300, 1.2e300, 1.3e300, 1.5e300, 1.6e300, 1.4e300,
                                 1.7e300]})
    with pytest.raises(ValueError, match="the variance of the outcome 'value' is not a finite "
                                         "positive number"):
        recipe_two_group_compare(df, v3_spec("two_group_compare", x="group", y="value",
                                             test="welch_t"))
    cd = pd.DataFrame({"a": [-8.6e-301, -2.8e-301, -1.1e-300, 5e-301], "b": [1, 2, 3, 4.5]})
    with pytest.raises(ValueError, match="the variance of the column 'a' is not a finite positive "
                                         "number"):
        recipe_correlation(cd, v3_spec("correlation", x="a", y="b"))


# F124: the data log of factorial_anova counts each reason on its own.
def test_factorial_anova_logs_each_reason_for_a_dropped_row():
    from pubready.recipes.factorial_anova import recipe_factorial_anova
    fa = pd.DataFrame({"g": ["wt"] * 8 + ["ko"] * 8, "t": ["v", "d"] * 8,
                       "y": ["5 mg", "2", "3", "4", "5", "6", "7", "8.5", "<0.1", "3", "4", "5",
                             "6", "7", "9", "10"]})
    fa.loc[2, "t"] = None
    res = recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y"))
    assert res["clean_steps"] == ["Coerced column 'y' to numeric; 2 value(s) became NA.",
                                  "2 row(s) have no value in 'y'.",
                                  "1 row(s) have no value in 't'.",
                                  "Dropped 3 row(s) with a missing value (16 -> 13 rows)."]


# F53, F92: factorial_anova checks its design before the fit.
def test_factorial_anova_stops_on_an_empty_cell_one_level_or_no_residual_df():
    from pubready.recipes.factorial_anova import recipe_factorial_anova
    fa = pd.DataFrame({"g": ["wt"] * 4 + ["het"] * 4 + ["ko"] * 4,
                       "t": ["v", "v", "d", "d", "v", "v", "v", "v", "v", "v", "d", "d"],
                       "y": [18, 17, 31, 26, 29, 29, 23, 29, 27, 25, 41, 37]})
    with pytest.raises(ValueError, match="no row holds g = 'het', t = 'd'"):
        recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y"))
    one = pd.DataFrame({"g": ["wt"] * 4 + ["ko"] * 4, "t": ["v"] * 8, "y": range(1, 9)})
    with pytest.raises(ValueError, match="'t' has 1 \\('v'\\)"):
        recipe_factorial_anova(one, v3_spec("factorial_anova", x="g", fill="t", y="y"))
    n1 = pd.DataFrame({"g": ["wt", "wt", "ko", "ko"], "t": ["v", "d", "v", "d"], "y": [1, 3, 2, 5]})
    with pytest.raises(ValueError, match="the 4 cells hold one value each"):
        recipe_factorial_anova(n1, v3_spec("factorial_anova", x="g", fill="t", y="y"))


# F37: under --test auto the Python factorial_anova checks the residuals, and
# stops where the R engine takes the aligned rank transform.
def test_python_factorial_anova_stops_when_the_residuals_fail_shapiro_wilk():
    from pubready.recipes.factorial_anova import recipe_factorial_anova
    y = [v + s for s in (0, 3, 6, 9) for v in list(range(1, 10)) + [100]]
    fa = pd.DataFrame({"g": ["wt"] * 20 + ["ko"] * 20, "t": (["v"] * 10 + ["d"] * 10) * 2, "y": y})
    with pytest.raises(ValueError, match="rerun with --engine r, or with --test anova"):
        recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y"))
    res = recipe_factorial_anova(fa, v3_spec("factorial_anova", x="g", fill="t", y="y",
                                             test="anova"))
    assert (res["stats"]["residual_shapiro_p"] < 0.05).all()
    assert res["test_meta"]["residual_normality"]["p"] < 0.05
