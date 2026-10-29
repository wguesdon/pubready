"""Benchmarks: every recipe on a public dataset, checked against a reference value
that an independent call to the canonical R function computed
(benchmarks/reference_values.csv, written by tools/make_benchmarks). ``figkit
test`` mounts the benchmarks folder; the tests skip when it is absent."""
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.io import read_tidy
from pubready.recipes.bland_altman import recipe_bland_altman
from pubready.recipes.correlation import recipe_correlation
from pubready.recipes.correlation_heatmap import recipe_correlation_heatmap
from pubready.recipes.cox_forest import recipe_cox_forest
from pubready.recipes.dose_response import recipe_dose_response
from pubready.recipes.embedding import recipe_embedding
from pubready.recipes.enrichment_dot import recipe_enrichment_dot
from pubready.recipes.factorial_anova import recipe_factorial_anova
from pubready.recipes.feature_panel import recipe_feature_panel
from pubready.recipes.forest_meta import recipe_forest_meta
from pubready.recipes.heatmap import recipe_heatmap
from pubready.recipes.ma_plot import recipe_ma_plot
from pubready.recipes.manhattan import recipe_manhattan
from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.recipes.paired_compare import recipe_paired_compare
from pubready.recipes.pca import recipe_pca
from pubready.recipes.proportions import recipe_proportions
from pubready.recipes.roc_curve import recipe_roc_curve
from pubready.recipes.survival_km import recipe_survival_km
from pubready.recipes.swimmer import recipe_swimmer
from pubready.recipes.timecourse import recipe_timecourse
from pubready.recipes.two_group_compare import recipe_two_group_compare
from pubready.recipes.upset import recipe_upset
from pubready.recipes.volcano import recipe_volcano
from pubready.spec import spec_from_args

BENCH = os.environ.get("PUBREADY_BENCHMARKS", "")
pytestmark = pytest.mark.skipif(not BENCH or not os.path.isdir(BENCH),
                                reason="benchmarks not mounted")


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


def expect_ref(actual, case, quantity):
    """Assert that a value matches its reference within the stated tolerance.

    Args:
        actual: The value the recipe computed.
        case: The benchmark case.
        quantity: The quantity name in reference_values.csv.
    """
    refs = pd.read_csv(os.path.join(BENCH, "reference_values.csv"))
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


def test_sleep_paired_t_matches_t_test():
    res = recipe_paired_compare(bench_data("sleep"), bench_spec(
        "paired_compare", x="group", y="extra", id="ID", test="t"))
    s = res["stats"].iloc[0]
    expect_ref(s["statistic"], "sleep_paired", "statistic")
    expect_ref(s["df"], "sleep_paired", "df")
    expect_ref(s["p_value"], "sleep_paired", "p_value")


def test_sleep_welch_matches_t_test_and_keeps_the_sign_of_d():
    res = recipe_two_group_compare(bench_data("sleep"), bench_spec(
        "two_group_compare", x="group", y="extra", test="welch_t"))
    s = res["stats"].iloc[0]
    expect_ref(s["statistic"], "sleep_welch", "statistic")
    expect_ref(s["df"], "sleep_welch", "df")
    expect_ref(s["p_value"], "sleep_welch", "p_value")
    assert np.sign(s["effect_size"]) == np.sign(s["statistic"])


def test_plant_growth_anova_and_tukey():
    res = recipe_multi_group_compare(bench_data("plant_growth"), bench_spec(
        "multi_group_compare", x="group", y="weight"))
    assert res["test_meta"]["name"] == "one-way ANOVA"
    expect_ref(res["test_meta"]["omnibus"]["statistic"], "plant_growth_anova", "statistic")
    expect_ref(res["test_meta"]["omnibus"]["p_value"], "plant_growth_anova", "p_value")
    expect_ref(pair_row(res["stats"], "trt1", "trt2")["p_adj"], "plant_growth_anova",
               "tukey_p_trt1_trt2")


def test_insect_sprays_kruskal():
    res = recipe_multi_group_compare(bench_data("insect_sprays"), bench_spec(
        "multi_group_compare", x="spray", y="count", test="kruskal"))
    expect_ref(res["test_meta"]["omnibus"]["statistic"], "insect_sprays_kruskal", "statistic")
    expect_ref(res["test_meta"]["omnibus"]["p_value"], "insect_sprays_kruskal", "p_value")


def test_tooth_growth_two_way_anova():
    res = recipe_factorial_anova(bench_data("tooth_growth"), bench_spec(
        "factorial_anova", x="dose", fill="supp", y="len", test="anova"))
    for term in ("dose", "supp", "dose:supp"):
        r = res["stats"].set_index("effect").loc[term]
        expect_ref(r["statistic"], "tooth_growth_twoway", f"F_{term}")
        expect_ref(r["p_value"], "tooth_growth_twoway", f"p_{term}")
    assert set(res["stats"]["effect_size_name"]) == {"partial eta squared"}


def test_lung_log_rank():
    res = recipe_survival_km(bench_data("lung"), bench_spec(
        "survival_km", time="time", event="event", x="sex"))
    expect_ref(res["test_meta"]["log_rank"]["chisq"], "lung_km", "log_rank_chisq")
    expect_ref(res["test_meta"]["log_rank"]["p_value"], "lung_km", "log_rank_p")


def test_lung_cox_hazard_ratios():
    res = recipe_cox_forest(bench_data("lung"), bench_spec(
        "cox_forest", time="time", event="event", covariates="age,sex,ph.ecog"))
    st = res["stats"].set_index("term")
    for term in ("age", "sexmale", "ph.ecog"):
        expect_ref(st.loc[term, "HR"], "lung_cox", f"HR_{term}")
        expect_ref(st.loc[term, "p_value"], "lung_cox", f"p_{term}")


def test_mtcars_correlation_and_matrix():
    res = recipe_correlation(bench_data("mtcars"), bench_spec(
        "correlation", x="wt", y="mpg", test="pearson"))
    expect_ref(res["stats"]["estimate"].iloc[0], "mtcars_correlation", "estimate")
    expect_ref(res["stats"]["p_value"].iloc[0], "mtcars_correlation", "p_value")
    hm = recipe_correlation_heatmap(bench_data("mtcars"), bench_spec("correlation_heatmap"))
    expect_ref(pair_row(hm["stats"], "mpg", "wt", "var1", "var2")["r"],
               "mtcars_correlation_heatmap", "r_mpg_wt")
    expect_ref(pair_row(hm["stats"], "disp", "cyl", "var1", "var2")["r"],
               "mtcars_correlation_heatmap", "r_disp_cyl")


def test_mtcars_heatmap_welch_tests():
    spec = bench_spec("heatmap", annotation=os.path.join(BENCH, "data", "mtcars_annotation.csv"))
    res = recipe_heatmap(bench_data("mtcars_matrix"), spec)
    st = res["stats"].set_index("feature")
    for f in ("mpg", "disp", "hp", "drat", "wt", "qsec"):
        expect_ref(st.loc[f, "p_value"], "mtcars_heatmap", f"p_{f}")


def test_airquality_spearman_drops_missing_values():
    res = recipe_correlation(bench_data("airquality"), bench_spec(
        "correlation", x="Temp", y="Ozone", test="spearman"))
    expect_ref(res["stats"]["estimate"].iloc[0], "airquality_spearman", "estimate")
    expect_ref(res["stats"]["n"].iloc[0], "airquality_spearman", "n")


def test_iris_pca_variance():
    res = recipe_pca(bench_data("iris"), bench_spec("pca", group="Species"))
    expect_ref(res["test_meta"]["pc1_variance_pct"], "iris_pca", "pc1_variance_pct")
    expect_ref(res["test_meta"]["pc2_variance_pct"], "iris_pca", "pc2_variance_pct")


def test_titanic_and_ucb_chi_squared():
    ti = recipe_proportions(bench_data("titanic"), bench_spec(
        "proportions", x="Class", y="Survived", test="chisq"))
    expect_ref(ti["test_meta"]["statistic"], "titanic_proportions", "statistic")
    expect_ref(ti["test_meta"]["p_value"], "titanic_proportions", "p_value")
    ucb = recipe_proportions(bench_data("ucb_admissions"), bench_spec(
        "proportions", x="Gender", y="Admit", test="chisq"))
    expect_ref(ucb["test_meta"]["statistic"], "ucb_proportions", "statistic")
    expect_ref(ucb["test_meta"]["p_value"], "ucb_proportions", "p_value")
    expect_ref(ucb["test_meta"]["odds_ratio"], "ucb_proportions", "odds_ratio")


def test_dnase_four_parameter_logistic():
    res = recipe_dose_response(bench_data("dnase_run1"), bench_spec(
        "dose_response", x="conc", y="density"))
    for q in ("bottom", "top", "ec50", "hill", "rss"):
        expect_ref(res["stats"][q].iloc[0], "dnase_4pl", q)
    sh = recipe_dose_response(bench_data("dnase"), bench_spec(
        "dose_response", x="conc", y="density", group="Run"))
    expect_ref(sh["test_meta"]["shared_ec50_test"]["statistic"], "dnase_shared", "statistic")
    expect_ref(sh["test_meta"]["shared_ec50_test"]["p_value"], "dnase_shared", "p_value")


def test_asah_roc_matches_proc():
    res = recipe_roc_curve(bench_data("asah"), bench_spec(
        "roc_curve", y="outcome", x="s100b,ndka"))
    st = res["stats"].set_index("marker")
    expect_ref(st.loc["s100b", "auc"], "asah_roc", "auc_s100b")
    expect_ref(st.loc["s100b", "auc_ci_low"], "asah_roc", "auc_low_s100b")
    expect_ref(st.loc["s100b", "auc_ci_high"], "asah_roc", "auc_high_s100b")
    expect_ref(st.loc["ndka", "auc"], "asah_roc", "auc_ndka")
    expect_ref(st.loc["ndka", "delong_p_vs_first"], "asah_roc", "delong_p")


def test_pefr_bland_altman_matches_the_paper():
    res = recipe_bland_altman(bench_data("pefr"), bench_spec(
        "bland_altman", x="wright_1", y="mini_1"))
    s = res["stats"].iloc[0]
    expect_ref(s["bias"], "pefr_bland_altman", "bias_published")
    expect_ref(s["sd_diff"], "pefr_bland_altman", "sd_published")
    expect_ref(s["bias"], "pefr_bland_altman", "bias")
    expect_ref(s["loa_low"], "pefr_bland_altman", "loa_low")
    expect_ref(s["loa_high"], "pefr_bland_altman", "loa_high")


def test_bcg_meta_analysis_matches_metafor():
    res = recipe_forest_meta(bench_data("bcg"), bench_spec(
        "forest_meta", y="risk_ratio", ci_low="ci_low", ci_high="ci_high", label="trial"))
    assert res["test_meta"]["measure"] == "ratio"
    st = res["stats"].set_index("label")
    re = st.loc["Random effects (DerSimonian-Laird)"]
    expect_ref(re["estimate"], "bcg_meta", "random_estimate")
    expect_ref(re["ci_low"], "bcg_meta", "random_ci_low")
    expect_ref(re["ci_high"], "bcg_meta", "random_ci_high")
    expect_ref(re["tau2"], "bcg_meta", "tau2")
    expect_ref(re["i2"], "bcg_meta", "i2")
    expect_ref(re["q"], "bcg_meta", "q")
    expect_ref(st.loc["Fixed effect", "estimate"], "bcg_meta", "fixed_estimate")


def check_mixed(res, case, group, time):
    """Compare a mixed ANOVA with its references.

    Args:
        res: The recipe result.
        case: The benchmark case.
        group: The name of the between-subjects effect.
        time: The name of the within-subjects effect.
    """
    st = res["stats"].set_index("effect")
    expect_ref(st.loc[group, "F"], case, "F_group")
    expect_ref(st.loc[group, "p_value"], case, "p_group")
    expect_ref(st.loc[time, "F"], case, "F_time")
    expect_ref(st.loc[f"{group}:{time}", "F"], case, "F_interaction")
    expect_ref(st.loc[f"{group}:{time}", "p_value"], case, "p_interaction")
    expect_ref(st.loc[time, "gg_epsilon"], case, "gg_epsilon")
    expect_ref(res["test_meta"]["sphericity"]["p_value"], case, "mauchly_p")


def test_orthodont_mixed_anova():
    res = recipe_timecourse(bench_data("orthodont"), bench_spec(
        "timecourse", x="age", y="distance", id="Subject", group="Sex"))
    check_mixed(res, "orthodont_timecourse", "Sex", "age")
    assert res["test_meta"]["correction_applied"] is False


def test_chick_weight_mixed_anova_drops_incomplete_chicks():
    res = recipe_timecourse(bench_data("chick_weight"), bench_spec(
        "timecourse", x="Time", y="weight", id="Chick", group="Diet"))
    expect_ref(res["test_meta"]["n_subjects"], "chick_weight_timecourse", "n_complete_chicks")
    check_mixed(res, "chick_weight_timecourse", "Diet", "Time")
    assert res["test_meta"]["correction_applied"] is True


def test_aml_swimmer_arm_summaries():
    res = recipe_swimmer(bench_data("aml"), bench_spec(
        "swimmer", id="patient", time="time", event="status", group="x"))
    arms = res["test_meta"]["arms"]
    expect_ref(arms["Maintained"]["median_duration"], "aml_swimmer", "median_maintained")
    expect_ref(arms["Nonmaintained"]["median_duration"], "aml_swimmer", "median_nonmaintained")
    expect_ref(arms["Maintained"]["ongoing"], "aml_swimmer", "ongoing_maintained")
    expect_ref(arms["Nonmaintained"]["ongoing"], "aml_swimmer", "ongoing_nonmaintained")


def test_airway_volcano_counts_match_deseq2():
    res = recipe_volcano(bench_data("airway_de"), bench_spec("volcano", y="padj"))
    for q in ("n_up", "n_down", "n_total"):
        expect_ref(res["test_meta"][q], "airway_volcano", q)


def test_airway_heatmap_welch_tests():
    spec = bench_spec("heatmap", annotation=os.path.join(BENCH, "data", "airway_annotation.csv"))
    res = recipe_heatmap(bench_data("airway_heatmap"), spec)
    st = res["stats"].set_index("feature")
    for g in ("SPARCL1", "CACNB2", "DUSP1", "SAMHD1", "MAOA"):
        expect_ref(st.loc[g, "p_value"], "airway_heatmap", f"p_{g}")


def test_airway_pca_variance_and_exact_permanova():
    res = recipe_pca(bench_data("airway_pca"), bench_spec("pca", group="dex"))
    expect_ref(res["test_meta"]["pc1_variance_pct"], "airway_pca", "pc1_variance_pct")
    expect_ref(res["test_meta"]["pc2_variance_pct"], "airway_pca", "pc2_variance_pct")
    # Eight samples in two groups of four give 70 labelings, so the test is exact.
    # The observed split and its mirror are the two most extreme labelings: p = 2/70.
    perm = res["test_meta"]["permanova"]
    assert perm["method"] == "exact"
    assert perm["p_value"] == pytest.approx(2 / 70)


def test_airway_enrichment_top_term():
    res = recipe_enrichment_dot(bench_data("airway_enrichment"), bench_spec("enrichment_dot"))
    top = res["stats"].iloc[0]
    expect_ref(top["count"], "airway_enrichment", "top_count")
    expect_ref(top["p_adjust"], "airway_enrichment", "top_p_adjust")
    expect_ref(top["gene_ratio"], "airway_enrichment", "top_gene_ratio")


def test_airway_upset_intersections():
    res = recipe_upset(bench_data("airway_upset"), bench_spec("upset"))
    for q in ("n_elements", "largest_intersection", "n_intersections"):
        expect_ref(res["test_meta"][q], "airway_upset", q)


def test_anorexia_family_therapy_paired_t():
    res = recipe_paired_compare(bench_data("anorexia_ft"), bench_spec(
        "paired_compare", x="condition", y="weight_lbs", id="subject", test="t"))
    s = res["stats"].iloc[0]
    expect_ref(s["statistic"], "anorexia_paired", "statistic")
    expect_ref(s["p_value"], "anorexia_paired", "p_value")


def test_veteran_log_rank_by_cell_type():
    res = recipe_survival_km(bench_data("veteran"), bench_spec(
        "survival_km", time="time", event="status", x="celltype"))
    expect_ref(res["test_meta"]["log_rank"]["chisq"], "veteran_km", "log_rank_chisq")
    expect_ref(res["test_meta"]["log_rank"]["p_value"], "veteran_km", "log_rank_p")


def test_veteran_cox_hazard_ratios():
    res = recipe_cox_forest(bench_data("veteran"), bench_spec(
        "cox_forest", time="time", event="status", covariates="trt,celltype,karno,age"))
    st = res["stats"].set_index("term")
    for term in ("trttest", "celltypelarge", "celltypesmallcell", "celltypesquamous",
                 "karno", "age"):
        expect_ref(st.loc[term, "HR"], "veteran_cox", f"HR_{term}")
        expect_ref(st.loc[term, "p_value"], "veteran_cox", f"p_{term}")


def test_state_x77_correlation_matrix_keeps_spaced_names():
    res = recipe_correlation_heatmap(bench_data("state_x77"), bench_spec("correlation_heatmap"))
    lm = pair_row(res["stats"], "Murder", "Life Exp", "var1", "var2")
    expect_ref(lm["r"], "state_correlation_heatmap", "r_murder_life_exp")
    expect_ref(lm["p_adj"], "state_correlation_heatmap", "padj_murder_life_exp")
    expect_ref(pair_row(res["stats"], "Illiteracy", "Frost", "var1", "var2")["r"],
               "state_correlation_heatmap", "r_illiteracy_frost")


def test_airway_ma_plot_counts():
    """Check the airway MA plot counts against the DESeq2 table."""
    res = recipe_ma_plot(bench_data("airway_de"), bench_spec("ma_plot"))
    tm = res["test_meta"]
    assert tm["mean_column"] == "baseMean"
    assert tm["p_column"] == "padj"
    for q in ("n_up", "n_down", "n_total"):
        expect_ref(tm[q], "airway_ma", q)
    assert len(res["stats"]) == tm["n_up"] + tm["n_down"]


def test_gwas_results_manhattan_lambda_and_counts():
    """Check lambda_GC and the counts of the qqman GWAS data."""
    res = recipe_manhattan(bench_data("gwas_results"), bench_spec(
        "manhattan", group="CHR", x="BP", y="P", label="SNP"))
    for q in ("lambda_gc", "n_genome_wide", "n_suggestive", "min_p"):
        expect_ref(res["test_meta"][q], "gwas_manhattan", q)
    genome = res["stats"][res["stats"]["row_type"] == "genome"].iloc[0]
    expect_ref(genome["lambda_gc"], "gwas_manhattan", "lambda_gc")
    expect_ref(genome["p_value"], "gwas_manhattan", "min_p")


IRIS_FEATURES = ("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")


def test_iris_feature_panel_kruskal_with_bh():
    """Check the iris Kruskal-Wallis tests and their BH adjustment."""
    res = recipe_feature_panel(bench_data("iris"), bench_spec("feature_panel", x="Species"))
    assert res["test_meta"]["test"] == "kruskal"
    st = res["stats"].set_index("feature")
    for f in IRIS_FEATURES:
        expect_ref(st.loc[f, "statistic"], "iris_feature_panel_kruskal", f"H_{f}")
        expect_ref(st.loc[f, "p_value"], "iris_feature_panel_kruskal", f"p_{f}")
        expect_ref(st.loc[f, "p_adj"], "iris_feature_panel_kruskal", f"p_adj_{f}")


def test_iris_feature_panel_anova_with_bh():
    """Check the iris one-way ANOVAs and their BH adjustment."""
    res = recipe_feature_panel(bench_data("iris"), bench_spec("feature_panel", x="Species",
                                                              test="anova"))
    st = res["stats"].set_index("feature")
    for f in IRIS_FEATURES:
        expect_ref(st.loc[f, "statistic"], "iris_feature_panel_anova", f"F_{f}")
        expect_ref(st.loc[f, "p_adj"], "iris_feature_panel_anova", f"p_adj_{f}")


def test_iris_embedding_silhouettes():
    """Check the iris silhouettes against cluster::silhouette."""
    res = recipe_embedding(bench_data("iris"), bench_spec(
        "embedding", x="Petal.Length", y="Petal.Width", group="Species"))
    st = res["stats"].set_index("cluster")
    for sp in ("setosa", "versicolor", "virginica"):
        expect_ref(st.loc[sp, "mean_silhouette"], "iris_embedding", f"silhouette_{sp}")
    expect_ref(res["test_meta"]["mean_silhouette"], "iris_embedding", "silhouette_overall")
