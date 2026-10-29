"""Behaviour of the four omics recipes (ma_plot, manhattan, feature_panel,
embedding) on small synthetic inputs. The numbers themselves are checked
against public benchmarks in test_benchmarks.py."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pubready.recipes.embedding import point_size, recipe_embedding, silhouette_widths
from pubready.recipes.feature_panel import bh_adjust, recipe_feature_panel
from pubready.recipes.ma_plot import recipe_ma_plot
from pubready.recipes.manhattan import DOUBLE_XMIN, chr_sort_key, recipe_manhattan, thin_mask
from pubready.spec import spec_from_args


def omics_spec(recipe, **kw):
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


def test_ma_plot_reads_limma_names_and_keeps_missing_p_as_not_significant():
    """Check that the MA plot reads limma names and keeps missing p as not significant."""
    df = pd.DataFrame({"ID": [f"g{i}" for i in range(1, 7)],
                       "logFC": [2, -3, 0.2, 1.5, -1.2, 4], "AveExpr": [5, 6, 7, 8, 9, 10],
                       "adj.P.Val": [1e-4, 1e-6, 0.5, 0.2, 0.01, np.nan]})
    res = recipe_ma_plot(df, omics_spec("ma_plot"))
    tm = res["test_meta"]
    assert tm["mean_column"] == "AveExpr"
    assert tm["p_column"] == "adj.P.Val"
    assert not tm["mean_log10"]
    assert (tm["n_up"], tm["n_down"], tm["n_not_significant"], tm["n_missing_p"]) == (1, 2, 3, 1)
    assert list(res["stats"]["gene"]) == ["g2", "g1", "g5"]
    assert list(res["stats"]["A"]) == [6, 5, 9]


def test_ma_plot_logs_basemean_and_drops_a_zero_mean():
    """Check that the MA plot logs basemean and drops a zero mean."""
    df = pd.DataFrame({"gene": ["a", "b", "c"], "baseMean": [100, 0, 1000],
                       "log2FoldChange": [1.5, 2, -0.1], "padj": [0.01, 0.01, 0.9]})
    res = recipe_ma_plot(df, omics_spec("ma_plot"))
    assert res["test_meta"]["mean_log10"]
    assert res["test_meta"]["n_total"] == 2
    assert res["stats"]["A"].tolist() == [2.0]
    assert any("baseMean of zero" in s for s in res["clean_steps"])


def test_ma_plot_needs_a_mean_column():
    """Check that the MA plot needs a mean column."""
    df = pd.DataFrame({"gene": ["a"], "log2FoldChange": [1.0], "padj": [0.01]})
    with pytest.raises(ValueError, match="no mean-expression column"):
        recipe_ma_plot(df, omics_spec("ma_plot"))


def test_manhattan_natural_order_and_one_lead_per_chromosome():
    """Check the natural chromosome order and one lead SNP per chromosome."""
    df = pd.DataFrame({"snp": [f"rs{i}" for i in range(1, 9)],
                       "chr": ["chrX", "chr10", "chr2", "chr2", "chrMT", "chr1", "chr10", "chr22"],
                       "pos": [5, 10, 20, 30, 1, 7, 40, 3],
                       "p": [0.5, 1e-9, 2e-10, 3e-12, 0.9, 0.04, 2e-9, 0.3]})
    res = recipe_manhattan(df, omics_spec("manhattan", group="chr", x="pos", y="p", label="snp"))
    assert res["test_meta"]["n_chromosomes"] == 6
    leads = res["stats"][res["stats"]["row_type"] == "lead"]
    assert list(leads["snp"]) == ["rs4", "rs2"]
    assert list(leads["chr"]) == ["2", "10"]
    assert sorted(["X", "10", "2", "MT", "1", "22", "Y"], key=chr_sort_key) == \
        ["1", "2", "10", "22", "X", "Y", "MT"]
    assert res["test_meta"]["n_genome_wide"] == 4
    assert any("'chr' prefix" in s for s in res["clean_steps"])


def test_manhattan_zero_p_and_p_above_one():
    """Check the handling of a p of zero and of a p above one."""
    df = pd.DataFrame({"chr": [1, 1, 2, 2], "pos": [1, 2, 3, 4], "p": [0, 0.5, 0.2, 0.9]})
    res = recipe_manhattan(df, omics_spec("manhattan", group="chr", x="pos", y="p"))
    assert res["test_meta"]["n_zero_p"] == 1
    assert res["test_meta"]["min_p"] == DOUBLE_XMIN
    bad = pd.DataFrame({"chr": [1, 1], "pos": [1, 2], "p": [0.5, 1.5]})
    with pytest.raises(ValueError, match="outside 0 to 1"):
        recipe_manhattan(bad, omics_spec("manhattan", group="chr", x="pos", y="p"))


def test_manhattan_thins_only_large_p_and_keeps_lambda_on_all():
    """Check that the Manhattan plot thins only large p and keeps lambda on all SNPs."""
    rng = np.random.default_rng(7)
    n = 60000
    p = rng.uniform(size=n)
    df = pd.DataFrame({"chr": np.arange(n) % 3 + 1, "pos": np.arange(1, n + 1), "p": p})
    res = recipe_manhattan(df, omics_spec("manhattan", group="chr", x="pos", y="p"))
    keep, k = thin_mask(p)
    assert keep[p < 0.01].all()
    assert k > 1 and res["test_meta"]["thinning_step"] == k
    assert res["test_meta"]["n_drawn"] <= 51000
    expected = np.median(stats.chi2.ppf(1 - p, 1)) / stats.chi2.ppf(0.5, 1)
    assert res["test_meta"]["lambda_gc"] == pytest.approx(expected, rel=1e-10)


def test_feature_panel_welch_t_and_bh():
    """Check the Welch t-test and the BH adjustment of the feature panel."""
    rng = np.random.default_rng(11)
    df = pd.DataFrame({"group": ["ctrl"] * 12 + ["trt"] * 12, "f1": rng.normal(size=24),
                       "f2": rng.normal(size=24) + np.repeat([0, 2], 12), "f3": rng.normal(size=24)})
    res = recipe_feature_panel(df, omics_spec("feature_panel", x="group", test="welch_t"))
    st = res["stats"].set_index("feature")
    tt = stats.ttest_ind(df["f2"][:12], df["f2"][12:], equal_var=False)
    assert st.loc["f2", "statistic"] == pytest.approx(tt.statistic, rel=1e-10)
    assert st.loc["f2", "effect_size_name"] == "Cohen's d"
    assert np.allclose(res["stats"]["p_adj"], bh_adjust(res["stats"]["p_value"]))
    assert list(res["stats"]["n_ctrl"]) == [12, 12, 12]


def test_bh_adjust_matches_the_rules_of_p_adjust():
    """Check that bh adjust matches the rules of p adjust."""
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    q = bh_adjust(p)
    assert np.isnan(q[2])
    assert q[[0, 1, 3, 4]] == pytest.approx([0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])


def test_feature_panel_long_data_with_the_exact_wilcoxon_rule():
    """Check long data and the exact Wilcoxon rule of the feature panel."""
    df = pd.DataFrame({"gene": ["A"] * 10 + ["B"] * 10,
                       "arm": (["x"] * 5 + ["y"] * 5) * 2,
                       "value": list(range(1, 11)) + [1, 3, 5, 7, 9, 2, 4, 6, 8, 10]})
    res = recipe_feature_panel(df, omics_spec("feature_panel", x="arm", y="value", facet="gene",
                                              test="wilcoxon"))
    exact = stats.mannwhitneyu([1, 3, 5, 7, 9], [2, 4, 6, 8, 10], method="exact").pvalue
    assert res["stats"].set_index("feature").loc["B", "p_value"] == pytest.approx(exact, rel=1e-12)
    assert list(res["stats"]["feature"]) == ["A", "B"]


def test_feature_panel_refuses_too_many_features_and_a_wrong_test():
    """Check that feature panel refuses too many features and a wrong test."""
    rng = np.random.default_rng(1)
    wide = pd.DataFrame(rng.normal(size=(6, 37)), columns=[f"V{i}" for i in range(37)])
    wide["group"] = ["a", "b"] * 3
    with pytest.raises(ValueError, match="heatmap"):
        recipe_feature_panel(wide, omics_spec("feature_panel", x="group"))
    three = pd.DataFrame({"group": ["a"] * 4 + ["b"] * 4 + ["c"] * 4, "v": rng.normal(size=12)})
    with pytest.raises(ValueError, match="exactly 2 groups"):
        recipe_feature_panel(three, omics_spec("feature_panel", x="group", test="welch_t"))


def test_silhouette_widths_follow_the_rules_of_cluster_silhouette():
    """Check that silhouette widths follow the rules of cluster silhouette."""
    xy = np.array([[0.0, 0.0], [0.0, 1.0], [5.0, 0.0], [9.0, 9.0]])
    codes = np.array([0, 0, 1, 2])
    s = silhouette_widths(xy, codes, 3)
    a0, b0 = 1.0, 5.0
    assert s[0] == pytest.approx((b0 - a0) / b0)
    assert s[2] == 0.0 and s[3] == 0.0
    assert silhouette_widths(xy, np.zeros(4, dtype=int), 1) is None


def test_embedding_point_size_matches_r():
    """Check that a point has the diameter of the ggplot2 point of the R engine.

    R draws shape 16 as a circle of radius 0.375 times size * 72.27 / 25.4 points.
    """
    rng = np.random.default_rng(7)
    d = pd.DataFrame({"x": rng.normal(size=200), "y": rng.normal(size=200),
                      "k": np.repeat(["a", "b"], 100)})
    res = recipe_embedding(d, omics_spec("embedding", x="x", y="y", group="k"))
    diameter = 0.75 * 72.27 / 25.4 * point_size(200)
    assert res["fig"].axes[0].collections[0].get_sizes()[0] == pytest.approx(diameter ** 2)


def test_embedding_ids_and_scores():
    """Check that whole-number ids are clusters and a score is continuous."""
    rng = np.random.default_rng(6)
    ids = pd.DataFrame({"x": rng.normal(size=200), "y": rng.normal(size=200),
                        "k": np.tile(np.arange(20), 10)})
    res = recipe_embedding(ids, omics_spec("embedding", x="x", y="y", group="k"))
    assert not res["test_meta"]["continuous"]
    assert list(res["stats"]["cluster"][:3]) == ["0", "1", "2"]
    assert res["stats"]["percent"].sum() == pytest.approx(100)
    score = pd.DataFrame({"x": rng.normal(size=200), "y": rng.normal(size=200),
                          "s": rng.uniform(size=200)})
    res2 = recipe_embedding(score, omics_spec("embedding", x="x", y="y", group="s"))
    assert res2["test_meta"]["continuous"]
    assert np.isnan(res2["stats"]["mean_silhouette"].iloc[0])


def test_embedding_subsamples_the_silhouette_above_5000_cells():
    """Check that embedding subsamples the silhouette above 5000 cells."""
    rng = np.random.default_rng(8)
    n = 12001
    df = pd.DataFrame({"x": rng.normal(size=n), "y": rng.normal(size=n),
                       "g": rng.choice(["p", "q"], size=n)})
    res = recipe_embedding(df, omics_spec("embedding", x="x", y="y", group="g"))
    assert res["test_meta"]["silhouette_step"] == 3
    assert res["test_meta"]["silhouette_n"] == len(range(0, n, 3))
    assert "one in every 3 cells" in res["methods"]


def test_kruskal_effect_size_is_eta_squared_h_in_both_recipes():
    from pubready.recipes.multi_group_compare import recipe_multi_group_compare
    df = pd.DataFrame({"group": ["a"] * 8 + ["b"] * 8 + ["c"] * 8,
                       "f1": list(range(1, 9)) + [v + 3 for v in range(1, 9)]
                       + [v + 6 for v in range(1, 9)]})
    h = stats.kruskal(*[df.loc[df["group"] == g, "f1"] for g in "abc"]).statistic
    want = (h - 3 + 1) / (24 - 3)
    fp = recipe_feature_panel(df, omics_spec("feature_panel", x="group", test="kruskal"))
    assert fp["stats"]["effect_size"].iloc[0] == pytest.approx(want, rel=1e-12)
    assert fp["stats"]["effect_size_name"].iloc[0] == "eta-squared (H)"
    mg = recipe_multi_group_compare(df.rename(columns={"f1": "value"}), spec_from_args(
        SimpleNamespace(recipe="multi_group_compare", data="d.csv", x="group", y="value",
                        test="kruskal")))
    assert mg["test_meta"]["omnibus"]["effect_size"]["value"] == pytest.approx(want, rel=1e-12)
