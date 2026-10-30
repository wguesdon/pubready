import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubready.recipes.heatmap import heatmap_stars, recipe_heatmap
from pubready.spec import spec_from_args


def test_heatmap_stars_follow_the_r_thresholds():
    assert heatmap_stars([1e-5, 0.005, 0.03, 0.2, np.nan]) == ["***", "**", "*", "", ""]


def test_heatmap_draws_stars_and_keeps_legends_clear_of_the_row_names(tmp_path):
    rng = np.random.RandomState(3)
    samples = [f"a{i}" for i in range(6)] + [f"b{i}" for i in range(6)]
    values = rng.normal(0, 1, (10, 12))
    values[:3, 6:] += 4  # three features differ between the groups
    df = pd.DataFrame(values, columns=samples)
    df.insert(0, "gene", [f"gene_long_name_{i:02d}" for i in range(10)])
    ann_path = tmp_path / "ann.csv"
    pd.DataFrame({"sample": samples, "group": ["a"] * 6 + ["b"] * 6}).to_csv(ann_path, index=False)
    a = SimpleNamespace(recipe="heatmap", data="matrix.csv", x=None, y=None,
                        annotation=str(ann_path), test="welch_t")
    res = recipe_heatmap(df, spec_from_args(a))

    fig = res["fig"]
    assert res["test_meta"]["engine_lib"] == "PyComplexHeatmap"
    stars = [t.get_text() for ax in fig.axes for t in ax.texts if set(t.get_text()) == {"*"}]
    assert len(stars) == int(res["stats"]["significant"].sum()) >= 3

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    heat = [ax for ax in fig.axes
            if any(t.get_text().startswith("gene_") for t in ax.get_yticklabels())]
    assert len(heat) == 1
    names_right = max(t.get_window_extent(renderer).x1 for t in heat[0].get_yticklabels())
    heat_right = heat[0].get_window_extent(renderer).x1
    # The legend axes hold a legend or the colour bar; the empty right axes of
    # PyComplexHeatmap hold nothing.
    legends = [ax.get_window_extent(renderer).x0 for ax in fig.axes
               if ax.get_window_extent(renderer).x0 > heat_right
               and (ax.get_legend() is not None or ax.collections or ax.images)]
    assert legends, "no legend axes to the right of the heatmap"
    assert min(legends) > names_right


# The heatmap leaf order is the same in both engines. The R suite
# (test-heatmap.R) uses the same matrix and expects the same order. g1 and g3
# form one branch and g2 and g4 the other. The branch g2, g4 has the higher
# mean, and inside each branch the gene with the higher mean comes first. The
# columns s1 and s3 are equal, as are s2 and s4, so the order inside each column
# branch falls to the tie rule: the earlier input column first.
def hm_matrix():
    """The four-gene matrix of the R suite.

    Returns:
        pandas.DataFrame: The matrix, with the gene names in the first column.
    """
    return pd.DataFrame({"gene": ["g1", "g2", "g3", "g4"],
                         "s1": [1, 10, 1.5, 10.2], "s2": [2, 11, 2.5, 11.2],
                         "s3": [1, 10, 1.5, 10.2], "s4": [2, 11, 2.5, 11.2]})


def hm_spec(cluster="both"):
    """A heatmap spec with unscaled values.

    Args:
        cluster: both, rows, columns or none.

    Returns:
        dict: The spec.
    """
    return spec_from_args(SimpleNamespace(recipe="heatmap", data="m.csv", x=None, y=None,
                                          scale="none", cluster=cluster))


def test_heatmap_orders_each_branch_by_the_mean_of_its_unscaled_values():
    res = recipe_heatmap(hm_matrix(), hm_spec())
    assert res["test_meta"]["row_order"] == ["g4", "g2", "g3", "g1"]
    assert res["test_meta"]["column_order"] == ["s2", "s4", "s1", "s3"]
    fig = res["fig"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    heat = [ax for ax in fig.axes if any(t.get_text() == "g1" for t in ax.get_yticklabels())]
    assert len(heat) == 1
    top_down = sorted((t for t in heat[0].get_yticklabels() if t.get_text()),
                      key=lambda t: -t.get_window_extent(renderer).y0)
    assert [t.get_text() for t in top_down] == ["g4", "g2", "g3", "g1"]


def test_an_axis_that_is_not_clustered_keeps_the_input_order():
    res = recipe_heatmap(hm_matrix(), hm_spec("rows"))
    assert res["test_meta"]["row_order"] == ["g4", "g2", "g3", "g1"]
    assert res["test_meta"]["column_order"] == ["s1", "s2", "s3", "s4"]


def test_the_emitted_script_carries_the_order_rule_and_runs(tmp_path):
    df = hm_matrix()
    res = recipe_heatmap(df, hm_spec())
    script = res["build_script"]("input.csv", "figure")
    assert "def ordered_linkage(" in script
    df.to_csv(tmp_path / "input.csv", index=False)
    (tmp_path / "script.py").write_text(script)
    run = subprocess.run([sys.executable, "script.py"], cwd=tmp_path, capture_output=True,
                         text=True, check=False)
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "figure.pdf").exists()


def test_heatmap_methods_name_the_groups_the_independence_and_the_count():
    from pubready.recipes.heatmap import _methods
    txt = _methods("row", "both", 40, 8, ("trt", "untrt"), 38, True)
    assert "compared between trt and untrt with a Welch two-sample t-test" in txt
    assert "treats the samples as independent" in txt
    assert "38 of the 40 features had an adjusted P below 0.05" in txt
    assert "feature(s)" not in _methods("row", "both", 40, 8, None, None, True)


def test_heatmap_legends_are_not_clipped_to_the_legend_axes(tmp_path):
    from matplotlib.legend import Legend
    samples = [f"a{i}" for i in range(4)] + [f"b{i}" for i in range(4)]
    rng = np.random.RandomState(1)
    df = pd.DataFrame(rng.normal(0, 1, (6, 8)), columns=samples)
    df.insert(0, "gene", [f"g{i}" for i in range(6)])
    ann_path = tmp_path / "ann.csv"
    pd.DataFrame({"sample": samples, "group": ["a"] * 4 + ["b"] * 4,
                  "donor": ["N052611", "N061011", "N080611", "N61311"] * 2}).to_csv(
        ann_path, index=False)
    a = SimpleNamespace(recipe="heatmap", data="matrix.csv", x=None, y=None,
                        annotation=str(ann_path))
    res = recipe_heatmap(df, spec_from_args(a))
    legends = [c for ax in res["fig"].axes for c in ax.get_children() if isinstance(c, Legend)]
    assert legends and not any(lg.get_clip_on() for lg in legends)


def hm_null(tmp_path):
    """A matrix and an annotation with no difference between the two groups.

    Args:
        tmp_path: The folder for the annotation file.

    Returns:
        tuple: The matrix and the path of the annotation.
    """
    df = pd.DataFrame({"gene": [f"g{i}" for i in range(1, 5)],
                       "a1": [5.1, 3.0, 7.2, 1.1], "a2": [4.9, 3.3, 6.8, 0.9],
                       "a3": [5.3, 2.8, 7.0, 1.2], "b1": [5.0, 3.1, 7.1, 1.0],
                       "b2": [5.2, 2.9, 6.9, 1.3], "b3": [4.8, 3.2, 7.3, 0.8]})
    ann_path = tmp_path / "ann.csv"
    pd.DataFrame({"sample": ["a1", "a2", "a3", "b1", "b2", "b3"],
                  "group": ["a"] * 3 + ["b"] * 3}).to_csv(ann_path, index=False)
    return df, ann_path


def test_a_null_heatmap_draws_no_star_column(tmp_path):
    """A null result gets no star column in either engine (F46)."""
    df, ann_path = hm_null(tmp_path)
    a = SimpleNamespace(recipe="heatmap", data="matrix.csv", x=None, y=None,
                        annotation=str(ann_path), test="welch_t")
    res = recipe_heatmap(df, spec_from_args(a))
    assert res["test_meta"]["n_significant"] == 0
    stars = [t for ax in res["fig"].axes for t in ax.texts if set(t.get_text()) == {"*"}]
    assert not stars


def test_the_heatmap_runs_no_test_by_default(tmp_path):
    """The Welch test runs only with --test welch_t (F97)."""
    import pytest
    df, ann_path = hm_null(tmp_path)
    a = SimpleNamespace(recipe="heatmap", data="matrix.csv", x=None, y=None,
                        annotation=str(ann_path))
    res = recipe_heatmap(df, spec_from_args(a))
    assert list(res["stats"].columns) == ["feature", "mean", "sd"]
    assert res["test_meta"]["per_feature_test"] is None
    assert "No statistical test was applied." in res["methods"]
    three = pd.read_csv(ann_path)
    three.loc[0, "group"] = "c"
    three.to_csv(ann_path, index=False)
    a.test = "welch_t"
    with pytest.raises(ValueError, match="compares two groups"):
        recipe_heatmap(df, spec_from_args(a))
    a.test = "anova"
    with pytest.raises(ValueError, match="unknown --test 'anova'"):
        recipe_heatmap(df, spec_from_args(a))
