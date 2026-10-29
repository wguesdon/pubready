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
                        annotation=str(ann_path))
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
