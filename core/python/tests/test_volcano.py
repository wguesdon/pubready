from types import SimpleNamespace

import pandas as pd

from pubready.recipes.volcano import recipe_volcano
from pubready.spec import spec_from_args


def test_volcano_counts_and_columns():
    df = pd.DataFrame({
        "gene": [f"g{i}" for i in range(6)],
        "log2FoldChange": [3.0, -3.0, 0.1, 2.5, -0.2, 0.0],
        "pvalue": [1e-8, 1e-9, 0.5, 1e-6, 0.3, 0.9],
    })
    a = SimpleNamespace(recipe="volcano", data="de.csv", x=None, y=None, label=None,
                        fc_cutoff=1.0, p_cutoff=0.05, top_n=15)
    spec = spec_from_args(a)

    res = recipe_volcano(df, spec)

    tm = res["test_meta"]
    assert tm["n_up"] == 2      # g0 (+3.0), g3 (+2.5)
    assert tm["n_down"] == 1    # g1 (-3.0)
    assert tm["n_significant"] == 3
    assert set(res["stats"].columns) >= {"gene", "log2FC", "p_value", "direction"}


def test_volcano_legend_skips_empty_classes_and_labels_do_not_overlap():
    # No feature passes the fold change alone, so there is no "FC" key; R drops
    # it too. The labels of the hits must not overlap each other.
    n = 40
    df = pd.DataFrame({
        "gene": [f"Gene{i:05d}" for i in range(n)],
        "log2FoldChange": [2.0 + 0.4 * i for i in range(20)] + [0.1] * 20,
        "pvalue": [10.0 ** -(20 + 3 * i) for i in range(20)] + [0.5] * 20,
    })
    a = SimpleNamespace(recipe="volcano", data="de.csv", x=None, y=None, label=None,
                        fc_cutoff=1.0, p_cutoff=0.05, top_n=15)
    res = recipe_volcano(df, spec_from_args(a))
    ax = res["fig"].axes[0]
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["NS", "FC and p"]
    res["fig"].canvas.draw()
    renderer = res["fig"].canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in ax.texts if t.get_text().startswith("Gene")]
    assert len(boxes) == 15
    for i, b1 in enumerate(boxes):
        for b2 in boxes[i + 1:]:
            assert not b1.overlaps(b2)
