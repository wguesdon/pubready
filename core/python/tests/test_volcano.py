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
    assert set(res["stats"].columns) >= {"gene", "log2FC", "pvalue", "direction"}


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


def test_volcano_methods_name_the_p_column_and_the_label_rule():
    df = pd.DataFrame({"gene": [f"g{i}" for i in range(6)],
                       "log2FoldChange": [3.0, -3.0, 0.1, 2.5, -0.2, 0.0],
                       "padj": [1e-8, 1e-9, 0.5, 1e-6, 0.3, 0.9]})
    a = SimpleNamespace(recipe="volcano", data="de.csv", x=None, y="padj", label=None,
                        fc_cutoff=None, p_cutoff=None, top_n=15)
    res = recipe_volcano(df, spec_from_args(a))
    assert "|log2 fold change| >= 1 and padj <= 0.05" in res["methods"]
    assert "The 3 significant features with the largest product" in res["methods"]


def test_volcano_keeps_zero_p_and_stops_on_p_above_one():
    """A p of 0 is kept as the smallest double; a p above 1 stops the run (F45)."""
    import pytest
    df = pd.DataFrame({"gene": [f"g{i}" for i in range(5)],
                       "log2FoldChange": [5.0, -4.0, 0.1, 2.0, -0.3],
                       "pvalue": [0.0, 0.0, 0.5, 1e-6, 0.4]})
    a = SimpleNamespace(recipe="volcano", data="de.csv", x=None, y=None, label=None,
                        fc_cutoff=1.0, p_cutoff=0.05, top_n=15)
    spec = spec_from_args(a)
    res = recipe_volcano(df, spec)
    tm = res["test_meta"]
    assert tm["n_total"] == 5 and tm["n_zero_p"] == 2
    assert tm["n_up"] == 2 and tm["n_down"] == 1
    assert "Replaced 2 p value(s) of exactly 0" in " ".join(res["clean_steps"])
    assert "drawn as triangles" in res["methods"]
    df.loc[2, "pvalue"] = 1.7
    with pytest.raises(ValueError, match="outside 0 to 1"):
        recipe_volcano(df, spec)


def test_volcano_stats_name_the_p_column_after_the_input():
    """The p column of the stats table keeps the input name (F123)."""
    df = pd.DataFrame({"gene": [f"g{i}" for i in range(1, 5)],
                       "log2FoldChange": [3.0, -3.0, 0.1, 2.0],
                       "padj": [1e-8, 1e-9, 0.5, 1e-6]})
    a = SimpleNamespace(recipe="volcano", data="de.csv", x=None, y="padj", label=None,
                        fc_cutoff=1.0, p_cutoff=0.05, top_n=15)
    res = recipe_volcano(df, spec_from_args(a))
    st = res["stats"]
    assert list(st.columns) == ["gene", "log2FC", "padj", "neg_log10_p", "direction"]
    assert st.loc[st["gene"] == "g1", "padj"].iloc[0] == 1e-8
