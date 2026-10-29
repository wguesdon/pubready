from types import SimpleNamespace

import pandas as pd

from pubready.recipes.enrichment_dot import recipe_enrichment_dot, wrap_label
from pubready.spec import spec_from_args


def test_enrichment_dot_ora_mode():
    df = pd.DataFrame({
        "Description": [f"term{i}" for i in range(5)],
        "GeneRatio": ["10/100", "8/100", "6/100", "4/100", "2/100"],
        "Count": [10, 8, 6, 4, 2],
        "p.adjust": [1e-6, 1e-5, 1e-4, 1e-3, 1e-2],
    })
    a = SimpleNamespace(recipe="enrichment_dot", data="e.csv", x=None, y=None, label=None, top_n=15)
    spec = spec_from_args(a)

    res = recipe_enrichment_dot(df, spec)

    assert res["test_meta"]["name"] == "ORA dot plot"
    assert len(res["stats"]) == 5
    assert "gene_ratio" in res["stats"].columns


def test_enrichment_dot_stats_follow_the_adjusted_p():
    # The gene ratio runs against the adjusted p here, so the plot order and the
    # table order differ; the R engine writes the same table order.
    df = pd.DataFrame({
        "Description": [f"term{i}" for i in range(1, 5)],
        "GeneRatio": ["8/100", "6/100", "4/100", "2/100"],
        "Count": [8, 6, 4, 2],
        "p.adjust": [1e-6, 1e-5, 1e-4, 1e-3],
    })
    a = SimpleNamespace(recipe="enrichment_dot", data="e.csv", x=None, y=None, label=None, top_n=15)

    res = recipe_enrichment_dot(df, spec_from_args(a))

    assert list(res["stats"]["term"]) == [f"term{i}" for i in range(1, 5)]
    assert list(res["stats"]["p_adjust"]) == [1e-6, 1e-5, 1e-4, 1e-3]


def test_wrap_label_balances_long_labels():
    # The R engine has the same cases in test-enrichment-dot.R.
    assert wrap_label("glycolytic process") == "glycolytic process"
    assert (wrap_label("regulation of transcription by RNA polymerase II")
            == "regulation of transcription\nby RNA polymerase II")
    assert (wrap_label("negative regulation of extrinsic apoptotic signaling "
                       "pathway via death domain receptors")
            == "negative regulation of extrinsic apoptotic\n"
               "signaling pathway via death domain receptors")


def test_enrichment_dot_legends_sit_right_of_the_panel():
    df = pd.DataFrame({
        "Description": [f"term{i}" for i in range(6)],
        "GeneRatio": [f"{10 - i}/100" for i in range(6)],
        "Count": [10, 8, 6, 4, 2, 1],
        "p.adjust": [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 2e-2],
    })
    a = SimpleNamespace(recipe="enrichment_dot", data="e.csv", x=None, y=None, label=None, top_n=15)
    res = recipe_enrichment_dot(df, spec_from_args(a))
    fig = res["fig"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panel = fig.axes[0]
    legend = panel.get_legend()
    assert legend is not None
    assert legend.get_window_extent(renderer).x0 > panel.get_window_extent(renderer).x1
