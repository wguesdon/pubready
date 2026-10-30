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


def enrich_area_d():
    """One term with an adjusted p of 0, one with none and one with no name, as in R.

    Returns:
        pandas.DataFrame: The table.
    """
    return pd.DataFrame({"Description": ["Pathway 0", "Pathway 1", None, "Pathway 3", "Pathway 4"],
                         "GeneRatio": ["5/200", "6/200", "7/200", "8/200", "9/200"],
                         "Count": [5, 6, 7, 8, 9],
                         "p.adjust": [0.0, float("nan"), 1e-4, 1e-5, 1e-6]})


def enrich_spec():
    """An enrichment_dot spec.

    Returns:
        dict: The spec.
    """
    return spec_from_args(SimpleNamespace(recipe="enrichment_dot", data="e.csv", x=None, y=None,
                                          label=None, top_n=15))


def test_enrichment_dot_keeps_a_zero_p_and_logs_the_dropped_rows():
    """An adjusted p of 0 is kept and marked; each dropped row is logged (F51)."""
    import numpy as np
    res = recipe_enrichment_dot(enrich_area_d(), enrich_spec())
    assert res["test_meta"]["n_terms_total"] == 3
    assert res["test_meta"]["n_zero_p"] == 1
    assert res["stats"]["term"].iloc[0] == "Pathway 0"
    assert res["stats"]["p_adjust"].iloc[0] == np.finfo(float).tiny
    steps = " | ".join(res["clean_steps"])
    assert "Dropped 2 row(s) with a missing term, adjusted p or x value (5 -> 3)." in steps
    assert "Replaced 1 adjusted p value(s) of exactly 0" in steps
    assert "drawn as triangles" in res["methods"]


def test_enrichment_dot_labels_a_raw_p_as_not_adjusted():
    """A table with only a raw p is labelled as not adjusted, as in R (F52)."""
    df = enrich_area_d().iloc[[0, 3, 4]].rename(columns={"p.adjust": "pvalue"})
    df["pvalue"] = [1e-3, 1e-5, 1e-6]
    res = recipe_enrichment_dot(df, enrich_spec())
    assert res["test_meta"]["p_adjusted"] is False
    assert "p_value" in res["stats"].columns and "p_adjust" not in res["stats"].columns
    assert "the raw p-value, which is not adjusted for multiple testing" in res["methods"]
    assert "raw p column 'pvalue'" in res["clean_steps"][0]
    titles = [ax.get_title(loc="left") for ax in res["fig"].axes]
    assert "-log10 p (not adjusted)" in titles
