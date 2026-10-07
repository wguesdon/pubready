from types import SimpleNamespace

import pandas as pd

from pubready.recipes.enrichment_dot import recipe_enrichment_dot
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
