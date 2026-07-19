from types import SimpleNamespace

import pandas as pd

from pubplot.recipes.volcano import recipe_volcano
from pubplot.spec import spec_from_args


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
