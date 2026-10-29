"""The interaction geom of factorial_anova and the label wrapping of its x axis."""
from types import SimpleNamespace

import numpy as np
import pandas as pd

from pubready.recipes.factorial_anova import fig_size, needs_wrap, recipe_factorial_anova, wrap_x
from pubready.spec import spec_from_args


def fa_spec(**kw):
    """Build a factorial_anova spec as the plot entry does.

    Args:
        **kw: The command-line options.

    Returns:
        dict: The spec.
    """
    base = {"recipe": "factorial_anova", "data": "d.csv", "x": None, "y": None, "test": "anova"}
    base.update(kw)
    return spec_from_args(SimpleNamespace(**base))


def fa_data():
    """A 3 x 2 design with long level names, five values per cell.

    Returns:
        pandas.DataFrame: The data.
    """
    rng = np.random.default_rng(9)
    rows = []
    for strain in ("C57BL/6J wild type", "NOD scid gamma (NSG)", "FVB/N"):
        for dose in ("Vehicle control", "High dose (10 mg/kg)"):
            for _ in range(5):
                rows.append({"strain": strain, "dose": dose,
                             "y": 10 + rng.normal() + (0 if dose == "Vehicle control" else 3)})
    return pd.DataFrame(rows)


def test_interaction_geom_builds_and_keeps_the_box_statistics():
    """Check that interaction geom builds and keeps the box statistics."""
    df = fa_data()
    box = recipe_factorial_anova(df, fa_spec(x="strain", fill="dose", y="y"))
    inter = recipe_factorial_anova(df, fa_spec(x="strain", fill="dose", y="y", geom="interaction"))
    pd.testing.assert_frame_equal(inter["stats"], box["stats"])
    ax = inter["fig"].axes[0]
    assert len(ax.lines) >= 2
    assert len(ax.collections) >= 2
    script = inter["build_script"]("input.csv", "figure")
    assert "0.25 / k" in script


def test_long_x_levels_wrap_at_twelve_and_the_figure_grows():
    """Check that long x levels wrap at twelve and the figure grows."""
    assert wrap_x("NOD scid gamma (NSG)") == "NOD scid\ngamma (NSG)"
    assert not needs_wrap(["wt", "ko"])
    assert fig_size(["wt", "ko"], ["vehicle", "low", "high"], "treatment", 1, False) == (5.2, 4.4)
    w, _ = fig_size([f"strain {i}" for i in range(4)], ["a", "b", "c"], "treatment", 1, False)
    assert w > 5.2
