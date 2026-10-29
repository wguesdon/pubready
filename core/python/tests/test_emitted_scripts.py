"""The emitted scripts of the comparison recipes draw the geom of the figure, and
they run. Each script runs in a temporary folder, as a user runs it from a
bundle. The R suite (test-emitted-scripts.R) runs the same cases."""
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.multi_group_compare import recipe_multi_group_compare
from pubready.recipes.two_group_compare import recipe_two_group_compare
from pubready.spec import spec_from_args

# The call that only the given geom emits.
MARKERS = {"box": "sns.boxplot(", "violin": "sns.violinplot(", "bar": "ax.bar(",
           "raincloud": "ax.violinplot("}

RECIPES = {"two_group_compare": (recipe_two_group_compare, 2),
           "multi_group_compare": (recipe_multi_group_compare, 3)}


def es_spec(recipe, geom):
    """Build a comparison spec as the plot entry does.

    Args:
        recipe: The recipe name.
        geom: The geom.

    Returns:
        dict: The spec.
    """
    return spec_from_args(SimpleNamespace(recipe=recipe, data="d.csv", x="group", y="value",
                                          test="auto", geom=geom))


def es_data(k):
    """Eight values in each of k groups, the same values as in the R suite.

    Args:
        k: The count of groups, 2 or 3.

    Returns:
        pandas.DataFrame: The tidy table.
    """
    i = np.arange(1, 9)
    values = np.concatenate([10 + np.sin(i), 14 + np.cos(i), 18 + np.sin(2 * i)])[: 8 * k]
    return pd.DataFrame({"group": np.repeat(["ctrl", "low", "high"][:k], 8), "value": values})


def run_script(script, df, folder):
    """Write the input and the script into a folder and run the script there.

    Args:
        script: The emitted script.
        df: The input table.
        folder: The folder.

    Returns:
        subprocess.CompletedProcess: The finished run.
    """
    df.to_csv(folder / "input.csv", index=False)
    (folder / "script.py").write_text(script)
    return subprocess.run([sys.executable, "script.py"], cwd=folder, capture_output=True,
                          text=True, check=False)


@pytest.mark.parametrize("geom", ["box", "violin", "bar", "raincloud"])
@pytest.mark.parametrize("name", ["two_group_compare", "multi_group_compare"])
def test_the_script_draws_the_geom_and_runs(name, geom, tmp_path):
    recipe, k = RECIPES[name]
    df = es_data(k)
    res = recipe(df, es_spec(name, geom))
    script = res["build_script"]("input.csv", "figure")
    assert MARKERS[geom] in script
    run = run_script(script, df, tmp_path)
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "figure.pdf").exists()


def test_a_script_without_points_draws_no_points():
    spec = es_spec("two_group_compare", "box")
    spec["appearance"]["show_points"] = False
    res = recipe_two_group_compare(es_data(2), spec)
    assert "stripplot" not in res["build_script"]("input.csv", "figure")
