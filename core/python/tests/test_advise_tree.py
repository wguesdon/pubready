"""Check tree integrity, documented leaves and the recipe test contracts."""

import copy
import json
import os
from pathlib import Path
import re
from types import SimpleNamespace

import pytest

from pubready.advise import load_tree, validate_tree
from pubready.io import read_tidy
from pubready.recipes import get_recipe
from pubready.spec import read_recipe_options, spec_from_args

TREE = load_tree()[0]
CORE = Path(__file__).resolve().parents[2]
FIXTURES = json.loads((CORE / "advise_test_fixtures.json").read_text())
LEAVES = {n["leaf"]["recipe"]: n["leaf"] for n in TREE["nodes"].values() if "leaf" in n}
TEST_CASES = [(name, method) for name, leaf in LEAVES.items()
              for method in leaf.get("test", {}).get("allowed", [])]


def test_tree_and_docs():
    """Require all recipes and every heading to remain documented."""
    table = read_recipe_options()
    assert set(LEAVES) == set(table) - {"always"}
    doc = Path("/opt/pubready/reference/decision_tree.md").read_text()
    headings = set(re.findall(r"^#+ (.+)$", doc, re.M))
    heading_recipes = {name for heading in headings
                       for name in re.findall(r"`([a-z_]+)`", heading)}
    assert heading_recipes <= set(LEAVES)
    parts = re.split(r"^#+ (.+)$", doc, flags=re.M)
    sections = dict(zip(parts[1::2], parts[2::2]))
    for node in TREE["nodes"].values():
        content = node.get("leaf", node)
        assert content["doc_heading"] in headings
    for name, leaf in LEAVES.items():
        assert callable(get_recipe(name))
        assert f"`{name}`" in doc
        if leaf.get("test", {}).get("allowed"):
            values = [v for v in leaf["test"]["allowed"] if v != "auto"]
            section = sections[leaf["doc_heading"]]
            forced = re.findall(r"Force with\s+`--test ([^`]+)`", section)
            assert forced == [" | ".join(values)], name
    assert set(FIXTURES) == {n for n, leaf in LEAVES.items() if leaf.get("test")}


def test_cycles_unreachable_nodes_and_missing_branches():
    """Reject unsafe or incomplete trees before a walk starts."""
    tree = copy.deepcopy(TREE)
    tree["nodes"]["start"]["branches"]["compare"] = "start"
    with pytest.raises(ValueError, match="cycle"):
        validate_tree(tree)
    tree = copy.deepcopy(TREE)
    tree["nodes"]["unused"] = {"out_of_scope": "Unused."}
    with pytest.raises(ValueError, match="unreachable"):
        validate_tree(tree)
    tree = copy.deepcopy(TREE)
    del tree["nodes"]["start"]["branches"]["compare"]
    with pytest.raises(ValueError, match="coverage"):
        validate_tree(tree)


def fixture_spec(recipe, method):
    """Load a shared fixture and build its recipe specification.

    Args:
        recipe: The recipe name.
        method: The canonical test value.

    Returns:
        tuple: Input data and specification.
    """
    fixture = FIXTURES[recipe]
    base = Path(os.environ["PUBREADY_BENCHMARKS"]) / "data"
    flags = dict(fixture["flags"])
    if "annotation" in flags:
        flags["annotation"] = str(base / flags["annotation"])
    args = dict(recipe=recipe, data=str(base / (fixture["dataset"] + ".csv")),
                x=None, y=None, test=method)
    args.update(flags)
    df = read_tidy(args["data"])
    if recipe == "feature_panel" and method in ("welch_t", "wilcoxon"):
        df = df[df["Species"].isin(df["Species"].unique()[:2])]
    return df, spec_from_args(SimpleNamespace(**args))


@pytest.mark.parametrize("recipe,method", TEST_CASES)
def test_canonical_test_contract(recipe, method):
    """Execute every advertised test and verify the exact reported outcome."""
    df, spec = fixture_spec(recipe, method)
    fn = get_recipe(recipe)
    if recipe == "factorial_anova" and method == "art":
        with pytest.raises(ValueError, match="R-only"):
            fn(df, spec)
        return
    result = fn(df, spec)
    allowed = {o["test"] for o in LEAVES[recipe]["test"]["outcomes"]}
    assert result["test_meta"]["name"] in allowed


@pytest.mark.parametrize("recipe", list(FIXTURES))
def test_unknown_tests_fail(recipe):
    """The engine must reject a test absent from its vocabulary."""
    df, spec = fixture_spec(recipe, "unknown_advise_test")
    with pytest.raises(ValueError):
        get_recipe(recipe)(df, spec)
