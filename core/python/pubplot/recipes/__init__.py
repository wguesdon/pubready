"""Recipe registry for the Python engine."""
import importlib


def get_recipe(name):
    try:
        mod = importlib.import_module(f"pubplot.recipes.{name}")
    except ModuleNotFoundError:
        raise ValueError(f"unknown recipe '{name}'")
    fn = getattr(mod, f"recipe_{name}", None)
    if fn is None:
        raise ValueError(f"recipe '{name}' has no recipe_{name} function")
    return fn
