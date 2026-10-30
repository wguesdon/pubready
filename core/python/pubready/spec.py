"""The figure spec: one dict that both `plot` (CLI args) and `render` (config)
build, then hand to a recipe. Mirrors the R spec structure."""
import os
import re

import yaml


def default_appearance():
    return {
        "theme": "pubready_house",
        "palette": None,
        "x_label": None,
        "y_label": None,
        "title": None,
        "y_limits": [None, None],
        "show_points": True,
        "geom": "box",
        "bracket": {"show": True, "label": "p.signif"},
        "scale": "row",
        "cluster": "both",
        "fc_cutoff": 1.0,
        "p_cutoff": 0.05,
        "top_n": 15,
        "pd_threshold": 20.0,
        "pr_threshold": -30.0,
    }


#: The text options of figkit plot and figkit diagnose, in the order of the R
#: list (spec.R plot_option_list).
PLOT_TEXT_OPTIONS = ("x", "y", "fill", "facet", "time", "event", "covariates", "annotation",
                     "id", "group", "label", "scale", "cluster")
#: The other options with a value, and their types.
PLOT_VALUE_OPTIONS = (("fc_cutoff", float), ("p_cutoff", float), ("top_n", int),
                      ("pd_threshold", float), ("pr_threshold", float), ("positive", str),
                      ("se", str), ("ci_low", str), ("ci_high", str), ("measure", str),
                      ("markers", str), ("ref", str), ("cause", str), ("levels", str),
                      ("visit_levels", str), ("features", str), ("test", str),
                      ("geom", str), ("theme", str), ("xlab", str),
                      ("ylab", str), ("title", str), ("palette", str), ("sheet", str),
                      ("stamp", str))


def add_plot_arguments(p):
    """Add the options of figkit plot and figkit diagnose to a parser.

    One list serves both commands, so diagnose accepts every option that plot
    accepts. An option that the user does not give stays None, and
    spec_from_args() fills its default. The R engine (spec.R) defines the same
    options.

    Args:
        p: An argparse.ArgumentParser.

    Returns:
        argparse.ArgumentParser: The parser.
    """
    p.add_argument("--recipe", required=True)
    p.add_argument("--data", required=True)
    for name in PLOT_TEXT_OPTIONS:
        p.add_argument(f"--{name}")
    for name, kind in PLOT_VALUE_OPTIONS:
        p.add_argument(f"--{name}", type=kind)
    p.add_argument("--paired", action="store_true")
    p.add_argument("--engine")  # consumed by the figkit dispatcher; ignored here
    p.add_argument("--out", default="pubready_output")
    return p


def split_list_option(x):
    """Split a list option such as --levels into its items.

    The items are split on "|" when the value holds one, so an item can hold a
    comma, for example --levels "edema, no diuretics|edema despite diuretics".
    Else they are split on ",". The R engine uses the same rule.

    Args:
        x: The option value, or None.

    Returns:
        list or None: The stripped items, or None when x is None or empty.
    """
    if not x:
        return None
    sep = "|" if "|" in x else ","
    return [s.strip() for s in x.split(sep)]


def spec_from_args(a):
    ap = default_appearance()
    if getattr(a, "geom", None):
        ap["geom"] = a.geom
    if getattr(a, "theme", None):
        ap["theme"] = a.theme
    if getattr(a, "xlab", None):
        ap["x_label"] = a.xlab
    if getattr(a, "ylab", None):
        ap["y_label"] = a.ylab
    if getattr(a, "title", None):
        ap["title"] = a.title
    if getattr(a, "palette", None):
        ap["palette"] = [s.strip() for s in a.palette.split(",")]
    if getattr(a, "scale", None):
        ap["scale"] = a.scale
    if getattr(a, "cluster", None):
        ap["cluster"] = a.cluster
    if getattr(a, "fc_cutoff", None) is not None:
        ap["fc_cutoff"] = a.fc_cutoff
    if getattr(a, "p_cutoff", None) is not None:
        ap["p_cutoff"] = a.p_cutoff
    if getattr(a, "top_n", None) is not None:
        ap["top_n"] = a.top_n
    if getattr(a, "pd_threshold", None) is not None:
        ap["pd_threshold"] = a.pd_threshold
    if getattr(a, "pr_threshold", None) is not None:
        ap["pr_threshold"] = a.pr_threshold
    covs = split_list_option(getattr(a, "covariates", None))
    return {
        "engine": "python",
        "recipe": a.recipe,
        "data": {
            "file": os.path.basename(a.data), "x": a.x, "y": a.y,
            "fill": getattr(a, "fill", None), "facet": getattr(a, "facet", None),
            "time": getattr(a, "time", None), "event": getattr(a, "event", None),
            "annotation": getattr(a, "annotation", None), "covariates": covs,
            "id": getattr(a, "id", None), "group": getattr(a, "group", None),
            "label": getattr(a, "label", None),
            "positive": getattr(a, "positive", None),
            "se": getattr(a, "se", None),
            "ci_low": getattr(a, "ci_low", None),
            "ci_high": getattr(a, "ci_high", None),
            "measure": getattr(a, "measure", None),
            "ref": getattr(a, "ref", None),
            "cause": getattr(a, "cause", None),
            "levels": split_list_option(getattr(a, "levels", None)),
            "visit_levels": split_list_option(getattr(a, "visit_levels", None)),
            "features": split_list_option(getattr(a, "features", None)),
            "markers": split_list_option(getattr(a, "markers", None)),
        },
        "test": {
            "method": getattr(a, "test", None) or "auto",
            "paired": bool(getattr(a, "paired", False)),
        },
        "appearance": ap,
    }


def spec_from_config(path):
    cfg = yaml.safe_load(open(path))
    ap = default_appearance()
    ap.update(cfg.get("appearance") or {})
    cfg["appearance"] = ap
    t = {"method": "auto", "paired": False}
    t.update(cfg.get("test") or {})
    # --p_adjust was removed: no recipe applied it. An old config drops the key.
    t.pop("p_adjust", None)
    cfg["test"] = t
    cfg.setdefault("engine", "python")
    cfg.setdefault("data", {})
    return cfg


def resolve_config_path(path, config):
    """The path of a file that a plot_config.yaml names.

    A relative path is found in the folder of the config first, so a bundle
    renders from any folder, and else in the current folder, as before. The R
    engine has the same rule.

    Args:
        path: The path in the config, or None.
        config: The path of the config.

    Returns:
        str or None: The path to read.
    """
    if not path or re.match(r"^(/|~|[A-Za-z]:)", path):
        return path
    near = os.path.join(os.path.dirname(config), path)
    return near if os.path.exists(near) else path


def write_config(spec, path):
    with open(path, "w") as f:
        yaml.safe_dump(spec, f, sort_keys=False, default_flow_style=False)
