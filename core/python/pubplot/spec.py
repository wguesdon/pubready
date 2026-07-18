"""The figure spec: one dict that both `plot` (CLI args) and `render` (config)
build, then hand to a recipe. Mirrors the R spec structure."""
import os
import yaml


def default_appearance():
    return {
        "theme": "pubplot_house",
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
    }


def spec_from_args(a):
    ap = default_appearance()
    if getattr(a, "geom", None):
        ap["geom"] = a.geom
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
    covs = [s.strip() for s in a.covariates.split(",")] if getattr(a, "covariates", None) else None
    return {
        "engine": "python",
        "recipe": a.recipe,
        "data": {
            "file": os.path.basename(a.data), "x": a.x, "y": a.y,
            "fill": getattr(a, "fill", None), "facet": getattr(a, "facet", None),
            "time": getattr(a, "time", None), "event": getattr(a, "event", None),
            "annotation": getattr(a, "annotation", None), "covariates": covs,
        },
        "test": {
            "method": getattr(a, "test", None) or "auto",
            "paired": bool(getattr(a, "paired", False)),
            "p_adjust": getattr(a, "p_adjust", None) or "none",
        },
        "appearance": ap,
    }


def spec_from_config(path):
    cfg = yaml.safe_load(open(path))
    ap = default_appearance()
    ap.update(cfg.get("appearance") or {})
    cfg["appearance"] = ap
    t = {"method": "auto", "paired": False, "p_adjust": "none"}
    t.update(cfg.get("test") or {})
    cfg["test"] = t
    cfg.setdefault("engine", "python")
    cfg.setdefault("data", {})
    return cfg


def write_config(spec, path):
    with open(path, "w") as f:
        yaml.safe_dump(spec, f, sort_keys=False, default_flow_style=False)
