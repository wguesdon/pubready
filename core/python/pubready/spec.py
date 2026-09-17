"""The figure spec: one dict that both `plot` (CLI args) and `render` (config)
build, then hand to a recipe. Mirrors the R spec structure."""
import os
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
    covs = [s.strip() for s in a.covariates.split(",")] if getattr(a, "covariates", None) else None
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
