"""The artifact bundle writer for the Python engine. Produces the same bundle
layout as the R engine, tagged engine = "python"."""
import copy
import datetime
import json
import os
import platform
import re
import shutil
from importlib.metadata import distributions
from importlib.metadata import version as _pkgver

import numpy as np

from .io import file_checksums
from .spec import write_config
from .version import PUBPLOT_VERSION


def _now():
    return datetime.datetime.now(datetime.timezone.utc)


def _ts(stamp=None):
    return stamp if stamp else _now().strftime("%Y_%m_%d_%H%M%S")


def _created():
    return os.environ.get("PUBPLOT_CREATED") or _now().strftime("%Y-%m-%dT%H:%M:%SZ")


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s).lower()).strip("_")


def pkg_versions(names):
    out = {}
    for n in names:
        try:
            out[n] = _pkgver(n)
        except Exception:
            out[n] = None
    return out


def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and np.isnan(o):
        return None
    return str(o)


def save_figure(fig, d, stub, w, h):
    fig.set_size_inches(w, h)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(d, f"{stub}.{ext}"),
                    dpi=300, bbox_inches="tight", facecolor="white")


def write_data_log(path, raw, in_name, sums, df_used, steps):
    types = "\n".join(f"- `{c}`: {df_used[c].dtype}" for c in df_used.columns)
    if isinstance(steps, (list, tuple)):
        steps = "\n".join(f"- {s}" for s in steps)
    else:
        steps = f"- {steps}"
    txt = (
        "# Data provenance log\n\n"
        f"- Original file: `{os.path.basename(raw)}`\n"
        f"- Stored in bundle as: `{in_name}`\n"
        f"- md5: `{sums['md5']}`\n"
        f"- sha256: `{sums['sha256']}`\n"
        f"- Analyzed rows: {len(df_used)}\n"
        f"- Columns: {df_used.shape[1]}\n\n"
        "## Column types (analyzed data)\n\n"
        f"{types}\n\n"
        "## Cleaning steps\n\n"
        f"{steps}\n"
    )
    open(path, "w").write(txt)


def write_session_info(path):
    lines = sorted(
        f"{d.metadata['Name']}=={d.version}"
        for d in distributions() if d.metadata.get("Name")
    )
    header = f"Python {platform.python_version()} ({platform.platform()})\n\nInstalled packages:\n"
    open(path, "w").write(header + "\n".join(lines) + "\n")


def build_manifest(spec, in_name, orig, sums, n_rows, test_meta, container, git_commit, created, fig_stub):
    return {
        "pubplot_version": PUBPLOT_VERSION,
        "pubplot_git_commit": git_commit,
        "created_utc": created,
        "engine": "python",
        "recipe": spec["recipe"],
        "recipe_version": "1",
        "arguments": {k: spec["data"].get(k) for k in
                      ("x", "y", "fill", "facet", "time", "event", "covariates")},
        "input": {"file": in_name, "original_name": orig, "md5": sums["md5"],
                  "sha256": sums["sha256"], "n_rows": n_rows},
        "statistical_test": test_meta,
        "container": container,
        "environment": {
            "language": "Python",
            "language_version": platform.python_version(),
            "packages": pkg_versions(["pandas", "numpy", "scipy", "matplotlib",
                                      "seaborn", "pingouin", "statannotations"]),
        },
        "outputs": [f"{fig_stub}.pdf", f"{fig_stub}.png", f"{fig_stub}.svg"],
    }


def reproduce_md(container, git_commit, script_name):
    img = container.get("image", "localhost/pubplot:0.4.1")
    return (
        "# Reproduce this figure\n\n"
        "This bundle is self-contained. To regenerate the figure in the exact\n"
        "environment it was made in:\n\n"
        "## 1. Get the pinned environment\n\n"
        f"- Container image: `{img}`\n"
        f"- Image id: `{container.get('image_id', 'unknown')}`\n"
        f"- pubplot commit: `{git_commit}`\n\n"
        "## 2. Run the standalone script\n\n"
        "From inside this folder:\n\n"
        "```bash\n"
        f'podman run --rm -v "$PWD":/work -w /work \\\n'
        f"  {img} python3 {script_name}\n"
        "```\n\n"
        "The script reads the bundled input copy, reruns the analysis, and redraws\n"
        "the figure. See `manifest_*.json` for package versions and\n"
        "`session_info.txt` for the full environment.\n"
    )


def write_bundle(spec, result, raw_input, out_root, container, git_commit, stamp=None):
    ts = _ts(stamp)
    created = _created()
    base = slugify(result.get("label") or spec["data"].get("y") or
                   spec["data"].get("time") or spec["data"].get("x") or spec["recipe"])
    bdir = os.path.join(out_root, f"{spec['recipe']}_{base}_{ts}")
    os.makedirs(bdir, exist_ok=True)

    fig_stub = f"figure_{base}_{ts}"
    save_figure(result["fig"], bdir, fig_stub, result.get("width", 3.8), result.get("height", 4.0))

    result["stats"].to_csv(os.path.join(bdir, f"stats_{ts}.csv"), index=False)

    ext = os.path.splitext(raw_input)[1].lstrip(".") or "csv"
    in_name = f"input_{base}.{ext}"
    shutil.copyfile(raw_input, os.path.join(bdir, in_name))
    sums = file_checksums(raw_input)

    for ex in (result.get("extra_inputs") or []):
        if ex and os.path.exists(ex):
            shutil.copyfile(ex, os.path.join(bdir, os.path.basename(ex)))

    write_data_log(os.path.join(bdir, "data_log.md"), raw_input, in_name, sums,
                   result["df_used"], result["clean_steps"])

    cfg = copy.deepcopy(spec)
    cfg["data"]["file"] = in_name
    write_config(cfg, os.path.join(bdir, "plot_config.yaml"))

    write_session_info(os.path.join(bdir, "session_info.txt"))

    manifest = build_manifest(spec, in_name, os.path.basename(raw_input), sums,
                              len(result["df_used"]), result["test_meta"], container,
                              git_commit, created, fig_stub)
    with open(os.path.join(bdir, f"manifest_{ts}.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=_json_default)

    open(os.path.join(bdir, f"methods_{ts}.md"), "w").write(result["methods"].rstrip() + "\n")

    script_name = f"script_{ts}.py"
    open(os.path.join(bdir, script_name), "w").write(result["build_script"](in_name, fig_stub))

    open(os.path.join(bdir, "REPRODUCE.md"), "w").write(reproduce_md(container, git_commit, script_name))

    return bdir
