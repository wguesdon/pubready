"""The artifact bundle writer for the Python engine. Produces the same bundle
layout as the R engine, tagged engine = "python"."""
import copy
import datetime
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from importlib.metadata import distributions
from importlib.metadata import version as _pkgver

import numpy as np
import pandas as pd

from .io import file_checksums
from .spec import write_config
from .version import PUBREADY_VERSION


def _now():
    return datetime.datetime.now(datetime.timezone.utc)


def _ts(stamp=None):
    return stamp if stamp else _now().strftime("%Y_%m_%d_%H%M%S")


def _created():
    return os.environ.get("PUBREADY_CREATED") or _now().strftime("%Y-%m-%dT%H:%M:%SZ")


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


def json_safe(o):
    """Make a value valid JSON, as the R engine writes it.

    NaN and infinity become None, which JSON writes as null, and the numpy
    scalars and arrays become the Python types.

    Args:
        o: The value: a dict, a list or a scalar.

    Returns:
        The value with each part made valid JSON.
    """
    if isinstance(o, dict):
        return {k: json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return json_safe(o.tolist())
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if o is None or isinstance(o, (int, str)):
        return o
    if o is pd.NA or o is pd.NaT:
        return None
    return _json_default(o)


def write_json(obj, path):
    """Write a value as indented, valid JSON, with null for NaN and infinity.

    Args:
        obj: The value.
        path: The output file.
    """
    with open(path, "w") as f:
        json.dump(json_safe(obj), f, indent=2, allow_nan=False, default=_json_default)


def save_figure(fig, d, stub, w, h):
    fig.set_size_inches(w, h)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(os.path.join(d, f"{stub}.{ext}"),
                    dpi=300, bbox_inches="tight", facecolor="white")


def _markdown_value(value):
    """Escape one Markdown table value.

    Args:
        value: The value to format.

    Returns:
        The escaped text.
    """
    if value is None:
        return ""
    return str(value).replace("|", r"\|").replace("\n", " ")


def integer_columns(frame):
    """Write each float column that holds only whole numbers as integers.

    pandas stores a count or a df column as float when the column holds a
    missing value, and to_csv then writes 16470.0 where R writes 16470. The
    nullable Int64 type keeps the missing value and drops the ".0".

    Args:
        frame: The stats table.

    Returns:
        pandas.DataFrame: A copy with the whole-number float columns as Int64.
    """
    out = frame.copy()
    for col in out.columns:
        values = out[col]
        if values.dtype.kind != "f":
            continue
        present = values.dropna().to_numpy()
        if (len(present) and np.all(np.isfinite(present)) and np.all(np.abs(present) < 1e15)
                and np.all(present == np.round(present))):
            out[col] = values.astype("Int64")
    return out


def dataframe_markdown(frame, limit=100):
    """Convert a data frame to a small Markdown table.

    Args:
        frame: The pandas data frame.
        limit: The maximum row count.

    Returns:
        A Markdown table.
    """
    shown = frame.head(limit)
    headers = [_markdown_value(column) for column in shown.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in shown.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(_markdown_value(value) for value in row) + " |")
    return "\n".join(lines)


def _git_dirty():
    """Read the recorded repository dirty state."""
    value = os.environ.get("PUBREADY_GIT_DIRTY", "not_recorded")
    if value in {"true", "false"}:
        return value == "true"
    return value


def _report_provenance(spec, container, git_commit, created, in_name, sums):
    """Build the shared report provenance record.

    Args:
        spec: The resolved recipe specification.
        container: The container record.
        git_commit: The source commit.
        created: The creation time.
        in_name: The bundled input name.
        sums: The input checksums.

    Returns:
        The provenance record.
    """
    return {
        "repository": "pubready",
        "repository_version": PUBREADY_VERSION,
        "git_commit": git_commit,
        "git_dirty": _git_dirty(),
        "llm_provider": os.environ.get("PUBREADY_LLM_PROVIDER", "none"),
        "llm_model": os.environ.get("PUBREADY_LLM_MODEL", "none"),
        "agent_harness": os.environ.get("PUBREADY_AGENT_HARNESS", "manual"),
        "agent_harness_version": os.environ.get(
            "PUBREADY_AGENT_HARNESS_VERSION", "not_recorded"
        ),
        "container_image": container.get("image", "unknown"),
        "container_digest": container.get("digest") or "not_recorded",
        "input": {"file": in_name, **sums},
        "command": f"figkit plot --recipe {spec['recipe']}",
        "analysis_options": spec,
        "random_seeds": spec.get("seed", "not_applicable"),
        "session_record": "session_info.txt",
        "created_utc": created,
        "review_state": "not_applicable",
        "reviewer": "none",
        "quarto_version": "1.11.4",
    }


def _provenance_frame(provenance):
    """Convert provenance to a two column data frame.

    Args:
        provenance: The report provenance record.

    Returns:
        A pandas data frame.
    """
    values = []
    for value in provenance.values():
        if isinstance(value, (dict, list)):
            values.append(json.dumps(json_safe(value), allow_nan=False, default=_json_default))
        else:
            values.append(value)
    return pd.DataFrame({"field": provenance.keys(), "value": values})


def write_quarto_report(bdir, spec, stats_df, methods_text, fig_stub, provenance):
    """Write and render one self contained Quarto report.

    Args:
        bdir: The bundle directory.
        spec: The resolved recipe specification.
        stats_df: The saved result values.
        methods_text: The exact methods paragraph.
        fig_stub: The figure file stem.
        provenance: The shared provenance record.

    Returns:
        The report output names.

    Raises:
        RuntimeError: If Quarto does not render the report.
    """
    write_json(provenance, os.path.join(bdir, "report_provenance.json"))
    model = json.dumps(provenance["llm_model"])
    harness = json.dumps(provenance["agent_harness"])
    if provenance["llm_model"] == "none":
        callout = "A manual command generated this report without an LLM model."
    else:
        callout = (
            f"An LLM generated this report through the `{provenance['agent_harness']}` "
            "agent harness."
        )
    qmd = f'''---
title: "Publication figure report"
author:
  - name: {model}
    affiliations:
      - name: {harness}
format:
  html:
    embed-resources: true
    toc: true
---

::: {{.callout-important}}
{callout}
No separate human scientific review is part of this fixed recipe.
The report records the selected recipe and does not claim a new statistical method.
:::

## Figure

![Publication figure]({fig_stub}.png)

## Statistical results

{dataframe_markdown(stats_df)}

## Methods

{methods_text.strip()}

## Plot configuration

The complete plot configuration is in `plot_config.yaml`.

## Provenance

{dataframe_markdown(_provenance_frame(provenance))}

The report source reads only files in this bundle.
'''
    with open(os.path.join(bdir, "report.qmd"), "w") as stream:
        stream.write(qmd)
    completed = subprocess.run(
        ["quarto", "render", "report.qmd", "--to", "html", "--output", "report.html"],
        cwd=bdir,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Quarto did not create report.html: {completed.stderr.strip()}")
    return ["report.qmd", "report.html", "report_provenance.json"]


def output_checksums(bdir, names):
    """Calculate SHA256 checksums for bundle outputs.

    Args:
        bdir: The bundle directory.
        names: The relative output names.

    Returns:
        Output checksum records.
    """
    records = []
    for name in names:
        path = os.path.join(bdir, name)
        if not os.path.isfile(path):
            continue
        digest = hashlib.sha256()
        with open(path, "rb") as stream:
            block = stream.read(1024 * 1024)
            while block:
                digest.update(block)
                block = stream.read(1024 * 1024)
        records.append({"file": name, "bytes": os.path.getsize(path), "sha256": digest.hexdigest()})
    return records


#: The read rule of io.py, as the data log states it. The R engine writes the
#: same text.
READ_RULE_TEXT = (
    "Each cell was read as text, and the spaces around each value were removed. "
    "Only an empty cell or \"NA\" was read as missing. A column was read as numeric "
    "when each present cell is a number with \".\" as the decimal mark, and as "
    "logical when each present cell is TRUE or FALSE."
)


def write_data_log(path, raw, in_name, sums, df_used, steps, read_steps=None):
    types = "\n".join(f"- `{c}`: {df_used[c].dtype}" for c in df_used.columns)
    if isinstance(steps, (list, tuple)):
        steps = "\n".join(f"- {s}" for s in steps)
    else:
        steps = f"- {steps}"
    read_txt = "\n".join(f"- {s}" for s in (read_steps or ["No value changed at read."]))
    txt = (
        "# Data provenance log\n\n"
        f"- Original file: `{os.path.basename(raw)}`\n"
        f"- Stored in bundle as: `{in_name}`\n"
        f"- md5: `{sums['md5']}`\n"
        f"- sha256: `{sums['sha256']}`\n"
        f"- Analyzed rows: {len(df_used)}\n"
        f"- Columns: {df_used.shape[1]}\n\n"
        "## Read rule\n\n"
        f"{READ_RULE_TEXT}\n\n"
        f"{read_txt}\n\n"
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


def _os_pretty():
    """The OS of the container, from /etc/os-release (e.g. 'Ubuntu 22.04.5 LTS')."""
    try:
        for line in open("/etc/os-release"):
            if line.startswith("PRETTY_NAME="):
                return line.split("=", 1)[1].strip().strip('"')
    except Exception:
        pass
    return None


def script_direct_dists(script_text):
    """Distributions the generated script imports directly, in import order.
    Maps import names to distribution names (e.g. sklearn -> scikit-learn)."""
    import re
    try:
        from importlib.metadata import packages_distributions
        mod_map = packages_distributions()
    except Exception:
        mod_map = {}
    mods = []
    for line in script_text.splitlines():
        s = line.strip()
        m = re.match(r"import\s+(.+)", s)
        if m and not s.startswith("from"):
            for part in m.group(1).split(","):
                name = part.strip().split(" as ")[0].split(".")[0].strip()
                if name:
                    mods.append(name)
        m = re.match(r"from\s+([A-Za-z0-9_.]+)\s+import", s)
        if m:
            mods.append(m.group(1).split(".")[0])
    out = []
    for mod in dict.fromkeys(mods):        # dedupe, keep order
        for dist in mod_map.get(mod, [mod]):
            if dist not in out:
                out.append(dist)
    return out


def write_packages(path, script_text=None):
    """The package versions that produced this figure, next to the figure and
    script: the packages the script imports directly first, then the complete set."""
    seen = {}
    for d in distributions():
        nm = d.metadata.get("Name")
        if nm:
            seen[nm] = d.version
    lines = ["# Package versions for this figure",
             f"# Python {platform.python_version()} | {_os_pretty() or 'unknown OS'}", ""]
    if script_text:
        direct = [n for n in script_direct_dists(script_text) if n in seen]
        if direct:
            lines.append("## Used by this script")
            lines += [f"{n}=={seen[n]}" for n in direct]
            lines.append("")
    lines += ["## Full environment (every installed package)", ""]
    lines += [f"{n}=={v}" for n, v in sorted(seen.items(), key=lambda t: t[0].lower())]
    open(path, "w").write("\n".join(lines) + "\n")


#: The plot options that the manifest records, in this order, with the spec
#: field that holds each one. Each option is recorded, with null when the run
#: did not set it. The R engine records the same options in the same order.
MANIFEST_DATA_ARGS = ("x", "y", "fill", "facet", "time", "event", "covariates",
                      "annotation", "id", "group", "label", "positive", "se",
                      "ci_low", "ci_high", "measure", "markers", "ref", "cause",
                      "levels", "visit_levels", "features")
MANIFEST_APPEARANCE_ARGS = (("scale", "scale"), ("cluster", "cluster"),
                            ("fc_cutoff", "fc_cutoff"), ("p_cutoff", "p_cutoff"),
                            ("top_n", "top_n"), ("pd_threshold", "pd_threshold"),
                            ("pr_threshold", "pr_threshold"), ("geom", "geom"),
                            ("theme", "theme"), ("xlab", "x_label"), ("ylab", "y_label"),
                            ("title", "title"), ("palette", "palette"))
MANIFEST_LIST_ARGS = ("covariates", "markers", "levels", "visit_levels", "features",
                      "palette")


def manifest_arguments(spec):
    """The arguments block of the manifest.

    Args:
        spec: The figure spec.

    Returns:
        dict: Every plot option, named as the command line names it, with None
        for an option that the run did not set.
    """
    data, test, ap = spec.get("data") or {}, spec.get("test") or {}, spec.get("appearance") or {}
    out = {k: data.get(k) for k in MANIFEST_DATA_ARGS}
    out["test"] = test.get("method") or "auto"
    out["paired"] = bool(test.get("paired", False))
    out["p_adjust"] = test.get("p_adjust") or "none"
    for k, field in MANIFEST_APPEARANCE_ARGS:
        out[k] = ap.get(field)
    for k in MANIFEST_LIST_ARGS:
        if out[k] is not None:
            out[k] = [str(v) for v in (out[k] if isinstance(out[k], (list, tuple)) else [out[k]])]
    return out


def build_manifest(spec, in_name, orig, sums, n_rows, test_meta, container, git_commit, created, fig_stub):
    return {
        "pubready_version": PUBREADY_VERSION,
        "pubready_git_commit": git_commit,
        "created_utc": created,
        "engine": "python",
        "recipe": spec["recipe"],
        "recipe_version": "1",
        "arguments": manifest_arguments(spec),
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
    img = container.get("image", "localhost/pubready:3.0.0")
    return (
        "# Reproduce this figure\n\n"
        "This bundle is self-contained. To regenerate the figure in the exact\n"
        "environment it was made in:\n\n"
        "## 1. Get the pinned environment\n\n"
        f"- Container image: `{img}`\n"
        f"- Image id: `{container.get('image_id', 'unknown')}`\n"
        f"- Image digest: `{container.get('digest') or 'n/a (local build)'}`\n"
        f"- pubready commit: `{git_commit}`\n\n"
        "## 2. Run the standalone script\n\n"
        "From inside this folder:\n\n"
        "```bash\n"
        'podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \\\n'
        f"  {img} python3 {script_name}\n"
        "```\n\n"
        "The script reads the bundled input copy, reruns the analysis, and redraws\n"
        "the figure. See `manifest_*.json` for package versions and\n"
        "`session_info.txt` for the full environment.\n"
    )


def write_bundle(spec, result, raw_input, out_root, container, git_commit, stamp=None,
                 read_steps=None):
    ts = _ts(stamp)
    created = _created()
    base = slugify(result.get("label") or spec["data"].get("y") or
                   spec["data"].get("time") or spec["data"].get("x") or spec["recipe"])
    bdir = os.path.join(out_root, f"{spec['recipe']}_{base}_{ts}")
    os.makedirs(bdir, exist_ok=True)

    fig_stub = f"figure_{base}_{ts}"
    save_figure(result["fig"], bdir, fig_stub, result.get("width", 3.8), result.get("height", 4.0))

    # Normality QC panel: a QQ per checked quantity, so the reader can judge the
    # assumption behind the test rather than trust a single Shapiro-Wilk p-value.
    if result.get("qc"):
        import matplotlib.pyplot as plt

        from .qc import qc_normality
        qfig, _ = qc_normality(result["qc"])
        qfig.savefig(os.path.join(bdir, f"qc_normality_{ts}.png"), dpi=200, facecolor="white")
        qfig.savefig(os.path.join(bdir, f"qc_normality_{ts}.pdf"), facecolor="white")
        plt.close(qfig)

    integer_columns(result["stats"]).to_csv(os.path.join(bdir, f"stats_{ts}.csv"), index=False)

    ext = os.path.splitext(raw_input)[1].lstrip(".") or "csv"
    in_name = f"input_{base}.{ext}"
    shutil.copyfile(raw_input, os.path.join(bdir, in_name))
    sums = file_checksums(raw_input)
    script_text = result["build_script"](in_name, fig_stub)

    for ex in (result.get("extra_inputs") or []):
        if ex and os.path.exists(ex):
            shutil.copyfile(ex, os.path.join(bdir, os.path.basename(ex)))

    write_data_log(os.path.join(bdir, "data_log.md"), raw_input, in_name, sums,
                   result["df_used"], result["clean_steps"], read_steps)

    cfg = copy.deepcopy(spec)
    cfg["data"]["file"] = in_name
    write_config(cfg, os.path.join(bdir, "plot_config.yaml"))

    write_session_info(os.path.join(bdir, "session_info.txt"))
    write_packages(os.path.join(bdir, f"packages_{ts}.txt"), script_text)

    manifest = build_manifest(spec, in_name, os.path.basename(raw_input), sums,
                              len(result["df_used"]), result["test_meta"], container,
                              git_commit, created, fig_stub)
    manifest["outputs"].append(f"packages_{ts}.txt")
    if result.get("qc"):
        manifest["outputs"] += [f"qc_normality_{ts}.png", f"qc_normality_{ts}.pdf"]
    provenance = _report_provenance(spec, container, git_commit, created, in_name, sums)
    report_outputs = write_quarto_report(
        bdir,
        spec,
        result["stats"],
        result["methods"],
        fig_stub,
        provenance,
    )
    manifest["outputs"] += report_outputs
    manifest["output_checksums"] = output_checksums(bdir, manifest["outputs"])
    write_json(manifest, os.path.join(bdir, f"manifest_{ts}.json"))

    open(os.path.join(bdir, f"methods_{ts}.md"), "w").write(result["methods"].rstrip() + "\n")

    script_name = f"script_{ts}.py"
    open(os.path.join(bdir, script_name), "w").write(script_text)

    open(os.path.join(bdir, "REPRODUCE.md"), "w").write(reproduce_md(container, git_commit, script_name))

    return bdir
