"""Drive one recipe end to end for the Python engine."""
import os

from . import bundle
from .io import read_tidy_logged
from .recipes import get_recipe


def container_info():
    return {
        "image": os.environ.get("PUBREADY_IMAGE", "localhost/pubready:3.0.0"),
        "image_id": os.environ.get("PUBREADY_IMAGE_ID", ""),
        "digest": os.environ.get("PUBREADY_IMAGE_DIGEST", ""),
        "podman_version": os.environ.get("PUBREADY_PODMAN_VERSION", ""),
    }


def run_recipe(spec, raw_input, out_root, sheet=None, stamp=None):
    df, read_steps = read_tidy_logged(raw_input, sheet)
    fn = get_recipe(spec["recipe"])
    result = fn(df, spec)
    bdir = bundle.write_bundle(spec, result, raw_input, out_root, container_info(),
                               os.environ.get("PUBREADY_GIT_COMMIT", "unknown"), stamp,
                               read_steps=read_steps)
    print(f"Wrote bundle: {bdir}")
    return bdir


def inspect_data(df):
    cols = []
    for nm in df.columns:
        s = df[nm]
        is_num = str(s.dtype).startswith(("int", "float"))
        info = {"name": nm, "type": str(s.dtype), "n_missing": int(s.isna().sum()),
                "n_unique": int(s.dropna().nunique())}
        if not is_num and 0 < info["n_unique"] <= 20:
            info["levels"] = sorted(map(str, s.dropna().unique()))
        if is_num and s.notna().any():
            info["min"] = float(s.min()); info["max"] = float(s.max()); info["mean"] = float(s.mean())
        cols.append(info)
    return {"n_rows": int(len(df)), "n_cols": int(df.shape[1]), "columns": cols}


#: A wide table shows this many columns in full in the text report.
INSPECT_MAX_COLS = 25
#: A candidate list shows this many names before it counts the rest.
INSPECT_MAX_NAMES = 20


def inspect_names(names):
    """Join column names, and count the ones past the display limit.

    Args:
        names: The column names.

    Returns:
        str: The joined names, as the R engine prints them.
    """
    if len(names) <= INSPECT_MAX_NAMES:
        return ", ".join(names)
    return f"{', '.join(names[:INSPECT_MAX_NAMES])}, and {len(names) - INSPECT_MAX_NAMES} more"


def print_inspect(info):
    """Print the text report of figkit inspect.

    A wide table shows its first columns in full and summarises the rest, so a
    gene matrix does not print hundreds of lines; --format json keeps every column.

    Args:
        info: The inspection result of inspect_data.
    """
    print(f"Rows: {info['n_rows']}   Columns: {info['n_cols']}\n")
    groups, nums = [], []
    for i, c in enumerate(info["columns"]):
        if "levels" in c and 2 <= c["n_unique"] <= 6:
            groups.append(c["name"])
        if "mean" in c:
            nums.append(c["name"])
        if i >= INSPECT_MAX_COLS:
            continue
        line = f"  {c['name']:<20} {c['type']:<10}  missing={c['n_missing']}  unique={c['n_unique']}"
        if "levels" in c:
            line += "  levels: " + ", ".join(c["levels"])
        if "mean" in c:
            line += f"  [{c['min']:.3g} .. {c['max']:.3g}, mean {c['mean']:.3g}]"
        print(line)
    hidden = len(info["columns"]) - INSPECT_MAX_COLS
    if hidden > 0:
        print(f"  ... and {hidden} more column(s); use --format json for every column.")
    print()
    if groups:
        print("Candidate grouping columns (x): ", inspect_names(groups))
    if nums:
        print("Candidate outcome columns  (y): ", inspect_names(nums))
