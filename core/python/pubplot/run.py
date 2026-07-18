"""Drive one recipe end to end for the Python engine."""
import os

from . import bundle
from .io import read_tidy
from .recipes import get_recipe


def container_info():
    return {
        "image": os.environ.get("PUBPLOT_IMAGE", "localhost/pubplot:0.1.0"),
        "image_id": os.environ.get("PUBPLOT_IMAGE_ID", ""),
        "digest": os.environ.get("PUBPLOT_IMAGE_DIGEST", ""),
        "podman_version": os.environ.get("PUBPLOT_PODMAN_VERSION", ""),
    }


def run_recipe(spec, raw_input, out_root, sheet=None, stamp=None):
    df = read_tidy(raw_input, sheet)
    fn = get_recipe(spec["recipe"])
    result = fn(df, spec)
    bdir = bundle.write_bundle(spec, result, raw_input, out_root, container_info(),
                               os.environ.get("PUBPLOT_GIT_COMMIT", "unknown"), stamp)
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


def print_inspect(info):
    print(f"Rows: {info['n_rows']}   Columns: {info['n_cols']}\n")
    groups, nums = [], []
    for c in info["columns"]:
        line = f"  {c['name']:<20} {c['type']:<10}  missing={c['n_missing']}  unique={c['n_unique']}"
        if "levels" in c:
            line += "  levels: " + ", ".join(c["levels"])
            if 2 <= c["n_unique"] <= 6:
                groups.append(c["name"])
        if "mean" in c:
            line += f"  [{c['min']:.3g} .. {c['max']:.3g}, mean {c['mean']:.3g}]"
            nums.append(c["name"])
        print(line)
    print()
    if groups:
        print("Candidate grouping columns (x): ", ", ".join(groups))
    if nums:
        print("Candidate outcome columns  (y): ", ", ".join(nums))
