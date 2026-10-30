"""Drive one recipe end to end for the Python engine."""
import os
import sys

import numpy as np
import pandas as pd

from . import bundle
from .io import read_tidy_logged
from .recipes import get_recipe
from .util import level_labels, level_order

#: The errors that a user can correct in the data or the command. An entry
#: script prints their message only, as R prints a stop() message; any other
#: exception keeps its traceback, because it is a defect of pubready.
USER_ERRORS = (ValueError, FileNotFoundError)


def run_entry(main):
    """Run the main function of an entry script.

    Args:
        main: The main function, which takes no argument.

    Raises:
        SystemExit: With status 1 and the message only, for a user error.
    """
    try:
        main()
    except USER_ERRORS as err:
        print(f"Error: {err}", file=sys.stderr)
        sys.exit(1)


def container_info():
    return {
        "image": os.environ.get("PUBREADY_IMAGE", "localhost/pubready:3.0.0"),
        "image_id": os.environ.get("PUBREADY_IMAGE_ID", ""),
        "digest": os.environ.get("PUBREADY_IMAGE_DIGEST", ""),
        "podman_version": os.environ.get("PUBREADY_PODMAN_VERSION", ""),
    }


def run_recipe(spec, raw_input, out_root, sheet=None, stamp=None):
    # The manifest and plot_config.yaml record the engine that ran, also when a
    # render reads a config of the other engine.
    spec["engine"] = "python"
    df, read_steps = read_tidy_logged(raw_input, sheet)
    fn = get_recipe(spec["recipe"])
    result = fn(df, spec)
    bdir = bundle.write_bundle(spec, result, raw_input, out_root, container_info(),
                               os.environ.get("PUBREADY_GIT_COMMIT", "unknown"), stamp,
                               read_steps=read_steps)
    print(f"Wrote bundle: {bdir}")
    return bdir


#: The inspect rules, the same in the R engine (run.R). A text or logical
#: column with 1 to 20 distinct values, and a numeric column with 2 to 6, lists
#: its levels in the order of level_order(), each with its count of rows. A
#: column with 2 to 6 levels is a candidate group, and a numeric one is marked
#: as a numeric code. A column is a possible subject id when each of its values
#: occurs at most once in each level of another column, half of its values or
#: more occur in two of those levels or more, and it holds more distinct values
#: than that column. The other column holds 2 to 20 levels and repeats each
#: level at least twice on average; the first 10 such columns are checked. A
#: numeric column with more than 6 distinct values that is not a possible
#: subject id is a candidate outcome.
INSPECT_GROUP_MAX = 6
INSPECT_LEVELS_MAX = 20
INSPECT_BY_MAX = 10


def inspect_type(s):
    """The type name of a column in the inspect report.

    Args:
        s: The column.

    Returns:
        str: numeric, text, logical, date or other, as the R engine names them.
    """
    if pd.api.types.is_bool_dtype(s):
        return "logical"
    if pd.api.types.is_numeric_dtype(s):
        return "numeric"
    if pd.api.types.is_datetime64_any_dtype(s):
        return "date"
    if s.dtype == object or isinstance(s.dtype, pd.CategoricalDtype) or \
            pd.api.types.is_string_dtype(s):
        present = s.dropna()
        if len(present) and all(isinstance(v, (bool, np.bool_)) for v in present):
            return "logical"
        return "text"
    return "other"


def inspect_labels(s):
    """The text labels of a column for the level and id rules.

    Args:
        s: The column.

    Returns:
        pandas.Series: level_labels() of each value, with NaN for a missing one.
    """
    if pd.api.types.is_datetime64_any_dtype(s):
        return s.dt.strftime("%Y-%m-%d").where(s.notna(), np.nan)
    if isinstance(s.dtype, pd.CategoricalDtype):
        s = s.astype(object)
    return level_labels(s)


def inspect_has_levels(kind, n_unique):
    """True when a column of the given type and distinct count lists its levels."""
    return ((kind in ("text", "logical") and 1 <= n_unique <= INSPECT_LEVELS_MAX) or
            (kind == "numeric" and 2 <= n_unique <= INSPECT_GROUP_MAX))


def inspect_column(name, s):
    """One column of the inspect report.

    Args:
        name: The column name.
        s: The column.

    Returns:
        dict: The type, the counts, the levels and the numeric summary.
    """
    kind = inspect_type(s)
    present = inspect_labels(s).dropna()
    n_unique = int(present.nunique())
    info = {"name": name, "type": kind, "n_missing": int(s.isna().sum()), "n_unique": n_unique}
    if inspect_has_levels(kind, n_unique):
        lv = level_order(list(present))
        counts = present.value_counts()
        info["levels"] = lv
        info["level_counts"] = [int(counts[v]) for v in lv]
    if kind == "numeric":
        values = s.to_numpy(dtype=float, na_value=np.nan)
        if np.isfinite(values).any():
            info["min"] = float(np.nanmin(values))
            info["max"] = float(np.nanmax(values))
            info["mean"] = float(np.nanmean(values))
    return info


def inspect_id_candidate(a_lab, b_lab, a_name, b_name):
    """The id record of column a across the levels of column b.

    Args:
        a_lab: The labels of column a.
        b_lab: The labels of column b.
        a_name: The name of column a.
        b_name: The name of column b.

    Returns:
        dict or None: The record, or None when a is not a possible subject id
        across b.
    """
    keep = a_lab.notna() & b_lab.notna()
    pair = pd.DataFrame({"a": a_lab[keep].to_numpy(), "b": b_lab[keep].to_numpy()})
    n_a = pair["a"].nunique()
    if n_a < 2 or n_a <= pair["b"].nunique():
        return None
    if pair.duplicated().any():
        return None
    tab = pd.crosstab(pair["a"], pair["b"])
    n_levels = (tab > 0).sum(axis=1)
    if (n_levels >= 2).mean() < 0.5:
        return None
    return {"column": a_name, "across": b_name, "n_subjects": int(tab.shape[0]),
            "complete": bool((n_levels == tab.shape[1]).all())}


def inspect_data(df):
    """Inspect a table for the figkit inspect command.

    Args:
        df: The table, as read_tidy reads it.

    Returns:
        dict: The row and column counts, one record for each column, and the
        candidate groups, outcomes and subject ids.
    """
    names = [str(c) for c in df.columns]
    cols = [inspect_column(nm, df.iloc[:, i]) for i, nm in enumerate(names)]
    labs = [inspect_labels(df.iloc[:, i]) for i in range(len(names))]
    by = [j for j, c in enumerate(cols)
          if 2 <= c["n_unique"] <= INSPECT_LEVELS_MAX and c["n_unique"] * 2 <= len(df)]
    by = by[:INSPECT_BY_MAX]
    ids = []
    for i in range(len(cols)):
        for j in by:
            if i == j:
                continue
            cand = inspect_id_candidate(labs[i], labs[j], names[i], names[j])
            if cand is not None:
                ids.append(cand)
                break
    id_names = {c["column"] for c in ids}
    groups = [c["name"] for c in cols
              if "levels" in c and 2 <= c["n_unique"] <= INSPECT_GROUP_MAX]
    outcomes = [c["name"] for c in cols if c["type"] == "numeric"
                and c["n_unique"] > INSPECT_GROUP_MAX and c["name"] not in id_names]
    return {"n_rows": int(len(df)), "n_cols": int(df.shape[1]), "columns": cols,
            "group_candidates": groups, "outcome_candidates": outcomes, "id_candidates": ids}


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
    for c in info["columns"][:INSPECT_MAX_COLS]:
        line = f"  {c['name']:<20} {c['type']:<8}  missing={c['n_missing']}  unique={c['n_unique']}"
        if "levels" in c:
            line += "  levels: " + ", ".join(
                f'"{lv}" ({n})' for lv, n in zip(c["levels"], c["level_counts"]))
        if "mean" in c:
            line += f"  [{_g3(c['min'])} .. {_g3(c['max'])}, mean {_g3(c['mean'])}]"
        print(line)
    hidden = len(info["columns"]) - INSPECT_MAX_COLS
    if hidden > 0:
        print(f"  ... and {hidden} more column(s); use --format json for every column.")
    print()
    types = {c["name"]: c["type"] for c in info["columns"]}
    groups = [f"{g} (numeric code)" if types.get(g) == "numeric" else g
              for g in info["group_candidates"]]
    if groups:
        print(f"{'Candidate grouping columns (x):':<33}{inspect_names(groups)}")
    if info["outcome_candidates"]:
        print(f"{'Candidate outcome columns (y):':<33}{inspect_names(info['outcome_candidates'])}")
    if info["id_candidates"]:
        ids = [f"{c['column']} (repeats across the levels of {c['across']})"
               for c in info["id_candidates"]]
        print(f"{'Possible subject id columns:':<33}{inspect_names(ids)}")


def _g3(x):
    """Format a number with 3 significant digits, as R sprintf("%.3g") does."""
    if np.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if np.isnan(x):
        return "NaN"
    return f"{x:.3g}"
