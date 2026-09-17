"""Recipe: UpSet plot of set intersections from a binary membership matrix.
Mirrors the R recipe (ComplexHeatmap UpSet), here via the upsetplot package.
Input is a CSV where each row is an element and each column is a set, with an
optional leading id column. The plot is descriptive (intersection sizes), so
there is no hypothesis test; the stats table lists set sizes and every
intersection size."""
import matplotlib.pyplot as plt
import pandas as pd
import upsetplot

from ..bundle import slugify


def _to_member(s):
    if pd.api.types.is_bool_dtype(s):
        return s.astype(bool)
    sn = pd.to_numeric(s, errors="coerce")
    if sn.notna().all() and set(sn.unique()) <= {0, 1}:
        return sn == 1
    return s.astype(str).str.lower().isin(["1", "true", "yes", "y", "t"])


def recipe_upset(df, spec):
    first = df.columns[0]
    first_is_id = not (pd.api.types.is_numeric_dtype(df[first])
                       or pd.api.types.is_bool_dtype(df[first]))
    setdf = df.drop(columns=[first]) if first_is_id else df
    if setdf.shape[1] < 2:
        raise ValueError("upset needs at least 2 set columns")

    bin_df = setdf.apply(_to_member)
    set_cols = list(bin_df.columns)
    n_elements = len(bin_df)
    in_any = bin_df.any(axis=1)
    dropped = int((~in_any).sum())
    bin_df = bin_df[in_any].reset_index(drop=True)
    clean_steps = ([f"Dropped {dropped} element(s) that belonged to no set "
                    f"({n_elements} -> {len(bin_df)})."] if dropped
                   else ["No cleaning applied; every element belongs to at least one set."])

    # Distinct-mode intersection sizes: count elements per unique membership tuple.
    combo = bin_df.groupby(set_cols, sort=False).size()
    n_inter = int((combo > 0).sum())
    w = max(6.0, 0.45 * n_inter + 2.5)
    h = max(3.8, 0.45 * len(set_cols) + 2.8)

    data = upsetplot.from_indicators(set_cols, bin_df)
    fig = plt.figure(figsize=(w, h))
    upsetplot.UpSet(data, subset_size="count", sort_by="cardinality",
                    show_counts=True).plot(fig=fig)

    set_rows = [{"type": "set", "sets": c, "degree": None, "size": int(bin_df[c].sum())}
                for c in set_cols]
    inter_rows = []
    for idx, size in combo.items():
        if size == 0:
            continue
        members = [set_cols[i] for i, v in enumerate(idx if isinstance(idx, tuple) else (idx,)) if v]
        inter_rows.append({"type": "intersection", "sets": " & ".join(members),
                           "degree": len(members), "size": int(size)})
    inter_rows.sort(key=lambda r: r["size"], reverse=True)
    stats_df = pd.DataFrame(set_rows + inter_rows)

    test_meta = {
        "name": "UpSet plot (upsetplot)",
        "n_sets": len(set_cols), "n_elements": len(bin_df),
        "n_intersections": n_inter,
        "largest_intersection": int(combo.max()) if n_inter else 0,
    }
    methods = _methods(len(set_cols), len(bin_df), n_inter, set_cols)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": pd.concat([pd.Series(range(len(bin_df)), name="element").astype(str),
                              bin_df.reset_index(drop=True)], axis=1),
        "clean_steps": clean_steps,
        "label": slugify(spec["data"]["file"].rsplit(".", 1)[0]),
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(in_name, fig_stub, set_cols, first_is_id),
    }


def _methods(n_sets, n_elem, n_inter, set_cols):
    from importlib.metadata import version as v
    import platform
    return (
        f"Set intersections among {n_sets} sets ({', '.join(set_cols)}) across {n_elem} "
        f"elements were drawn as an UpSet plot with the upsetplot package. The plot shows the "
        f"{n_inter} observed intersections and each set's total size; it is descriptive, so no "
        "statistical test is applied."
        f" Rendered in Python {platform.python_version()} with upsetplot {v('upsetplot')}."
    )


def _script(in_name, fig_stub, set_cols, first_is_id):
    idcol = ", index_col=0" if first_is_id else ""
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, upsetplot

raw = pd.read_csv("{in_name}"{idcol})
set_cols = {set_cols}
bin_df = raw[set_cols].apply(lambda s: pd.to_numeric(s, errors="coerce") == 1)
bin_df = bin_df[bin_df.any(axis=1)].reset_index(drop=True)
data = upsetplot.from_indicators(set_cols, bin_df)
fig = plt.figure(figsize=(7, 4.5))
upsetplot.UpSet(data, subset_size="count", sort_by="cardinality", show_counts=True).plot(fig=fig)
fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
