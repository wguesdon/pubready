"""Recipe: UpSet plot of set intersections from a binary membership matrix.
Mirrors the R recipe (ComplexHeatmap UpSet). Input is a CSV where each row is an
element and each column is a set, with an optional leading id column. The plot
is descriptive (intersection sizes), so there is no hypothesis test; the stats
table lists set sizes and every intersection size.

The figure is drawn with matplotlib in the layout of ComplexHeatmap: the
intersection bars on top, the membership matrix below them with the set names on
the left, and the set-size bars on the right. It shows the ``top_n`` largest
intersections (15 by default) and says in its title how many it hides. The
order of the intersections and of the sets is the order of the R engine.
"""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from ..bundle import slugify

ABSENT = "#CCCCCC"   # the grey of an absent set, as in ComplexHeatmap
SHADE = "#F0F0F0"    # the background of every second set row


def _to_member(s):
    if pd.api.types.is_bool_dtype(s):
        return s.astype(bool)
    sn = pd.to_numeric(s, errors="coerce")
    if sn.notna().all() and set(sn.unique()) <= {0, 1}:
        return sn == 1
    return s.astype(str).str.lower().isin(["1", "true", "yes", "y", "t"])


def combo_table(bin_df, set_cols):
    """The observed intersections in the display order of both engines.

    The intersections are sorted by size, largest first, then by degree, then by
    their set names in code-point order, which is the C-locale order of R.

    Args:
        bin_df: The boolean membership table, one column per set.
        set_cols: The set columns.

    Returns:
        pandas.DataFrame: One row per intersection with the columns ``sets``,
        ``degree``, ``size`` and ``pattern``, the membership flags in the order
        of ``set_cols``.
    """
    combo = bin_df.groupby(set_cols, sort=False).size()
    rows = []
    for idx, size in combo.items():
        key = idx if isinstance(idx, tuple) else (idx,)
        if size == 0:
            continue
        flags = [bool(v) for v in key]
        members = [c for c, v in zip(set_cols, flags) if v]
        rows.append({"sets": " & ".join(members), "degree": len(members), "size": int(size),
                     "pattern": tuple(flags)})
    table = pd.DataFrame(rows, columns=["sets", "degree", "size", "pattern"])
    return table.sort_values(["size", "degree", "sets"], ascending=[False, True, True],
                             kind="stable").reset_index(drop=True)


def subtitle_text(n_shown, n_total):
    """The subtitle that says how many intersections the figure hides.

    Args:
        n_shown: The number of intersections drawn.
        n_total: The number of observed intersections.

    Returns:
        str or None: The text, or None when every intersection is drawn.
    """
    if n_shown >= n_total:
        return None
    return (f"The {n_shown} largest of {n_total} intersections; "
            f"{n_total - n_shown} smaller ones are not shown")


def label_inches(set_cols):
    """The width of the set-name column in inches, the estimate of the R engine.

    Args:
        set_cols: The set names.

    Returns:
        float: The width.
    """
    return 0.3 + 0.08 * max([len(str(c)) for c in set_cols] + [4])


def draw_upset(shown, set_cols, set_order, set_sizes, title):
    """Draw the UpSet plot in the layout of ComplexHeatmap.

    Args:
        shown: The intersections to draw, in order, from ``combo_table``.
        set_cols: The set columns, the order of each membership pattern.
        set_order: The set names from top to bottom.
        set_sizes: The size of each set, indexed by set name.
        title: The title text, or None.

    Returns:
        tuple: The figure, its width and its height in inches.
    """
    n, k = len(shown), len(set_order)
    lab_in = label_inches(set_order)
    w = max(6.0, lab_in + 0.32 * n + 1.9)
    h = max(3.8, 0.3 * k + 2.9)
    left = max(lab_in, 1.0)
    bottom, mat_h = 0.75, 0.3 * k
    top_pad = 0.45 if title else 0.15
    size_w = 1.0
    mat_w = w - left - size_w - 0.35
    bar_y0 = bottom + mat_h + 0.08
    fig = plt.figure(figsize=(w, h))
    ax_m = fig.add_axes([left / w, bottom / h, mat_w / w, mat_h / h])
    ax_t = fig.add_axes([left / w, bar_y0 / h, mat_w / w, (h - bar_y0 - top_pad) / h])
    ax_s = fig.add_axes([(left + mat_w + 0.1) / w, bottom / h, size_w / w, mat_h / h])

    ypos = {s: k - 1 - i for i, s in enumerate(set_order)}
    # ComplexHeatmap shades the bottom set row and every second row above it.
    for s in set_order:
        if ypos[s] % 2 == 0:
            ax_m.axhspan(ypos[s] - 0.5, ypos[s] + 0.5, color=SHADE, zorder=0, lw=0)
    xs = np.arange(n)
    ys = np.array([ypos[s] for s in set_order])
    col_of = {s: i for i, s in enumerate(set_cols)}
    for j, pattern in enumerate(shown["pattern"]):
        member = np.array([pattern[col_of[s]] for s in set_order], dtype=bool)
        ax_m.scatter(np.full(int((~member).sum()), j), ys[~member], s=60, color=ABSENT,
                     zorder=2, lw=0)
        ax_m.scatter(np.full(int(member.sum()), j), ys[member], s=60, color="black",
                     zorder=3, lw=0)
        if member.sum() > 1:
            ax_m.vlines(j, ys[member].min(), ys[member].max(), color="black", lw=1.6, zorder=2)
    ax_m.set_xlim(-0.5, n - 0.5)
    ax_m.set_ylim(-0.5, k - 0.5)
    ax_m.set_yticks([ypos[s] for s in set_order])
    ax_m.set_yticklabels(set_order, fontsize=10)
    ax_m.set_xticks([])
    ax_m.tick_params(axis="y", length=0)
    for sp in ax_m.spines.values():
        sp.set_visible(False)

    sizes = shown["size"].to_numpy()
    ax_t.bar(xs, sizes, width=0.7, color="black")
    top = float(sizes.max()) if n else 1.0
    for x_, s_ in zip(xs, sizes):
        ax_t.text(x_, s_ + 0.02 * top, f"{s_:d}", ha="center", va="bottom", fontsize=8)
    ax_t.set_xlim(-0.5, n - 0.5)
    ax_t.set_ylim(0, top * 1.15)
    ax_t.set_xticks([])
    ax_t.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10], integer=True))
    ax_t.tick_params(axis="y", labelsize=8)
    ax_t.set_ylabel("Intersection\nsize", fontsize=9, rotation=0, ha="right", va="center")
    for sp in ("top", "right", "bottom"):
        ax_t.spines[sp].set_visible(False)

    tot = np.array([set_sizes[s] for s in set_order])
    ax_s.barh([ypos[s] for s in set_order], tot, height=0.6, color="black")
    smax = float(tot.max())
    for s, v in zip(set_order, tot):
        ax_s.text(v + 0.03 * smax, ypos[s], f"{int(v):d}", ha="left", va="center", fontsize=8)
    ax_s.set_ylim(-0.5, k - 0.5)
    ax_s.set_xlim(0, smax * 1.35)
    ax_s.set_yticks([])
    ax_s.xaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    ax_s.tick_params(axis="x", labelsize=8)
    ax_s.set_xlabel("Set size", fontsize=9)
    for sp in ("top", "right", "left"):
        ax_s.spines[sp].set_visible(False)
    if title:
        fig.text((left + mat_w / 2) / w, 1 - 0.12 / h, title, ha="center", va="top",
                 fontsize=10, color="0.3")
    return fig, w, h


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
    table = combo_table(bin_df, set_cols)
    n_inter = len(table)
    top_n = max(1, int(spec["appearance"].get("top_n") or 15))
    shown = table.head(top_n)
    set_sizes = pd.Series({c: int(bin_df[c].sum()) for c in set_cols})
    set_order = list(set_sizes.sort_values(ascending=False, kind="stable").index)
    subtitle = subtitle_text(len(shown), n_inter)
    title = spec["appearance"].get("title")
    title = "\n".join(t for t in (title, subtitle) if t) or None
    fig, w, h = draw_upset(shown, set_cols, set_order, set_sizes, title)

    set_rows = [{"type": "set", "sets": c, "degree": None, "size": int(set_sizes[c])}
                for c in set_cols]
    inter_rows = [{"type": "intersection", "sets": r.sets, "degree": int(r.degree),
                   "size": int(r.size)} for r in table.itertuples(index=False)]
    stats_df = pd.DataFrame(set_rows + inter_rows)

    test_meta = {
        "name": "UpSet plot (matplotlib, ComplexHeatmap layout)",
        "n_sets": len(set_cols), "n_elements": len(bin_df),
        "n_intersections": n_inter, "n_intersections_shown": int(len(shown)),
        "largest_intersection": int(table["size"].max()) if n_inter else 0,
    }
    methods = _methods(len(set_cols), len(bin_df), n_inter, len(shown), set_cols)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": pd.concat([pd.Series(range(len(bin_df)), name="element").astype(str),
                              bin_df.reset_index(drop=True)], axis=1),
        "clean_steps": clean_steps,
        "label": slugify(spec["data"]["file"].rsplit(".", 1)[0]),
        "width": w, "height": h,
        "build_script": functools.partial(_script, set_cols, first_is_id, top_n),
    }


def _methods(n_sets, n_elem, n_inter, n_shown, set_cols):
    from importlib.metadata import version as v
    import platform
    if n_shown < n_inter:
        shown_txt = (f"The plot shows the {n_shown} largest of the {n_inter} observed "
                     "intersections and each set's total size; the stats table lists them all; ")
    else:
        shown_txt = f"The plot shows the {n_inter} observed intersections and each set's total size; "
    return (
        f"Set intersections among {n_sets} sets ({', '.join(set_cols)}) across {n_elem} "
        f"elements were drawn as an UpSet plot in the layout of ComplexHeatmap. " + shown_txt +
        "it is descriptive, so no statistical test is applied."
        f" Rendered in Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(set_cols, first_is_id, top_n, in_name, fig_stub):
    """Write the standalone script that redraws the UpSet plot.

    Args:
        set_cols: The set columns.
        first_is_id: Whether the first input column names the elements.
        top_n: The number of intersections to draw.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script text.
    """
    idcol = ", index_col=0" if first_is_id else ""
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container
# (see REPRODUCE.md) from this bundle folder.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}"{idcol})
set_cols = {set_cols!r}
bin_df = raw[set_cols].apply(pd.to_numeric, errors="coerce") == 1
bin_df = bin_df[bin_df.any(axis=1)].reset_index(drop=True)
combo = bin_df.groupby(set_cols, sort=False).size().reset_index(name="size")
combo["sets"] = [" & ".join(c for c in set_cols if r[c]) for _, r in combo.iterrows()]
combo["degree"] = combo[set_cols].sum(axis=1)
combo = combo.sort_values(["size", "degree", "sets"], ascending=[False, True, True], kind="stable")
shown = combo.head({top_n}).reset_index(drop=True)
totals = bin_df.sum().sort_values(ascending=False, kind="stable")
order = list(totals.index); k, n = len(order), len(shown)

fig = plt.figure(figsize=(max(6.0, 0.32 * n + 3.5), max(3.8, 0.3 * k + 2.9)))
gs = fig.add_gridspec(2, 2, width_ratios=[max(n, 1), 3], height_ratios=[2, max(k, 1) * 0.6],
                      wspace=0.05, hspace=0.05)
ax_t, ax_m, ax_s = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
for j in range(n):
    member = shown.loc[j, order].to_numpy(dtype=bool)
    ys = np.arange(k)[::-1]
    ax_m.scatter([j] * k, ys, s=60, color=np.where(member, "black", "#CCCCCC"), lw=0)
    if member.sum() > 1:
        ax_m.vlines(j, ys[member].min(), ys[member].max(), color="black", lw=1.6)
ax_m.set_yticks(np.arange(k)[::-1]); ax_m.set_yticklabels(order); ax_m.set_xticks([])
ax_m.set_xlim(-0.5, n - 0.5); ax_m.set_ylim(-0.5, k - 0.5)
ax_t.bar(range(n), shown["size"], width=0.7, color="black"); ax_t.set_xlim(-0.5, n - 0.5)
for j, s in enumerate(shown["size"]):
    ax_t.text(j, s, str(s), ha="center", va="bottom", fontsize=8)
ax_t.set_xticks([]); ax_t.set_ylabel("Intersection size")
ax_s.barh(np.arange(k)[::-1], totals[order], height=0.6, color="black")
for y, v in zip(np.arange(k)[::-1], totals[order]):
    ax_s.text(v, y, " %d" % v, va="center", fontsize=8)
ax_s.set_yticks([]); ax_s.set_ylim(-0.5, k - 0.5); ax_s.set_xlabel("Set size")
fig.savefig("{fig_stub}.pdf", bbox_inches="tight")
print("Reproduced figure.")
'''
