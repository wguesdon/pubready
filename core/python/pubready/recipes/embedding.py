"""Recipe: UMAP or t-SNE scatter coloured by cluster.

Mirrors the R recipe. The coordinates come from upstream (Seurat, scanpy);
pubready does not compute the embedding, because the two engines cannot give
the same UMAP. --x and --y are the two coordinates and --group the cluster.
Each cluster is labelled at its median position, and the stats give the size,
the median position and the mean silhouette width of each cluster on the 2D
coordinates, with the rules of R's cluster::silhouette. A numeric --group with
more than 12 values (and not a set of at most 50 whole-number cluster ids) is a
continuous score, coloured on the viridis scale without cluster labels.
"""
import functools
import os

import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from ..bundle import slugify
from ..theme import R_POINT_PT, apply_pub_style, new_fig, palette

MAX_LEGEND = 12
MAX_INT_IDS = 50
MAX_SIL = 5000
MAX_POINTS = 50000
#: Rows of the distance matrix computed at a time.
CHUNK = 1000


def is_continuous(v):
    """Decide whether a group column is a continuous score.

    Args:
        v: The group values, as a pandas Series without missing values.

    Returns:
        bool: True for a numeric column with more than MAX_LEGEND values that
        is not a set of at most MAX_INT_IDS whole numbers.
    """
    if not pd.api.types.is_numeric_dtype(v):
        return False
    u = np.unique(v.to_numpy(dtype=float))
    if len(u) <= MAX_LEGEND:
        return False
    return not (np.all(u == np.round(u)) and len(u) <= MAX_INT_IDS)


def group_label(v):
    """Text of one cluster value, as R's as.character writes it.

    Args:
        v: The value.

    Returns:
        str: The label; a whole number has no decimal part.
    """
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        return str(int(v))
    return str(v)


def cluster_levels(v):
    """Cluster levels: numbers in numeric order, other names in code-point order.

    Args:
        v: The group values, as a pandas Series without missing values.

    Returns:
        list: The level labels.
    """
    if pd.api.types.is_numeric_dtype(v):
        return [group_label(u) for u in sorted(pd.unique(v))]
    u = list(pd.unique(v.astype(str)))
    num = pd.to_numeric(pd.Series(u), errors="coerce")
    if num.notna().all():
        return [u[i] for i in np.argsort(num.to_numpy(), kind="stable")]
    return sorted(u)


def silhouette_widths(xy, codes, k):
    """Silhouette width of each point, with the rules of cluster::silhouette.

    a(i) is the mean distance to the other members of the own cluster, b(i)
    the smallest mean distance to another cluster, and s(i) = (b - a) /
    max(a, b). A point alone in its cluster, or with a(i) = b(i), gets 0.

    Args:
        xy: The coordinates, n by 2.
        codes: The cluster code of each point, 0 to k - 1.
        k: The count of clusters.

    Returns:
        numpy.ndarray: The widths, or None when fewer than 2 clusters are
        present or every point is its own cluster.
    """
    n = len(codes)
    present = np.unique(codes)
    if len(present) <= 1 or len(present) >= n:
        return None
    onehot = np.zeros((n, k))
    onehot[np.arange(n), codes] = 1.0
    counts = onehot.sum(0)
    sums = np.empty((n, k))
    for i0 in range(0, n, CHUNK):
        dx = xy[i0:i0 + CHUNK, 0][:, None] - xy[None, :, 0]
        dy = xy[i0:i0 + CHUNK, 1][:, None] - xy[None, :, 1]
        sums[i0:i0 + CHUNK] = np.sqrt(dx * dx + dy * dy) @ onehot
    own = counts[codes]
    a = np.where(own > 1, sums[np.arange(n), codes] / np.maximum(own - 1, 1), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_other = sums / counts[None, :]
    mean_other[np.arange(n), codes] = np.inf
    mean_other[:, counts == 0] = np.inf
    b = mean_other.min(1)
    s = np.zeros(n)
    ok = (own > 1) & (a != b)
    s[ok] = (b[ok] - a[ok]) / np.maximum(a[ok], b[ok])
    return s


def silhouette_summary(xy, codes, k):
    """Mean silhouette width per cluster and overall, on a subsample when large.

    Above MAX_SIL cells, every k-th cell in input order is used, with
    k = ceil(n / MAX_SIL), as the R engine does.

    Args:
        xy: The coordinates, n by 2.
        codes: The cluster code of each cell.
        k: The count of clusters.

    Returns:
        dict: The per-cluster means, the overall mean, the sample size and the step.
    """
    n = len(codes)
    step = int(np.ceil(n / MAX_SIL)) if n > MAX_SIL else 1
    idx = np.arange(0, n, step)
    out = {"per_cluster": np.full(k, np.nan), "overall": np.nan, "n": int(len(idx)), "step": step}
    s = silhouette_widths(xy[idx], codes[idx], k)
    if s is None:
        return out
    for c in range(k):
        m = codes[idx] == c
        if m.any():
            out["per_cluster"][c] = float(s[m].mean())
    out["overall"] = float(s.mean())
    return out


def fig_dims(x, y, legend, continuous, lv):
    """The figure size from the aspect ratio of the data, as the R engine sets it.

    Args:
        x: The first coordinates.
        y: The second coordinates.
        legend: Whether a legend is drawn.
        continuous: Whether the colour is a continuous score.
        lv: The cluster names.

    Returns:
        tuple: The width and the height in inches.
    """
    r = min(2.0, max(0.5, float(np.ptp(y)) / max(float(np.ptp(x)), 1e-12)))
    panel_w = min(4.4, 5.6 / r)
    if not legend:
        leg_w = 0.0
    elif continuous:
        leg_w = 1.2
    else:
        leg_w = 0.8 + 0.085 * max(len(s) for s in lv)
    return panel_w + 0.7 + leg_w, panel_w * r + 1.1


def point_size(n):
    """The point size in R units (mm), shrinking with n, as the R engine sets it.

    Args:
        n: The count of drawn points.

    Returns:
        float: The size.
    """
    return min(1.6, max(0.3, 1.6 * np.sqrt(500 / n)))


def point_alpha(n):
    """The point opacity, shrinking with n, as the R engine sets it.

    Args:
        n: The count of drawn points.

    Returns:
        float: The alpha.
    """
    return 0.9 if n <= 1000 else max(0.5, 0.9 * (1000 / n) ** 0.2)


def fmt_count(n):
    """Format a count with thousands separators.

    Args:
        n: The count.

    Returns:
        str: The formatted count.
    """
    return f"{int(n):,}"


def recipe_embedding(df, spec):
    """Draw the embedding coloured by cluster with labels at the medians.

    Args:
        df: One row per cell.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing.
    """
    dat = spec["data"]
    xc, yc, gc = dat.get("x"), dat.get("y"), dat.get("group")
    if not xc or not yc:
        raise ValueError("embedding needs --x and --y (the two embedding coordinates, for example "
                         "UMAP_1 and UMAP_2)")
    if not gc:
        raise ValueError("embedding needs --group (the cluster column)")
    for c in (xc, yc, gc):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    x = pd.to_numeric(df[xc], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[yc], errors="coerce").to_numpy(dtype=float)
    g = df[gc]
    keep = np.isfinite(x) & np.isfinite(y) & g.notna().to_numpy()
    if pd.api.types.is_numeric_dtype(g):
        keep &= np.isfinite(g.to_numpy(dtype=float))
    x, y, g = x[keep], y[keep], g[keep].reset_index(drop=True)
    n = len(x)
    if n < 2:
        raise ValueError("embedding needs at least 2 cells with both coordinates and a group")
    steps = ([f"Dropped {int((~keep).sum())} row(s) with a missing coordinate or group "
              f"({n0} -> {n})."] if int((~keep).sum()) else ["No cleaning applied; input used as-is."])

    continuous = is_continuous(g)
    step_draw = int(np.ceil(n / MAX_POINTS)) if n > MAX_POINTS else 1
    draw = np.arange(0, n, step_draw)
    size, alpha = point_size(len(draw)), point_alpha(len(draw))
    s_pt = (size * R_POINT_PT) ** 2  # scatter takes the marker area in points^2
    ap = spec["appearance"]
    x_lab, y_lab = ap.get("x_label") or xc, ap.get("y_label") or yc

    if continuous:
        stats_df = pd.DataFrame({"cluster": ["all"], "n": [n], "percent": [100.0],
                                 "median_x": [float(np.median(x))], "median_y": [float(np.median(y))],
                                 "mean_silhouette": [np.nan]})
        sil = {"overall": np.nan, "n": 0, "step": 1}
        lv = []
        legend = True
        sub = f"n = {fmt_count(n)} cells; colour = {gc}"
    else:
        lv = cluster_levels(g)
        gl = np.array([group_label(v) for v in g])
        index = {c: i for i, c in enumerate(lv)}
        codes = np.array([index[c] for c in gl])
        counts = np.bincount(codes, minlength=len(lv))
        mx = np.array([np.median(x[codes == c]) for c in range(len(lv))])
        my = np.array([np.median(y[codes == c]) for c in range(len(lv))])
        if len(lv) >= 2:
            sil = silhouette_summary(np.column_stack([x, y]), codes, len(lv))
        else:
            sil = {"per_cluster": np.array([np.nan]), "overall": np.nan, "n": 0, "step": 1}
        stats_df = pd.DataFrame({"cluster": lv, "n": counts, "percent": 100 * counts / n,
                                 "median_x": mx, "median_y": my,
                                 "mean_silhouette": sil["per_cluster"]})
        legend = len(lv) <= MAX_LEGEND
        sil_txt = "NA" if np.isnan(sil["overall"]) else f"{sil['overall']:.2f}"
        sub = f"n = {fmt_count(n)} cells; {len(lv)} clusters; mean silhouette = {sil_txt}"

    w, h = fig_dims(x, y, legend, continuous, lv)
    fig, ax = new_fig(w, h)
    if continuous:
        sc = ax.scatter(x[draw], y[draw], c=g.to_numpy(dtype=float)[draw], cmap="viridis", s=s_pt,
                        alpha=alpha, linewidths=0)
        # The colour bar shows the scale without the point transparency.
        cb = fig.colorbar(ScalarMappable(norm=sc.norm, cmap=sc.cmap), ax=ax, shrink=0.45,
                          aspect=12, pad=0.03)
        cb.ax.set_title(gc, fontsize=11, loc="left", pad=8)
        cb.outline.set_visible(False)
        cb.locator = MaxNLocator(nbins=4)
        cb.update_ticks()
        cb.ax.tick_params(labelsize=10)
    else:
        pal = palette(len(lv), ap.get("palette"))
        ax.scatter(x[draw], y[draw], c=[pal[c] for c in codes[draw]], s=s_pt, alpha=alpha,
                   linewidths=0)
        texts = [ax.text(mx[c], my[c], lv[c], ha="center", va="center",
                         fontsize=7.5 if len(lv) > 20 else 9, zorder=5,
                         bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none",
                               "alpha": 0.9})
                 for c in range(len(lv))]
        if len(texts) > 1:
            try:
                from adjustText import adjust_text
                # A fixed iteration count: without one, adjustText stops after one
                # second, so the labels would depend on the speed of the machine.
                adjust_text(texts, ax=ax, expand=(1.05, 1.15), force_text=(0.2, 0.3),
                            ensure_inside_axes=True, iter_lim=2000)
            except Exception:
                pass
        if legend:
            # The R legend overrides the key size to 3.
            handles = [Line2D([0], [0], marker="o", ls="", markersize=3 * R_POINT_PT,
                              markeredgewidth=0, markerfacecolor=pal[c], label=lv[c])
                       for c in range(len(lv))]
            ax.legend(handles=handles, title=gc, loc="center left", bbox_to_anchor=(1.02, 0.5),
                      frameon=False, fontsize=10, title_fontsize=11, handletextpad=0.3,
                      alignment="left")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_title(sub, fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    fig.tight_layout()

    test_meta = {"name": "Embedding scatter", "n_cells": n, "n_clusters": len(lv),
                 "continuous": bool(continuous), "mean_silhouette": sil["overall"],
                 "silhouette_n": sil["n"], "silhouette_step": sil["step"],
                 "n_drawn": int(len(draw)), "drawing_step": step_draw}
    df_used = pd.DataFrame({xc: x, yc: y, gc: g})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(test_meta, xc, yc, gc), "df_used": df_used, "clean_steps": steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]), "width": w, "height": h,
        "build_script": functools.partial(_script, xc, yc, gc, continuous),
    }


def _methods(tm, xc, yc, gc):
    """Write the methods paragraph, as the R engine words it.

    Args:
        tm: The test metadata.
        xc: The first coordinate column.
        yc: The second coordinate column.
        gc: The group column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    head = (f"The embedding coordinates {xc} and {yc} of {fmt_count(tm['n_cells'])} cells were "
            "computed upstream and are drawn as given; pubready does not compute the embedding. ")
    draw = ""
    if tm["drawing_step"] > 1:
        draw = (f"For the size of the vector outputs, one in every {tm['drawing_step']} cells in "
                f"input order is drawn ({fmt_count(tm['n_drawn'])} of {fmt_count(tm['n_cells'])}); "
                "the counts and the medians use all cells. ")
    if tm["continuous"]:
        return (head + f"Points are coloured by the continuous value {gc} on the viridis scale; no "
                f"cluster statistic applies. {draw}Drawn in Python {platform.python_version()} with "
                f"matplotlib {v('matplotlib')}.")
    sub = ""
    if tm["silhouette_step"] > 1:
        sub = (f"Because there were {fmt_count(tm['n_cells'])} cells, the silhouette was computed on "
               f"one in every {tm['silhouette_step']} cells in input order "
               f"({fmt_count(tm['silhouette_n'])} cells). ")
    overall = ("not defined" if np.isnan(tm["mean_silhouette"])
               else f"{tm['mean_silhouette']:.3f}")
    return (
        head + f"Points are coloured by {gc} ({tm['n_clusters']} clusters), and each cluster is "
        "labelled at its median position. The mean silhouette width of each cluster was computed "
        "from the Euclidean distances between cells on the two embedding coordinates (the rules of "
        "cluster::silhouette, restated with numpy). It describes how well the clusters separate in "
        "the embedding, not in the original feature space, because an embedding such as UMAP or "
        f"t-SNE distorts distances. {sub}The overall mean silhouette width was {overall}. {draw}"
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')} and "
        f"adjustText {v('adjustText')}."
    )


def _script(xc, yc, gc, continuous, in_name, fig_stub):
    """Write the standalone script that redraws the embedding.

    Args:
        xc: The first coordinate column.
        yc: The second coordinate column.
        gc: The group column.
        continuous: Whether the group is a continuous score.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    colour = (f'c=raw["{gc}"], cmap="viridis"' if continuous
              else f'c=pd.factorize(raw["{gc}"].astype(str))[0], cmap="tab20"')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

raw = pd.read_csv("{in_name}").dropna(subset=["{xc}", "{yc}", "{gc}"])
fig, ax = plt.subplots(figsize=(6.6, 5.2))
ax.scatter(raw["{xc}"], raw["{yc}"], s={(0.6 * R_POINT_PT) ** 2:.2f}, linewidths=0, {colour})
if not {continuous}:
    for name, sub in raw.groupby(raw["{gc}"].astype(str)):
        ax.text(sub["{xc}"].median(), sub["{yc}"].median(), name, ha="center",
                bbox={{"boxstyle": "round", "fc": "white", "ec": "none"}})
ax.set_aspect("equal", adjustable="datalim"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("{xc}"); ax.set_ylabel("{yc}")
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
