"""Recipe: clustered heatmap of a features-by-samples matrix. Tries
PyComplexHeatmap (to mirror the R ComplexHeatmap output); falls back to a seaborn
clustermap if PyComplexHeatmap raises. Row z-score and clustering by default;
a per-feature Welch t-test (BH) is computed with --test welch_t when the
annotation has two groups, and its stars are drawn between the heatmap and the
row names, as in R. The default, --test auto, runs no test.

Both engines draw the same leaf order (ordered_linkage): a complete-linkage tree
on the Euclidean distance of the drawn values, and at each merge the branch with
the higher mean of the unscaled values first. On a tie, the branch that holds the
earlier input row or column comes first."""
import functools
import inspect
import os

import matplotlib.pyplot as plt
from matplotlib.legend import Legend
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from statsmodels.stats.multitest import multipletests

from ..bundle import slugify
from ..io import read_tidy
from ..theme import PALETTE, diverging_cmap
from ..util import level_labels, level_order


def ordered_linkage(values, weights):
    """Build the tree of one heatmap axis with the leaf order of the R engine.

    The tree is complete linkage on the Euclidean distance of the rows of values.
    At each merge the branch with the higher mean weight comes first, and on a
    tie the branch that holds the earlier input index. scipy draws the first
    column of a linkage row on the left, so the swap sets the drawn order, as
    .hm_order_tree() does on the merge matrix in R.

    Args:
        values: The drawn values, one row per leaf.
        weights: The weight of each leaf, the mean of its unscaled values.

    Returns:
        numpy.ndarray: The linkage matrix.
    """
    z = hierarchy.linkage(np.asarray(values, dtype=float), method="complete", metric="euclidean")
    n = len(weights)
    sum_w = np.concatenate([np.asarray(weights, dtype=float), np.zeros(n - 1)])
    size = np.concatenate([np.ones(n), np.zeros(n - 1)])
    first = np.concatenate([np.arange(n, dtype=float), np.zeros(n - 1)])
    for i in range(n - 1):
        a, b = int(z[i, 0]), int(z[i, 1])
        ma, mb = sum_w[a] / size[a], sum_w[b] / size[b]
        if mb > ma or (mb == ma and first[b] < first[a]):
            z[i, 0], z[i, 1] = b, a
        sum_w[n + i] = sum_w[a] + sum_w[b]
        size[n + i] = size[a] + size[b]
        first[n + i] = min(first[a], first[b])
    return z


def color_limits(disp, scale_mode):
    """Colour limits that match the R engine.

    Args:
        disp: The displayed matrix.
        scale_mode: row, column or none.

    Returns:
        tuple: The lower limit, the centre and the upper limit.
    """
    if scale_mode == "none":
        q = np.quantile(disp.to_numpy(), [0.02, 0.5, 0.98])
        return float(q[0]), float(q[1]), float(q[2])
    return -2.0, 0.0, 2.0


def level_colours(series):
    """House-palette colours for the sorted levels of one annotation column.

    Args:
        series: The annotation column.

    Returns:
        dict: One colour for each level, as the R engine assigns them.
    """
    return {lvl: PALETTE[i % len(PALETTE)] for i, lvl in enumerate(sorted(series.unique()))}


def heatmap_stars(padj):
    """Significance stars of the BH-adjusted P values, as the R engine draws them.

    Args:
        padj: The adjusted P values, one per feature.

    Returns:
        list: "***" below 0.001, "**" below 0.01, "*" below 0.05, else "".
    """
    out = []
    for p in padj:
        if p is None or np.isnan(p):
            out.append("")
        elif p < 0.001:
            out.append("***")
        elif p < 0.01:
            out.append("**")
        elif p < 0.05:
            out.append("*")
        else:
            out.append("")
    return out


def welch_requested(method, group_col, n_groups):
    """Whether to run the per-feature Welch test, as .hm_welch_requested() in R.

    --test auto, the default, draws no test; --test welch_t runs it and needs an
    annotation whose first column holds two groups. Any other value stops the run.

    Args:
        method: The --test value.
        group_col: The first annotation column, or None without an annotation.
        n_groups: The count of distinct values in that column.

    Returns:
        bool: True for --test welch_t.

    Raises:
        ValueError: If the value is unknown or the annotation does not hold two
            groups.
    """
    method = method or "auto"
    if method not in ("auto", "welch_t"):
        raise ValueError(f"heatmap: unknown --test '{method}'; use auto (no test) or welch_t")
    if method == "auto":
        return False
    if group_col is None:
        raise ValueError("heatmap: --test welch_t needs --annotation, a table whose first column "
                         "after the sample id holds two groups")
    if n_groups != 2:
        raise ValueError(f"heatmap: --test welch_t compares two groups, but the annotation column "
                         f"'{group_col}' holds {n_groups}")
    return True


def welch_p(row, ca, cb):
    """The Welch t-test P value of one feature between two sample sets.

    The P is NaN when a group has fewer than two values or when the standard
    error of the difference is zero at double precision, the rule of R's
    t.test ("data are essentially constant"), as .hm_welch_p() in R (F87).

    Args:
        row: The values of one feature, indexed by sample.
        ca: The samples of the first group.
        cb: The samples of the second group.

    Returns:
        float: The two-sided P value, or NaN.
    """
    x = np.asarray(row[ca], dtype=float)
    y = np.asarray(row[cb], dtype=float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    se = np.sqrt(np.var(x, ddof=1) / len(x) + np.var(y, ddof=1) / len(y))
    if not np.isfinite(se) or se <= 10 * np.finfo(float).eps * max(abs(x.mean()), abs(y.mean())):
        return float("nan")
    return float(stats.ttest_ind(x, y, equal_var=False).pvalue)


def welch_groups(g, group_col):
    """The two groups of the Welch test, as .hm_test_groups() in R.

    Args:
        g: The group of each sample, indexed by sample, missing kept as NaN.
        group_col: The annotation column, for the message.

    Returns:
        list: The two groups in the order of level_order().

    Raises:
        ValueError: If a sample has no group or a group has fewer than two samples.
    """
    empty = [str(s) for s, v in g.items() if pd.isna(v) or not str(v).strip()]
    if empty:
        raise ValueError(f"heatmap: --test welch_t needs a group for each sample, but the "
                         f"annotation column '{group_col}' is empty for "
                         + ", ".join(f"'{s}'" for s in empty))
    lab = level_labels(g)
    lv = level_order(lab)
    small = [(lvl, int((lab == lvl).sum())) for lvl in lv if int((lab == lvl).sum()) < 2]
    if small:
        raise ValueError("heatmap: --test welch_t needs at least 2 samples in each group; "
                         + ", ".join(f"group '{lvl}' has {n}" for lvl, n in small))
    return lv


def bh_tested(p):
    """Benjamini-Hochberg adjustment over the P values that are not NaN.

    R's p.adjust() leaves a missing P out of the count, so a feature that was
    not tested does not change the adjusted P of the others.

    Args:
        p: The P values, NaN for a feature that was not tested.

    Returns:
        numpy.ndarray: The adjusted P values, NaN where p is NaN.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


# The gap in points between the heatmap and the row names when stars sit there.
STAR_PAD = 20.0


def add_row_stars(ax, star_of):
    """Write the significance stars between the heatmap and its row names.

    The row names move right by ``STAR_PAD`` points to make room, as the
    ComplexHeatmap row annotation of the R engine does.

    Args:
        ax: The heatmap axes that carries the row names on its right.
        star_of: The star text of each feature, keyed by feature name.
    """
    ax.yaxis.set_tick_params(pad=STAR_PAD)
    for y_, lab in zip(ax.get_yticks(), ax.get_yticklabels()):
        star = star_of.get(lab.get_text(), "")
        if star:
            ax.annotate(star, xy=(1.0, y_), xycoords=("axes fraction", "data"),
                        xytext=(3, 0), textcoords="offset points", ha="left",
                        va="center", fontsize=9, annotation_clip=False)


def label_overhang(fig, axes):
    """How far the row names reach past the right edge of the heatmap.

    Args:
        fig: The figure, drawn here so that the tick labels have their size.
        axes: The heatmap axes that carry the row names.

    Returns:
        float: The overhang as a fraction of the figure width, at least 0.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    right = max(a.get_window_extent(renderer).x1 for a in axes)
    ends = [t.get_window_extent(renderer).x1 for a in axes for t in a.get_yticklabels()
            if t.get_text() and t.get_visible()]
    ends += [c.get_window_extent(renderer).x1 for a in axes for c in a.texts]
    reach = max(ends) - right if ends else 0.0
    return max(reach, 0.0) / fig.get_window_extent(renderer).width


def unclip_legends(fig):
    """Let each legend draw past the axes that PyComplexHeatmap made for it.

    PyComplexHeatmap sizes the legend axes from its own width estimate, and
    add_artist clips each legend to that axes. A long level name then runs past
    the estimate, the legend is cut, and the tight bounding box of the saved
    figure stops at the same edge.

    Args:
        fig: The figure.
    """
    for ax in fig.axes:
        for child in ax.get_children():
            if isinstance(child, Legend):
                child.set_clip_on(False)


def row_label_axes(plotter):
    """The heatmap axes of a PyComplexHeatmap plotter that carry the row names.

    Args:
        plotter: The ClusterMapPlotter.

    Returns:
        list: The right-most column of heatmap axes.
    """
    axes = getattr(plotter, "heatmap_axes", None)
    if axes is None:
        return [plotter.ax_heatmap]
    return list(np.ravel(np.asarray(axes)[:, -1]))


def recipe_heatmap(df, spec):
    feat = df.iloc[:, 0].astype(str)
    mat = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    mat.index = list(feat)
    n0 = len(mat)
    keep = mat.notna().all(axis=1) & (mat.std(axis=1) > 0)
    mat = mat[keep]
    clean_steps = (f"Dropped {n0 - len(mat)} feature(s) with missing values or zero variance "
                   f"({n0} -> {len(mat)})." if len(mat) < n0 else "No cleaning applied; input used as-is.")
    df_used = df[df.iloc[:, 0].astype(str).isin(mat.index)].reset_index(drop=True)

    ap = spec["appearance"]
    scale_mode = ap.get("scale") or "row"
    cluster = ap.get("cluster") or "both"
    if scale_mode == "row":
        disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0)
    elif scale_mode == "column":
        disp = (mat - mat.mean()) / mat.std()
    else:
        disp = mat.copy()
    disp = disp.fillna(0.0)
    legend_name = "value" if scale_mode == "none" else "z-score"
    cluster_rows = cluster in ("both", "rows")
    cluster_cols = cluster in ("both", "columns")
    row_z = (ordered_linkage(disp.to_numpy(), mat.mean(axis=1).to_numpy())
             if cluster_rows and disp.shape[0] > 1 else None)
    col_z = (ordered_linkage(disp.to_numpy().T, mat.mean(axis=0).to_numpy())
             if cluster_cols and disp.shape[1] > 1 else None)
    row_order = (list(disp.index[hierarchy.leaves_list(row_z)]) if row_z is not None
                 else list(disp.index))
    col_order = (list(disp.columns[hierarchy.leaves_list(col_z)]) if col_z is not None
                 else list(disp.columns))

    ann, group_col, n_groups = None, None, 0
    if spec["data"].get("annotation"):
        araw = read_tidy(spec["data"]["annotation"])
        ids = level_labels(araw.iloc[:, 0])
        miss = [c for c in mat.columns if c not in set(ids)]
        if miss:
            raise ValueError(f"heatmap: the annotation {os.path.basename(spec['data']['annotation'])} "
                             "has no row for the sample(s) "
                             + ", ".join(f"'{m}'" for m in miss) + " of the matrix")
        araw = araw.iloc[:, 1:].set_axis(list(ids), axis=0)
        ann = araw[~araw.index.duplicated()].reindex(mat.columns)
        group_col = ann.columns[0]
        n_groups = int(ann[group_col].dropna().nunique())
        group_raw = ann[group_col]
        ann = ann.astype(str)

    # Per-feature differential test, only with --test welch_t. On genes that a DE
    # test chose, a second test on the same samples is circular, so the default
    # draws no test (F97).
    run_test = welch_requested(spec["test"].get("method"), group_col, n_groups)
    star_of = {}
    n_untested = 0
    if run_test:
        g = level_labels(group_raw)
        lv = welch_groups(g, group_col)
        ca = [c for c in mat.columns if g.loc[c] == lv[0]]
        cb = [c for c in mat.columns if g.loc[c] == lv[1]]
        # A feature with no variance in either group has no Welch test. Its P is
        # NaN, and the BH adjustment runs over the tested features only (F87).
        pvals = mat.apply(welch_p, axis=1, args=(ca, cb)).to_numpy(dtype=float)
        n_untested = int(np.isnan(pvals).sum())
        padj = bh_tested(pvals)
        star_of = dict(zip(mat.index, heatmap_stars(padj)))
        stats_df = pd.DataFrame({
            "feature": mat.index, f"mean_{lv[0]}": mat[ca].mean(axis=1).to_numpy(),
            f"mean_{lv[1]}": mat[cb].mean(axis=1).to_numpy(),
            "p_value": pvals, "p_adj": padj,
            "significant": np.where(np.isnan(padj), None, padj < 0.05),
        })
        test_label = f"per-feature Welch t-test ({lv[0]} vs {lv[1]}), BH-adjusted"
        test_levels = (lv[0], lv[1])
        n_sig = int(np.sum(padj < 0.05))
    else:
        stats_df = pd.DataFrame({"feature": mat.index, "mean": mat.mean(axis=1).to_numpy(),
                                 "sd": mat.std(axis=1).to_numpy()})
        test_label, test_levels, n_sig = None, None, None

    w = max(6.5, 3.2 + 0.24 * mat.shape[1])
    h = max(5.0, 1.5 + 0.16 * mat.shape[0])

    # The colour scale matches the R engine: z-scores from -2 to 2 with white at
    # zero, or the 2nd, 50th and 98th percentiles of unscaled values. Clustering
    # is complete linkage on Euclidean distance, the ComplexHeatmap default.
    vmin, vcenter, vmax = color_limits(disp, scale_mode)
    used_pch, fallback_reason = True, None
    try:
        import PyComplexHeatmap as pch
        col_ha = None
        if ann is not None:
            col_ha = pch.HeatmapAnnotation(
                **{c: pch.anno_simple(ann[c], colors=level_colours(ann[c]), add_text=False,
                                      legend=True)
                   for c in ann.columns},
                axis=1, verbose=0)
        plt.figure(figsize=(w, h))
        # The legends are drawn after the row names are measured. PyComplexHeatmap
        # measures them before the first draw, when they are still empty, and so
        # puts the legends on top of them.
        plotter = pch.ClusterMapPlotter(
            data=disp, top_annotation=col_ha,
            row_cluster=row_z is not None, col_cluster=col_z is not None,
            row_cluster_method="complete", col_cluster_method="complete",
            row_dendrogram=row_z is not None, col_dendrogram=col_z is not None,
            row_dendrogram_kws={"linkage": row_z}, col_dendrogram_kws={"linkage": col_z},
            cmap=diverging_cmap(), vmin=vmin, vmax=vmax, center=vcenter, label=legend_name,
            show_rownames=True, show_colnames=True,
            xticklabels=True, yticklabels=True,
            xticklabels_kws={"labelsize": 8}, yticklabels_kws={"labelsize": 8},
            row_names_side="right", col_names_side="bottom", plot_legend=False, verbose=0)
        fig = plt.gcf()
        label_axes = row_label_axes(plotter)
        if any(star_of.values()):
            for a in label_axes:
                add_row_stars(a, star_of)
        gap = 4.0 / 25.4 * fig.dpi / fig.get_window_extent().width
        plotter.legend_delta_x = label_overhang(fig, label_axes) + gap
        plotter.plot_legends(ax=plotter.ax_heatmap)
        unclip_legends(fig)
    except Exception as err:  # noqa: BLE001 - any PyComplexHeatmap failure falls back
        plt.close("all")
        used_pch = False
        fallback_reason = f"{type(err).__name__}: {err}"
        col_colors = None
        if ann is not None:
            col_colors = pd.DataFrame({c: ann[c].map(level_colours(ann[c])) for c in ann.columns})
        grid = sns.clustermap(
            disp, cmap=diverging_cmap(), vmin=vmin, vmax=vmax, center=vcenter, method="complete",
            row_cluster=row_z is not None, col_cluster=col_z is not None,
            row_linkage=row_z, col_linkage=col_z,
            col_colors=col_colors, figsize=(w, h), cbar_kws={"label": legend_name},
            xticklabels=True, yticklabels=True)
        fig = grid.figure
        if any(star_of.values()):
            add_row_stars(grid.ax_heatmap, star_of)

    test_meta = {
        "name": "Clustered heatmap", "engine_lib": "PyComplexHeatmap" if used_pch else "seaborn.clustermap",
        "fallback_reason": fallback_reason,
        "scale": scale_mode, "clustering": cluster,
        "n_features": int(mat.shape[0]), "n_samples": int(mat.shape[1]),
        "per_feature_test": test_label, "n_significant": n_sig, "n_untested": n_untested,
        "row_order": [str(v) for v in row_order], "column_order": [str(v) for v in col_order],
    }
    methods = _methods(scale_mode, cluster, mat.shape[0], mat.shape[1], test_levels, n_sig, used_pch,
                       n_untested)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": clean_steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "extra_inputs": [spec["data"].get("annotation")],
        "width": w, "height": h,
        "build_script": functools.partial(_script, spec, scale_mode, cluster_rows, cluster_cols,
                                          legend_name),
    }


def _methods(scale_mode, cluster, n_feat, n_samp, test_levels, n_sig, used_pch, n_untested=0):
    from importlib.metadata import version as v
    import platform
    lib = "PyComplexHeatmap" if used_pch else "seaborn.clustermap"
    lib_ver = v("PyComplexHeatmap") if used_pch else v("seaborn")
    scale_txt = {"row": "Values were z-scored across samples within each feature (row scaling).",
                 "column": "Values were z-scored within each sample (column scaling).",
                 "none": "Raw values were shown without scaling."}.get(scale_mode)
    clust_txt = {"both": "Rows and columns were hierarchically clustered.",
                 "rows": "Rows were hierarchically clustered; column order was preserved.",
                 "columns": "Columns were hierarchically clustered; row order was preserved.",
                 }.get(cluster, "No clustering was applied.")
    if cluster in ("both", "rows", "columns"):
        clust_txt += (" Clustering used complete linkage on the Euclidean distance of the displayed "
                      "values, and at each merge the branch with the higher mean of the unscaled "
                      "values was drawn first.")
    test_txt = " No statistical test was applied."
    if test_levels:
        n_tested = n_feat - n_untested
        untested_txt = (f" {n_untested} feature(s) had no variance in either group, so they were "
                        "not tested and were left out of the adjustment." if n_untested else "")
        sig_txt = (f"{n_sig} of the {n_tested} features had an adjusted P below 0.05 and are "
                   "marked with stars." if n_sig else
                   "no feature had an adjusted P below 0.05, so no star is drawn.")
        test_txt = (f" Each feature was compared between {test_levels[0]} and {test_levels[1]} with a "
                    "Welch two-sample t-test, which treats the samples as independent, and the P "
                    f"values were adjusted across the {n_tested} tested features with the "
                    f"Benjamini-Hochberg method; {sig_txt}{untested_txt}")
    return (
        f"A heatmap of {n_feat} features across {n_samp} samples was drawn with {lib}. "
        + scale_txt + " " + clust_txt + test_txt +
        f" Rendered in Python {platform.python_version()} with {lib} {lib_ver}."
    )


def _script(spec, scale_mode, cluster_rows, cluster_cols, legend_name, in_name, fig_stub):
    ann = spec["data"].get("annotation")
    ann_name = os.path.basename(ann) if ann else None
    disp_line = {"row": "disp = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1), axis=0).fillna(0)",
                 "column": "disp = ((mat - mat.mean()) / mat.std()).fillna(0)",
                 "none": "disp = mat"}.get(scale_mode)
    ann_lines = (f'ann = pd.read_csv("{ann_name}", index_col=0).astype(str).reindex(mat.columns)\n'
                 f'gc = ann.columns[0]\n'
                 f'colors = {{lvl: c for lvl, c in zip(sorted(ann[gc].unique()), ["#3B6DB3", "#C1432B", "#2E8B57"])}}\n'
                 f'col_colors = ann[gc].map(colors)') if ann_name else "col_colors = None"
    row_line = ("row_z = ordered_linkage(disp.to_numpy(), mat.mean(axis=1).to_numpy()) "
                "if len(disp) > 1 else None" if cluster_rows else "row_z = None")
    col_line = ("col_z = ordered_linkage(disp.to_numpy().T, mat.mean(axis=0).to_numpy()) "
                "if disp.shape[1] > 1 else None" if cluster_cols else "col_z = None")
    return f'''#!/usr/bin/env python3
# Standalone reproduction (seaborn clustermap). Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import numpy as np
import pandas as pd, seaborn as sns
from scipy.cluster import hierarchy


{inspect.getsource(ordered_linkage)}

raw = pd.read_csv("{in_name}", index_col=0)
mat = raw.apply(pd.to_numeric, errors="coerce")
mat = mat[mat.notna().all(axis=1) & (mat.std(axis=1) > 0)]
{disp_line}
{ann_lines}
{row_line}
{col_line}
g = sns.clustermap(disp, cmap="RdBu_r", center=0, row_cluster=row_z is not None,
                   col_cluster=col_z is not None, row_linkage=row_z, col_linkage=col_z,
                   col_colors=col_colors, cbar_kws={{"label": "{legend_name}"}},
                   xticklabels=True, yticklabels=True)
g.figure.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
