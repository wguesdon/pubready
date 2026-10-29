"""Recipe: MA plot of differential-expression results.

Mirrors the R recipe. M is the log2 fold change and A the mean expression of
each gene. DESeq2 (log2FoldChange, padj, baseMean), limma (logFC, adj.P.Val,
AveExpr) and edgeR (logFC, FDR, logCPM) names are detected; --x sets the mean
column and --y the p column. DESeq2's baseMean is drawn on the log10 scale.
Genes that pass both cutoffs are red (up) or blue (down); the rest are grey and
drawn first. No test is run here: the differential test is upstream.
"""
import functools
import os

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig

FC_COLS = ["log2FoldChange", "logFC", "log2FC", "log2fc", "avg_log2FC"]
PADJ_COLS = ["padj", "adj.P.Val", "FDR", "p_val_adj", "qvalue"]
RAWP_COLS = ["pvalue", "P.Value", "PValue", "p_value", "pval", "p.value"]
MEAN_COLS = ["baseMean", "AveExpr", "logCPM"]
MEAN_TITLES = {"baseMean": "log10 mean of normalized counts",
               "AveExpr": "Average log2 expression", "logCPM": "log2 CPM"}

#: Above this many grey points, the grey points that share a grid cell are drawn once.
THIN_ABOVE = 20000
#: The grid of the thinning, in cells per axis.
GRID_BINS = 600
#: The colours of the three categories, as in the R engine (grey70 for NS).
COLOURS = {"NS": "#B3B3B3", "Down": "#3B6DB3", "Up": "#C1432B"}
#: A gene label longer than this is cut and ends in an ellipsis on the figure.
LABEL_MAX = 18


def short_label(s, max_chars=LABEL_MAX):
    """Cut a long label for the figure; the stats table keeps the full name.

    Args:
        s: The label.
        max_chars: The longest label drawn in full.

    Returns:
        str: The label, or its first characters and an ellipsis.
    """
    s = str(s)
    return s if len(s) <= max_chars else s[:max_chars - 1] + "…"


def label_points(ax, xs, ys, labels, fontsize):
    """Draw repelled labels next to points, as ggrepel does in the R engine.

    Args:
        ax: The axes.
        xs: The x positions of the labelled points.
        ys: The y positions of the labelled points.
        labels: The label texts.
        fontsize: The font size in points.
    """
    texts = [ax.text(x, y, lab, fontsize=fontsize, zorder=5) for x, y, lab in zip(xs, ys, labels)]
    if not texts:
        return
    try:
        from adjustText import adjust_text
        adjust_text(texts, x=np.asarray(xs), y=np.asarray(ys), ax=ax,
                    expand=(1.15, 1.35), force_text=(0.4, 0.6), iter_lim=2000,
                    arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.5})
    except Exception:
        pass


def first_hit(cands, names):
    """Return the first candidate column that the table holds.

    Args:
        cands: The candidate names, in order of preference.
        names: The column names of the table.

    Returns:
        str or None: The first match.
    """
    for c in cands:
        if c in names:
            return c
    return None


def ma_columns(df, spec):
    """Resolve the fold-change, p, mean and label columns.

    Args:
        df: The results table.
        spec: The figure specification.

    Returns:
        dict: The column names; a missing one is None.
    """
    names = list(df.columns)
    dat = spec["data"]
    p_col = dat.get("y") or first_hit(PADJ_COLS + RAWP_COLS, names)
    a_col = dat.get("x") or first_hit(MEAN_COLS, names)
    lab = dat.get("label")
    if not lab:
        txt = [c for c in names if df[c].dtype == object]
        lab = txt[0] if txt else None
    return {"fc": first_hit(FC_COLS, names), "p": p_col, "a": a_col, "lab": lab}


def grid_first(x, y, nbin):
    """Flag the first point in each cell of an nbin by nbin grid over the range.

    Args:
        x: The x values.
        y: The y values.
        nbin: The cells per axis.

    Returns:
        numpy.ndarray: A boolean mask, True for the points to draw.
    """
    span_x = max(float(np.max(x) - np.min(x)), 1e-12)
    span_y = max(float(np.max(y) - np.min(y)), 1e-12)
    bx = np.floor((x - np.min(x)) / span_x * nbin)
    by = np.floor((y - np.min(y)) / span_y * nbin)
    _, first = np.unique(bx * (nbin + 1) + by, return_index=True)
    keep = np.zeros(len(x), dtype=bool)
    keep[first] = True
    return keep


def recipe_ma_plot(df, spec):
    """Draw M against A with the significant genes coloured and labelled.

    Args:
        df: The differential-expression results table.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing or a p value lies outside 0 to 1.
    """
    cols = ma_columns(df, spec)
    if cols["fc"] is None:
        raise ValueError("ma_plot: no log2 fold-change column found; name it one of "
                         + ", ".join(FC_COLS))
    if cols["p"] is None:
        raise ValueError("ma_plot: no p-value column found; set it with --y")
    if cols["a"] is None:
        raise ValueError("ma_plot: no mean-expression column (baseMean, AveExpr, logCPM) "
                         "found; set it with --x")
    for c in [cols["p"], cols["a"]] + ([cols["lab"]] if cols["lab"] else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    fc = pd.to_numeric(df[cols["fc"]], errors="coerce").to_numpy(dtype=float)
    p = pd.to_numeric(df[cols["p"]], errors="coerce").to_numpy(dtype=float)
    a_raw = pd.to_numeric(df[cols["a"]], errors="coerce").to_numpy(dtype=float)
    log_a = cols["a"] == "baseMean"
    if log_a:
        with np.errstate(divide="ignore", invalid="ignore"):
            a = np.where(np.isfinite(a_raw) & (a_raw > 0), np.log10(np.maximum(a_raw, 1e-300)), np.nan)
    else:
        a = a_raw
    lab = (df[cols["lab"]].astype(str).to_numpy() if cols["lab"]
           else np.array([str(i + 1) for i in range(len(df))]))
    if np.any((p < 0) | (p > 1)):
        raise ValueError(f"ma_plot: column '{cols['p']}' holds values outside 0 to 1")

    n0 = len(df)
    keep = np.isfinite(fc) & np.isfinite(a)
    fc, p, a, lab = fc[keep], p[keep], a[keep], lab[keep]
    if len(fc) < 1:
        raise ValueError("ma_plot: no usable rows after cleaning")
    n_missing_p = int(np.isnan(p).sum())
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing or non-finite fold change "
                     f"or mean{', or a baseMean of zero' if log_a else ''} ({n0} -> {len(fc)}).")
    if n_missing_p:
        steps.append(f"Kept {n_missing_p} row(s) with a missing p as not significant.")
    if log_a:
        steps.append("Plotted log10(baseMean) as the mean expression A.")
    if not steps:
        steps = ["No cleaning applied; input used as-is."]

    ap = spec["appearance"]
    fcc = ap.get("fc_cutoff", 1.0)
    pcc = ap.get("p_cutoff", 0.05)
    topn = int(ap.get("top_n", 15))

    with np.errstate(invalid="ignore"):
        sig = ~np.isnan(p) & (p <= pcc)
    up = sig & (fc >= fcc)
    down = sig & (fc <= -fcc)
    cat = np.where(up, "Up", np.where(down, "Down", "NS"))
    n_up, n_down = int(up.sum()), int(down.sum())
    n_ns = len(fc) - n_up - n_down

    hit = np.where(up | down)[0]
    hit = hit[np.lexsort((hit, -np.abs(fc[hit]), p[hit]))]
    stats_df = pd.DataFrame({"gene": lab[hit], "A": a[hit], "log2FC": fc[hit], "p_value": p[hit],
                             "direction": np.where(fc[hit] > 0, "up", "down")})
    top = hit[:topn]

    ns_idx = np.where(cat == "NS")[0]
    thinned = n_ns > THIN_ABOVE
    if thinned:
        ns_idx = ns_idx[grid_first(a[ns_idx], fc[ns_idx], GRID_BINS)]

    plab = "p" if cols["p"] in RAWP_COLS else cols["p"]
    a_title = MEAN_TITLES.get(cols["a"], cols["a"])
    # The y range covers the data and both cutoff lines, with a 6% margin.
    ylo, yhi = min(float(fc.min()), -1.5 * fcc), max(float(fc.max()), 1.5 * fcc)
    ylim = (ylo - 0.06 * (yhi - ylo), yhi + 0.06 * (yhi - ylo))

    w, h = 6.0, 5.4
    fig, ax = new_fig(w, h)
    ax.scatter(a[ns_idx], fc[ns_idx], s=1.6, c=COLOURS["NS"], linewidths=0, zorder=1)
    ax.axhline(0, color="#262626", lw=0.9, zorder=2)
    for yv in (-fcc, fcc):
        ax.axhline(yv, color="#666666", lw=0.7, ls="--", zorder=2)
    order = hit[::-1]
    if len(order):
        ax.scatter(a[order], fc[order], s=9, alpha=0.85, linewidths=0, zorder=3,
                   c=[COLOURS[c] for c in cat[order]])
    ax.set_ylim(*ylim)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]))
    label_points(ax, a[top], fc[top], [short_label(s) for s in lab[top]],
                 7.5 if len(top) > 10 else 8.5)

    ax.set_xlabel(ap.get("x_label") or a_title)
    ax.set_ylabel(ap.get("y_label") or "log2 fold change")
    ax.set_title(f"{n_up} up / {n_down} down ({plab} <= {pcc:g}, |log2FC| >= {fcc:g})",
                 fontsize=10, color="#4d4d4d")
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    handles = [Line2D([0], [0], marker="o", ls="", markersize=6, markeredgewidth=0,
                      markerfacecolor=COLOURS[k], label=lbl)
               for k, lbl in (("NS", "Not significant"), ("Down", "Down"), ("Up", "Up"))]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
              frameon=False, fontsize=10, handletextpad=0.2, columnspacing=1.2)
    fig.tight_layout()

    test_meta = {"name": "MA plot", "fc_cutoff": fcc, "p_cutoff": pcc, "p_column": cols["p"],
                 "mean_column": cols["a"], "mean_log10": log_a, "n_total": int(len(fc)),
                 "n_up": n_up, "n_down": n_down, "n_not_significant": int(n_ns),
                 "n_missing_p": n_missing_p, "n_grey_drawn": int(len(ns_idx))}
    df_used = pd.DataFrame({(cols["lab"] or "label"): lab, "A": a, cols["fc"]: fc, cols["p"]: p})
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta,
        "methods": _methods(test_meta, plab, a_title, len(top), thinned),
        "df_used": df_used, "clean_steps": steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": functools.partial(_script, cols, log_a, fcc, pcc, topn, a_title),
    }


def _methods(tm, plab, a_title, n_labels, thinned):
    """Write the methods paragraph, as the R engine words it.

    Args:
        tm: The test metadata.
        plab: The name of the p column in the text.
        a_title: The description of the mean expression.
        n_labels: The count of labelled genes.
        thinned: Whether the grey points were thinned.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    txt = (
        f"Differential-expression results for {tm['n_total']} genes were shown as an MA plot of "
        f"the log2 fold change (M) against the mean expression (A, {a_title}). "
        f"Genes with {plab} <= {tm['p_cutoff']:g} and |log2 fold change| >= {tm['fc_cutoff']:g} "
        f"were called significant: {tm['n_up']} were up-regulated, {tm['n_down']} down-regulated "
        f"and {tm['n_not_significant']} not significant. "
    )
    if tm["n_missing_p"]:
        txt += f"{tm['n_missing_p']} genes without a p value were counted as not significant. "
    if n_labels:
        txt += f"The {n_labels} significant genes with the smallest p are labelled. "
    if thinned:
        txt += (f"Grey points that fall in the same cell of a {GRID_BINS} by {GRID_BINS} grid were "
                f"drawn once ({tm['n_grey_drawn']} of {tm['n_not_significant']} drawn); every "
                "count uses all genes. ")
    return txt + (
        "The differential test itself was performed upstream; this figure summarises its output. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')} and "
        f"adjustText {v('adjustText')}."
    )


def _script(cols, log_a, fcc, pcc, topn, a_title, in_name, fig_stub):
    """Write the standalone script that redraws the MA plot.

    Args:
        cols: The resolved column names.
        log_a: Whether A is log10 of baseMean.
        fcc: The fold-change cutoff.
        pcc: The p cutoff.
        topn: The count of labels.
        a_title: The x axis title.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    lab_line = (f'lab = raw["{cols["lab"]}"].astype(str).to_numpy()' if cols["lab"]
                else "lab = (raw.index + 1).astype(str).to_numpy()")
    a_line = (f'a = np.log10(pd.to_numeric(raw["{cols["a"]}"], errors="coerce").to_numpy(float))'
              if log_a else f'a = pd.to_numeric(raw["{cols["a"]}"], errors="coerce").to_numpy(float)')
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
m = pd.to_numeric(raw["{cols["fc"]}"], errors="coerce").to_numpy(float)
p = pd.to_numeric(raw["{cols["p"]}"], errors="coerce").to_numpy(float)
{a_line}
{lab_line}
keep = np.isfinite(m) & np.isfinite(a)
m, p, a, lab = m[keep], p[keep], a[keep], lab[keep]
fcc, pcc = {fcc}, {pcc}
sig = ~np.isnan(p) & (p <= pcc)
up, down = sig & (m >= fcc), sig & (m <= -fcc)
hit = np.where(up | down)[0]
hit = hit[np.lexsort((hit, -np.abs(m[hit]), p[hit]))]
fig, ax = plt.subplots(figsize=(6, 5.4))
ns = ~(up | down)
ax.scatter(a[ns], m[ns], s=1.6, c="#B3B3B3", linewidths=0)
ax.axhline(0, color="black", lw=0.9)
ax.axhline(fcc, color="grey", ls="--", lw=0.7); ax.axhline(-fcc, color="grey", ls="--", lw=0.7)
ax.scatter(a[up], m[up], s=9, c="#C1432B", linewidths=0, label="Up")
ax.scatter(a[down], m[down], s=9, c="#3B6DB3", linewidths=0, label="Down")
for i in hit[:{topn}]:
    ax.annotate(lab[i], (a[i], m[i]), fontsize=8)
lo, hi = min(m.min(), -1.5 * fcc), max(m.max(), 1.5 * fcc)
ax.set_ylim(lo - 0.06 * (hi - lo), hi + 0.06 * (hi - lo))
ax.set_xlabel("{a_title}"); ax.set_ylabel("log2 fold change")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
