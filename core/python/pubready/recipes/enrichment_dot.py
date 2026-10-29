"""Recipe: dot plot of GSEA / ORA enrichment results, clusterProfiler style.
Mirrors the R recipe. An NES column -> GSEA dot plot (x = NES); otherwise an ORA
dot plot (x = gene ratio). Color = -log10 adjusted p, size = gene count. Descriptive; no
test is run here, the enrichment is computed upstream."""
import functools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..bundle import slugify
from ..theme import apply_pub_style, significance_cmap

_TERM = ["Description", "ID", "pathway", "term", "Term", "NAME"]
_PADJ = ["p.adjust", "padj", "qvalue", "q_value", "FDR", "pvalue", "p_value"]
_COUNT = ["Count", "setSize", "count", "size"]
_RATIO = ["GeneRatio", "gene_ratio", "generatio"]
_NES = ["NES", "nes", "enrichmentScore", "ES"]


def _pick(df, cands, override=None):
    if override:
        return override
    for c in cands:
        if c in df.columns:
            return c
    return None


def _parse_ratio(series):
    out = []
    for s in series.astype(str):
        try:
            if "/" in s:
                a, b = s.split("/")[:2]
                out.append(float(a) / float(b))
            else:
                out.append(float(s))
        except Exception:
            out.append(np.nan)
    return np.array(out)


def greedy_lines(words, width):
    """Fill words greedily into lines of at most ``width`` characters.

    Args:
        words: The words of one label.
        width: The longest line in characters; a longer word keeps its own line.

    Returns:
        list: The lines.
    """
    lines = [words[0]]
    for wd in words[1:]:
        if len(lines[-1]) + 1 + len(wd) <= width:
            lines[-1] = lines[-1] + " " + wd
        else:
            lines.append(wd)
    return lines


def wrap_label(s, width=45):
    """Wrap one term label at spaces into balanced lines of at most ``width`` characters.

    The label takes the fewest lines that the width allows, each about the same
    length, so no short word is left alone on a line. The R engine wraps the
    same way, so both engines break a long GO term at the same words.

    Args:
        s: The term label.
        width: The longest line in characters.

    Returns:
        str: The label with a newline between lines.
    """
    words = str(s).split()
    joined = " ".join(words)
    if len(words) < 2 or len(joined) <= width:
        return joined
    n_lines = len(greedy_lines(words, width))
    for w2 in range(int(np.ceil(len(joined) / n_lines)), width + 1):
        lines = greedy_lines(words, w2)
        if len(lines) <= n_lines:
            return "\n".join(lines)
    return "\n".join(greedy_lines(words, width))


def figure_size(labels):
    """The figure size from the wrapped labels, the same sizes as the R engine.

    Each term row grows with its lines, and the label column with its longest
    line.

    Args:
        labels: The wrapped labels.

    Returns:
        tuple: The width, the height and the label column width, in inches.
    """
    parts = [lab.split("\n") for lab in labels]
    longest = max([len(p) for lines in parts for p in lines] + [4])
    lab_in = 0.3 + 0.075 * longest
    w = max(7.2, lab_in + 3.2 + 1.6)
    h = max(4.0, 1.6 + sum(0.2 + 0.16 * len(lines) for lines in parts))
    return w, h, lab_in


def size_breaks(lo, hi):
    """Round count values for the size legend, as ggplot2 picks them.

    Args:
        lo: The smallest count.
        hi: The largest count.

    Returns:
        list: The values inside [lo, hi], at least one.
    """
    from matplotlib.ticker import MaxNLocator
    if hi <= lo:
        return [lo]
    ticks = MaxNLocator(nbins=6, steps=[1, 2, 5, 10]).tick_values(lo, hi)
    inside = [float(v) for v in ticks if lo <= v <= hi]
    return inside or [lo, hi]


def dot_area(ct, lo, rng):
    """The marker area of a count, 30 to 230 square points.

    Args:
        ct: The count, or an array of counts.
        lo: The smallest count.
        rng: The count range, max minus min.

    Returns:
        The area in square points.
    """
    arr = np.asarray(ct, dtype=float)
    if not rng:
        return np.full_like(arr, 100.0)
    return 30 + 200 * (arr - lo) / rng


def recipe_enrichment_dot(df, spec):
    term_col = _pick(df, _TERM, spec["data"].get("label"))
    padj_col = _pick(df, _PADJ)
    if term_col is None:
        raise ValueError("enrichment_dot: no term/Description column found; set it with --label")
    if padj_col is None:
        raise ValueError("enrichment_dot: no adjusted-p column found")
    count_col = _pick(df, _COUNT)
    ratio_col = _pick(df, _RATIO)
    nes_col = _pick(df, _NES)

    term = df[term_col].astype(str).to_numpy()
    padj = pd.to_numeric(df[padj_col], errors="coerce").to_numpy()
    count = (pd.to_numeric(df[count_col], errors="coerce").to_numpy()
             if count_col else np.full(len(term), np.nan))
    gsea = nes_col is not None
    if gsea:
        xval = pd.to_numeric(df[nes_col], errors="coerce").to_numpy(); xlab = "NES"
    elif ratio_col:
        xval = _parse_ratio(df[ratio_col]); xlab = "gene ratio"
    else:
        xval = -np.log10(padj); xlab = "-log10 adjusted p"

    keep = np.isfinite(padj) & np.isfinite(xval)
    term, padj, count, xval = term[keep], padj[keep], count[keep], xval[keep]
    if len(term) < 1:
        raise ValueError("enrichment_dot: no usable rows")
    df_used = pd.DataFrame({"term": term, "x": xval, "p_adjust": padj, "count": count})

    ap = spec["appearance"]
    topn = int(ap.get("top_n", 15))
    # Stable sorts break ties in the row order, as the R engine does.
    idx = np.argsort(padj, kind="stable")[:topn]
    # The stats table keeps the p.adjust order of the R engine; the plot runs
    # from the smallest x at the bottom to the largest at the top.
    top = {"term": term[idx], "x": xval[idx], "p_adjust": padj[idx], "count": count[idx]}
    yo = np.argsort(top["x"], kind="stable")
    t, x, pj, ct = top["term"][yo], top["x"][yo], top["p_adjust"][yo], top["count"][yo]
    has_size = not np.all(np.isnan(ct))

    labels = [wrap_label(s) for s in t]
    w, h, lab_in = figure_size(labels)
    top_in = 0.6 if ap.get("title") else 0.3
    bottom_in = 0.75
    leg_in = 1.6
    panel_w = w - lab_in - leg_in - 0.2
    panel_h = h - top_in - bottom_in
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([lab_in / w, bottom_in / h, panel_w / w, panel_h / h])
    ypos = np.arange(len(t))
    if has_size:
        lo_ct = float(np.nanmin(ct))
        rng = float(np.nanmax(ct)) - lo_ct
        sizes = dot_area(ct, lo_ct, rng)
    else:
        sizes = 70
    # Colour by -log10 of the adjusted p: the values span orders of magnitude, so a
    # linear scale would paint almost every term the same colour.
    sc = ax.scatter(x, ypos, c=-np.log10(pj), s=sizes, cmap=significance_cmap(), edgecolors="none")
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.6, len(t) - 0.4)
    ax.set_xlabel(ap.get("x_label") or xlab)
    if gsea:
        ax.axvline(0, ls="--", lw=0.6, color="grey")
    if ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    ax.tick_params(axis="y", labelsize=10)

    # The legends sit in their own column right of the panel, as in R: the count
    # sizes on top and the colour bar below them.
    leg_x = (lab_in + panel_w + 0.25) / w
    panel_top = bottom_in + panel_h
    cb_top = panel_top
    if has_size:
        brk = size_breaks(lo_ct, lo_ct + rng)
        handles = [ax.scatter([], [], s=float(dot_area(b, lo_ct, rng)), color="black",
                              edgecolors="none") for b in brk]
        ax.legend(handles, [f"{b:g}" for b in brk], title="count", loc="upper left",
                  bbox_to_anchor=(leg_x, panel_top / h), bbox_transform=fig.transFigure,
                  frameon=False, fontsize=9, title_fontsize=10, borderaxespad=0,
                  alignment="left", labelspacing=0.9)
        cb_top = panel_top - (0.4 + 0.24 * len(brk))
    cb_h = max(0.8, min(1.4, cb_top - bottom_in - 0.3))
    cax = fig.add_axes([leg_x + 0.05 / w, (cb_top - 0.3 - cb_h) / h, 0.18 / w, cb_h / h])
    fig.colorbar(sc, cax=cax)
    cax.set_title("-log10 p.adjust", fontsize=10, loc="left", pad=8)
    cax.tick_params(labelsize=9)

    stats_df = pd.DataFrame({
        "term": top["term"],
        ("NES" if gsea else "gene_ratio" if ratio_col else "neg_log10_padj"): top["x"],
        "p_adjust": top["p_adjust"], "count": top["count"]})
    test_meta = {"name": "GSEA dot plot" if gsea else "ORA dot plot",
                 "n_terms_total": len(term), "n_terms_shown": len(t), "x_axis": xlab}
    methods = _methods(gsea, len(term), len(t), xlab)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": ["No cleaning applied; enrichment table used as-is."],
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": functools.partial(_script, term_col, padj_col, count_col, ratio_col,
                                          nes_col, topn, w, h),
    }


def _methods(gsea, n_total, n_shown, xlab):
    from importlib.metadata import version as v
    import platform
    kind = "GSEA" if gsea else "over-representation (ORA)"
    return (
        f"The top {n_shown} of {n_total} enriched terms from a {kind} analysis were shown as a "
        f"dot plot. The x-axis is {xlab}, dot color the -log10 of the adjusted p-value, and "
        "dot size the gene count. The enrichment analysis was performed upstream; this figure "
        "summarises its output. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(term_col, padj_col, count_col, ratio_col, nes_col, topn, w, h, in_name, fig_stub):
    """Write the standalone script that redraws the dot plot.

    Args:
        term_col: The term column.
        padj_col: The adjusted-p column.
        count_col: The count column, or None.
        ratio_col: The gene-ratio column, or None.
        nes_col: The NES column, or None.
        topn: The number of terms to show.
        w: The figure width in inches.
        h: The figure height in inches.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script text.
    """
    if nes_col:
        xline = f'x = pd.to_numeric(raw["{nes_col}"], errors="coerce").to_numpy(); xlab = "NES"'
    elif ratio_col:
        xline = (f'x = np.array([float(s.split("/")[0]) / float(s.split("/")[1]) if "/" in s else float(s) '
                 f'for s in raw["{ratio_col}"].astype(str)]); xlab = "gene ratio"')
    else:
        xline = f'x = -np.log10(pd.to_numeric(raw["{padj_col}"], errors="coerce").to_numpy()); xlab = "-log10 adjusted p"'
    count_line = (f'ct = pd.to_numeric(raw["{count_col}"], errors="coerce").to_numpy()'
                  if count_col else "ct = np.full(len(raw), 60.0)")
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
term = raw["{term_col}"].astype(str).to_numpy()
padj = pd.to_numeric(raw["{padj_col}"], errors="coerce").to_numpy()
{count_line}
{xline}
idx = np.argsort(padj, kind="stable")[:{topn}]
t, x, pj, ct = term[idx], x[idx], padj[idx], ct[idx]
yo = np.argsort(x, kind="stable"); t, x, pj, ct = t[yo], x[yo], pj[yo], ct[yo]
# Wrap each term at about 45 characters, as the figure does.
t = [textwrap.fill(s, 45, break_long_words=False, break_on_hyphens=False) for s in t]
fig, ax = plt.subplots(figsize=({w:.2f}, {h:.2f}))
sc = ax.scatter(x, np.arange(len(t)), c=-np.log10(pj), s=70, cmap="coolwarm", edgecolors="none")
ax.set_yticks(np.arange(len(t))); ax.set_yticklabels(t)
ax.set_xlabel(xlab)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.colorbar(sc, ax=ax, label="-log10 p.adjust", fraction=0.046, pad=0.04)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
