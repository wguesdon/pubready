"""Recipe: volcano plot of differential-expression results, EnhancedVolcano style.
Mirrors the R recipe. Auto-detects DESeq2 / limma / edgeR columns, four-color
scheme, cutoff lines, top hits labeled with adjustText. No test is run here; the
differential test is upstream, so the recipe is descriptive of that result."""
import functools
import os

import numpy as np
import pandas as pd

from ..bundle import slugify
from ..theme import apply_pub_style, new_fig

_FC = ["log2FoldChange", "logFC", "log2FC", "log2fc", "avg_log2FC"]
_P = ["pvalue", "P.Value", "PValue", "p_value", "pval", "p.value"]
#: The smallest positive normalised double, R's .Machine$double.xmin.
DOUBLE_XMIN = float(np.finfo(float).tiny)


def _columns(df, spec):
    def pick(cands, override):
        if override:
            return override
        for c in cands:
            if c in df.columns:
                return c
        return None
    fc = pick(_FC, spec["data"].get("x"))
    p = pick(_P, spec["data"].get("y"))
    lab = spec["data"].get("label")
    if not lab:
        obj = [c for c in df.columns if df[c].dtype == object]
        lab = obj[0] if obj else None
    return fc, p, lab


def place_labels(ax, fc, neglog, lab, top):
    """Label the top hits and move the labels clear of each other and of the points.

    adjustText repels the labels from every point of the plot, as ggrepel does in
    R, and draws a grey segment to each labelled point. A fixed iteration limit
    keeps the result the same on every run; the default time limit would not.

    Args:
        ax: The volcano axes, in its final layout.
        fc: The log2 fold changes of all features.
        neglog: The -log10 p of all features.
        lab: The feature labels.
        top: The indices of the features to label.

    Returns:
        list: The label text objects.
    """
    texts = [ax.text(fc[i], neglog[i], lab[i], fontsize=8, zorder=5) for i in top]
    if not texts:
        return texts
    try:
        from adjustText import adjust_text
    except ImportError:
        return texts
    adjust_text(texts, x=fc, y=neglog, ax=ax, target_x=fc[top], target_y=neglog[top],
                expand=(1.3, 1.6), force_text=(0.4, 0.8), force_static=(0.3, 0.6),
                force_explode=(0.3, 0.8), max_move=(20, 20), iter_lim=600, min_arrow_len=0,
                arrowprops=dict(arrowstyle="-", color="0.35", lw=0.5))
    return texts


def recipe_volcano(df, spec):
    fc_col, p_col, lab_col = _columns(df, spec)
    if fc_col is None:
        raise ValueError("volcano: no log2 fold-change column found; set it with --x")
    if p_col is None:
        raise ValueError("volcano: no p-value column found; set it with --y")

    fc = pd.to_numeric(df[fc_col], errors="coerce").to_numpy()
    p = pd.to_numeric(df[p_col], errors="coerce").to_numpy()
    lab = (df[lab_col].astype(str).to_numpy() if lab_col
           else np.array([str(i) for i in range(len(df))]))
    n0 = len(df)
    bad = int(np.sum(np.isfinite(p) & ((p < 0) | (p > 1))))
    if bad:
        raise ValueError(f"volcano: p must lie in 0 to 1; column '{p_col}' holds {bad} "
                         "value(s) outside 0 to 1")
    mask = np.isfinite(fc) & np.isfinite(p)
    fc, p, lab = fc[mask], p[mask], lab[mask]
    if len(fc) < 1:
        raise ValueError("volcano: no usable rows after cleaning")
    # A p of exactly 0 is the strongest hit, not a missing value. It is set to the
    # smallest positive double, as manhattan does, and drawn as a triangle.
    is_zero = p == 0
    n_zero = int(is_zero.sum())
    p = np.where(is_zero, DOUBLE_XMIN, p)
    dropped = int((~mask).sum())
    clean_steps = []
    if dropped:
        clean_steps.append(f"Dropped {dropped} row(s) with missing or non-finite fold change or p "
                           f"({n0} -> {len(fc)}).")
    if n_zero:
        clean_steps.append(f"Replaced {n_zero} p value(s) of exactly 0 with the smallest positive "
                           f"double ({DOUBLE_XMIN:.3g}), so that -log10 p is finite.")
    if not clean_steps:
        clean_steps = ["No cleaning applied; input used as-is."]
    df_used = pd.DataFrame({(lab_col or "label"): lab, fc_col: fc, p_col: p})

    ap = spec["appearance"]
    fcc = ap.get("fc_cutoff", 1.0)
    pcc = ap.get("p_cutoff", 0.05)
    topn = int(ap.get("top_n", 15))

    neglog = -np.log10(p)
    passfc = np.abs(fc) >= fcc
    passp = p <= pcc
    cat = np.where(passfc & passp, "FC and p", np.where(passp, "p", np.where(passfc, "FC", "NS")))
    colors = {"NS": "#B3B3B3", "FC": "#2E8B57", "p": "#3B6DB3", "FC and p": "#C1432B"}
    # The legend and the y axis name the p column the threshold reads, as in R.
    plab = "p" if p_col in _P else p_col
    legend = {"NS": "NS", "FC": "FC", "p": plab, "FC and p": f"FC and {plab}"}

    ydisp = display_y(neglog, is_zero, pcc)
    w, h = 6.0, 5.6
    fig, ax = new_fig(w, h)
    for c in ("NS", "FC", "p", "FC and p"):
        m = cat == c
        if m.any():
            ax.scatter(fc[m & ~is_zero], ydisp[m & ~is_zero], s=10, alpha=0.7, c=colors[c],
                       label=legend[c], edgecolors="none")
            if (m & is_zero).any():
                ax.scatter(fc[m & is_zero], ydisp[m & is_zero], s=14, alpha=0.7, c=colors[c],
                           marker="^", edgecolors="none")
    ax.axvline(-fcc, ls="--", lw=0.6, color="grey")
    ax.axvline(fcc, ls="--", lw=0.6, color="grey")
    ax.axhline(-np.log10(pcc), ls="--", lw=0.6, color="grey")

    sig_idx = np.where(cat == "FC and p")[0]
    if len(sig_idx):
        order = sig_idx[np.argsort(-(np.abs(fc[sig_idx]) * neglog[sig_idx]))]
    else:
        order = np.argsort(-neglog)
    top = order[:topn]

    ax.set_xlabel(ap.get("x_label") or "log2 fold change")
    ax.set_ylabel(ap.get("y_label") or f"-log10 {plab}")
    if n_zero:
        ax.set_title(f"{n_zero} feature(s) with p = 0 drawn as triangles at the top",
                     fontsize=10, color="#4d4d4d")
        if ap.get("title"):
            fig.suptitle(ap["title"], fontweight="bold")
    elif ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False,
                    fontsize=9, markerscale=2.0)
    for handle in leg.legend_handles:
        handle.set_alpha(1.0)
    fig.tight_layout()
    # The labels are placed last, once the layout is final: adjustText works in
    # display units, so a later change of the axes size would undo its spacing.
    place_labels(ax, fc, ydisp, lab, top)

    n_up = int(np.sum((fc >= fcc) & (p <= pcc)))
    n_down = int(np.sum((fc <= -fcc) & (p <= pcc)))
    stats_df = pd.DataFrame({
        "gene": lab[top], "log2FC": fc[top], "p_value": p[top],
        "neg_log10_p": neglog[top], "direction": np.where(fc[top] > 0, "up", "down"),
    })
    test_meta = {"name": "Volcano plot", "fc_cutoff": fcc, "p_cutoff": pcc,
                 "n_total": len(fc), "n_up": n_up, "n_down": n_down,
                 "n_significant": n_up + n_down, "n_zero_p": n_zero}
    methods = _methods(len(fc), fcc, pcc, n_up, n_down, plab, len(top), n_zero)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": clean_steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": functools.partial(_script, fc_col, p_col, lab_col, fcc, pcc, topn),
    }


def display_y(neglog, is_zero, pcc):
    """-log10 p for drawing; a p of exactly 0 sits at a ceiling.

    The ceiling is 12% above the largest other value, or above the p cutoff
    line, as in the R engine and in manhattan.

    Args:
        neglog: The -log10 p of every feature, zeros already set to the
            smallest positive double.
        is_zero: True where the input p was exactly 0.
        pcc: The p cutoff.

    Returns:
        numpy.ndarray: The values to draw.
    """
    y = np.array(neglog, dtype=float)
    top = float(y[~is_zero].max()) if (~is_zero).any() else 0.0
    y[is_zero] = max(top, -np.log10(pcc)) * 1.12
    return y


def _methods(n, fcc, pcc, n_up, n_down, plab, n_lab, n_zero=0):
    from importlib.metadata import version as v
    import platform
    if n_up + n_down > 0:
        lab_txt = (f"The {n_lab} significant features with the largest product of "
                   f"|log2 fold change| and -log10 {plab} are labelled. ")
    else:
        lab_txt = (f"No feature was significant, so the {n_lab} features with the smallest "
                   f"{plab} are labelled. ")
    zero_txt = (f"{n_zero} {plab} value(s) of exactly 0 were set to the smallest positive double "
                f"({DOUBLE_XMIN:.3g}) and are drawn as triangles at the top of the plot. "
                if n_zero else "")
    return (
        f"Differential-expression results for {n} features were shown as a volcano plot. "
        f"Features with |log2 fold change| >= {float(fcc):g} and {plab} <= {float(pcc):g} were called "
        f"significant; {n_up} were up-regulated and {n_down} down-regulated. {zero_txt}{lab_txt}"
        "The differential test itself "
        "was performed upstream; this figure summarises its output. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')} and adjustText {v('adjustText')}."
    )


def _script(fc_col, p_col, lab_col, fcc, pcc, topn, in_name, fig_stub):
    lab_line = (f'lab = raw["{lab_col}"].astype(str).to_numpy()' if lab_col
                else "lab = raw.index.astype(str).to_numpy()")
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
fc = pd.to_numeric(raw["{fc_col}"], errors="coerce").to_numpy()
p = pd.to_numeric(raw["{p_col}"], errors="coerce").to_numpy()
{lab_line}
m = np.isfinite(fc) & np.isfinite(p)
fc, p, lab = fc[m], np.where(p[m] == 0, np.finfo(float).tiny, p[m]), lab[m]
neglog = -np.log10(p)
fcc, pcc = {fcc}, {pcc}
cat = np.where((np.abs(fc) >= fcc) & (p <= pcc), "FC and p",
      np.where(p <= pcc, "p", np.where(np.abs(fc) >= fcc, "FC", "NS")))
colors = {{"NS": "grey", "FC": "#2E8B57", "p": "#3B6DB3", "FC and p": "#C1432B"}}
fig, ax = plt.subplots(figsize=(6, 5.6))
for c in ("NS", "FC", "p", "FC and p"):
    mm = cat == c
    if mm.any():
        ax.scatter(fc[mm], neglog[mm], s=10, alpha=0.7, c=colors[c], label=c, edgecolors="none")
ax.axvline(-fcc, ls="--", lw=0.6, color="grey"); ax.axvline(fcc, ls="--", lw=0.6, color="grey")
ax.axhline(-np.log10(pcc), ls="--", lw=0.6, color="grey")
ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 p")
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, frameon=False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
