"""Recipe: dot plot of GSEA / ORA enrichment results, clusterProfiler style.
Mirrors the R recipe. An NES column -> GSEA dot plot (x = NES); otherwise an ORA
dot plot (x = gene ratio). Color = adjusted p, size = gene count. Descriptive; no
test is run here, the enrichment is computed upstream."""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..bundle import slugify
from ..theme import apply_pub_style

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
    idx = np.argsort(padj)[:topn]
    t, x, pj, ct = term[idx], xval[idx], padj[idx], count[idx]
    yo = np.argsort(x)
    t, x, pj, ct = t[yo], x[yo], pj[yo], ct[yo]
    has_size = not np.all(np.isnan(ct))

    w = 7.2
    h = max(4.0, 0.34 * len(t) + 1.6)
    fig, ax = plt.subplots(figsize=(w, h))
    ypos = np.arange(len(t))
    if has_size:
        rng = np.nanmax(ct) - np.nanmin(ct)
        sizes = 30 + 200 * ((ct - np.nanmin(ct)) / rng if rng else np.zeros_like(ct))
    else:
        sizes = 70
    sc = ax.scatter(x, ypos, c=pj, s=sizes, cmap="coolwarm_r", edgecolors="none")
    ax.set_yticks(ypos); ax.set_yticklabels(t, fontsize=9)
    ax.set_xlabel(ap.get("x_label") or xlab)
    if gsea:
        ax.axvline(0, ls="--", lw=0.6, color="grey")
    if ap.get("title"):
        ax.set_title(ap["title"], fontweight="bold")
    apply_pub_style(ax, ap.get("theme"))
    fig.colorbar(sc, ax=ax, label="p.adjust", fraction=0.046, pad=0.04)
    if has_size:
        handles, labels = sc.legend_elements(prop="sizes", num=3, alpha=0.6,
                                             func=lambda s: (s - 30) / 200 * rng + np.nanmin(ct))
        ax.legend(handles, labels, title="count", loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()

    stats_df = pd.DataFrame({
        "term": t, ("NES" if gsea else "gene_ratio" if ratio_col else "neg_log10_padj"): x,
        "p_adjust": pj, "count": ct})
    test_meta = {"name": "GSEA dot plot" if gsea else "ORA dot plot",
                 "n_terms_total": len(term), "n_terms_shown": len(t), "x_axis": xlab}
    methods = _methods(gsea, len(term), len(t), xlab)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df_used, "clean_steps": ["No cleaning applied; enrichment table used as-is."],
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(in_name, fig_stub, term_col, padj_col, count_col, ratio_col, nes_col, topn),
    }


def _methods(gsea, n_total, n_shown, xlab):
    from importlib.metadata import version as v
    import platform
    kind = "GSEA" if gsea else "over-representation (ORA)"
    return (
        f"The top {n_shown} of {n_total} enriched terms from a {kind} analysis were shown as a "
        f"dot plot. The x-axis is {xlab}, dot color the adjusted p-value, and dot size the gene "
        "count. The enrichment analysis was performed upstream; this figure summarises its output. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(in_name, fig_stub, term_col, padj_col, count_col, ratio_col, nes_col, topn):
    if nes_col:
        xline = f'x = pd.to_numeric(raw["{nes_col}"], errors="coerce").to_numpy(); xlab = "NES"'
    elif ratio_col:
        xline = (f'x = raw["{ratio_col}"].astype(str).map(lambda s: float(s.split("/")[0])/float(s.split("/")[1]) '
                 f'if "/" in s else float(s)).to_numpy(); xlab = "gene ratio"')
    else:
        xline = f'x = -np.log10(pd.to_numeric(raw["{padj_col}"], errors="coerce").to_numpy()); xlab = "-log10 adjusted p"'
    count_line = (f'ct = pd.to_numeric(raw["{count_col}"], errors="coerce").to_numpy()'
                  if count_col else "ct = np.full(len(raw), 60.0)")
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
term = raw["{term_col}"].astype(str).to_numpy()
padj = pd.to_numeric(raw["{padj_col}"], errors="coerce").to_numpy()
{count_line}
{xline}
idx = np.argsort(padj)[:{topn}]
t, x, pj, ct = term[idx], x[idx], padj[idx], ct[idx]
yo = np.argsort(x); t, x, pj, ct = t[yo], x[yo], pj[yo], ct[yo]
fig, ax = plt.subplots(figsize=(7.2, max(4.0, 0.34*len(t)+1.6)))
sc = ax.scatter(x, np.arange(len(t)), c=pj, s=70, cmap="coolwarm_r", edgecolors="none")
ax.set_yticks(np.arange(len(t))); ax.set_yticklabels(t)
ax.set_xlabel(xlab)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.colorbar(sc, ax=ax, label="p.adjust", fraction=0.046, pad=0.04)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
