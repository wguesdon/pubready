"""Recipe: Manhattan plot and QQ plot of GWAS summary statistics.

Mirrors the R recipe. Input is one row per SNP: the chromosome (--group), the
base-pair position (--x), the p value (--y) and, optionally, the SNP name
(--label). The chromosomes are put in natural order (a "chr" prefix is dropped;
1..22, then X, Y, XY, MT) and each SNP is drawn at its cumulative genomic
position. The genome-wide line is at 5e-8 and the suggestive line at 1e-5; the
lead SNP of each chromosome that passes 5e-8 is labelled. The QQ panel carries
a 95% band from the beta distribution of the uniform order statistics and the
genomic inflation factor lambda_GC. Every statistic uses all SNPs; the SNPs with
p >= 0.01 are thinned deterministically when there are more than 50,000.
"""
import functools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from scipy import stats

from ..bundle import slugify
from ..theme import apply_pub_style

GW = 5e-8
SUGG = 1e-5
MAX_POINTS = 50000
KEEP_BELOW = 0.01
MIN_REST = 10000
GAP_FRAC = 0.004
LABEL_SEP = 0.03
#: The grid cell of the QQ thinning, in -log10 p on both axes.
QQ_CELL = 0.005
SPECIAL = ["X", "Y", "XY", "MT", "M"]
BLUES = ["#3B6DB3", "#1F3A5F"]
#: The smallest positive normalised double, R's .Machine$double.xmin.
DOUBLE_XMIN = float(np.finfo(float).tiny)


def chr_name(v):
    """Strip the whitespace and a leading "chr" from chromosome names.

    Args:
        v: The chromosome names, as a pandas Series.

    Returns:
        pandas.Series: The names, missing values kept as NaN.
    """
    s = v.astype(str).str.strip().str.replace(r"^(?i:chr)", "", regex=True)
    return s.where(v.notna())


def chr_sort_key(name):
    """Sort key of the natural chromosome order.

    Args:
        name: One chromosome name, without the "chr" prefix.

    Returns:
        tuple: The tier (number, special, other), the rank and the name.
    """
    if name.isdigit():
        return (1, float(name), "")
    up = name.upper()
    if up in SPECIAL:
        return (2, float(SPECIAL.index(up) + 1), "")
    return (3, 0.0, name)


def thin_mask(p, max_points=MAX_POINTS, keep_below=KEEP_BELOW, min_rest=MIN_REST):
    """Keep every point below keep_below and every k-th of the rest.

    Args:
        p: The p values, in drawing order.
        max_points: The approximate count of points to draw.
        keep_below: The p below which every point is kept.
        min_rest: The smallest budget for the other points.

    Returns:
        tuple: The boolean mask and the step k (1 keeps all).
    """
    n = len(p)
    if n <= max_points:
        return np.ones(n, dtype=bool), 1
    low = p < keep_below
    rest = np.where(~low)[0]
    k = int(np.ceil(len(rest) / max(max_points - int(low.sum()), min_rest)))
    keep = low.copy()
    keep[rest[np.arange(len(rest)) % k == 0]] = True
    return keep, k


def label_keep(centres, min_sep):
    """Keep a chromosome label only when it is far enough from the last one.

    Args:
        centres: The label positions, in order.
        min_sep: The smallest distance between two kept labels.

    Returns:
        numpy.ndarray: A boolean mask of the kept labels.
    """
    keep = np.zeros(len(centres), dtype=bool)
    last = -np.inf
    for i, c in enumerate(centres):
        if c - last >= min_sep:
            keep[i] = True
            last = c
    return keep


def display_y(p, is_zero):
    """-log10 p for drawing; a p of exactly 0 sits at a ceiling.

    The ceiling is 12% above the largest other value, or above the genome-wide
    line, as in the R engine.

    Args:
        p: The p values, zeros already set to the smallest positive double.
        is_zero: True where the input p was exactly 0.

    Returns:
        numpy.ndarray: The values to draw.
    """
    y = -np.log10(p)
    top = float(y[~is_zero].max()) if (~is_zero).any() else 0.0
    y[is_zero] = max(top, -np.log10(GW)) * 1.12
    return y


def short_label(s):
    """Cut a SNP label longer than 20 characters; the stats keep the full name.

    Args:
        s: The label.

    Returns:
        str: The label, or its first 19 characters and an ellipsis.
    """
    s = str(s)
    return s if len(s) <= 20 else s[:19] + "…"


def qq_data(p, is_zero):
    """The QQ points and the 95% band of the order statistics, thinned as drawn.

    The sorted points lie on one curve, so every point with p < 0.01 is drawn
    and, for the others, the first point in each cell of a QQ_CELL grid in
    -log10 p on both axes, as in the R engine.

    Args:
        p: The p values of every SNP.
        is_zero: True where the input p was exactly 0.

    Returns:
        dict: The expected and observed -log10 p, the zero flags and the band
        at the drawn points.
    """
    n = len(p)
    o = np.argsort(p, kind="stable")
    ps, zs = p[o], is_zero[o]
    i = np.arange(1, n + 1)
    expd = -np.log10((i - 0.5) / n)
    obs = display_y(ps, zs)
    key = np.floor(expd / QQ_CELL) * 1e6 + np.floor(obs / QQ_CELL)
    keep = ps < KEEP_BELOW
    keep[np.unique(key, return_index=True)[1]] = True
    j = i[keep]
    return {"exp": expd[keep], "obs": obs[keep], "zero": zs[keep],
            "lo": -np.log10(stats.beta.ppf(0.975, j, n - j + 1)),
            "hi": -np.log10(stats.beta.ppf(0.025, j, n - j + 1)),
            "exp_max": float(expd.max())}


def fmt_count(n):
    """Format a count with thousands separators, as R's big.mark does.

    Args:
        n: The count.

    Returns:
        str: The formatted count.
    """
    return f"{int(n):,}"


def recipe_manhattan(df, spec):
    """Draw the Manhattan plot and the QQ plot side by side.

    Args:
        df: The summary statistics, one row per SNP.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing or a p value lies outside 0 to 1.
    """
    dat = spec["data"]
    chr_col, pos_col, p_col, lab_col = dat.get("group"), dat.get("x"), dat.get("y"), dat.get("label")
    if not chr_col:
        raise ValueError("manhattan needs --group (the chromosome column)")
    if not pos_col:
        raise ValueError("manhattan needs --x (the base-pair position column)")
    if not p_col:
        raise ValueError("manhattan needs --y (the p-value column)")
    for c in [chr_col, pos_col, p_col] + ([lab_col] if lab_col else []):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    raw_chr = df[chr_col].astype(str).str.strip().where(df[chr_col].notna())
    chr_s = chr_name(df[chr_col])
    pos = pd.to_numeric(df[pos_col], errors="coerce").to_numpy(dtype=float)
    p = pd.to_numeric(df[p_col], errors="coerce").to_numpy(dtype=float)
    bad = int(np.sum((p < 0) | (p > 1)))
    if bad:
        raise ValueError(f"manhattan: p must lie in (0, 1]; column '{p_col}' holds {bad} "
                         "value(s) outside 0 to 1")
    keep = (chr_s.notna() & (chr_s.str.len() > 0)).to_numpy() & np.isfinite(pos) & np.isfinite(p)
    had_prefix = bool(raw_chr[keep].str.match(r"(?i)^chr").any())
    chrs = chr_s.to_numpy()[keep].astype(str)
    pos, p = pos[keep], p[keep]
    snp = (df[lab_col].astype(str).to_numpy()[keep] if lab_col
           else np.array([f"{c}:{x:.0f}" for c, x in zip(chrs, pos)]))
    if len(p) < 2:
        raise ValueError("manhattan needs at least 2 SNPs with a chromosome, a position and a p value")
    is_zero = p == 0
    n_zero = int(is_zero.sum())
    p = np.where(is_zero, DOUBLE_XMIN, p)
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing chromosome, position or p "
                     f"({n0} -> {len(p)}).")
    if had_prefix:
        steps.append("Read the chromosome names without their 'chr' prefix.")
    if n_zero:
        steps.append(f"Replaced {n_zero} p value(s) of exactly 0 with the smallest positive double "
                     f"({DOUBLE_XMIN:.3g}), so that -log10 p is finite.")
    if not steps:
        steps = ["No cleaning applied; input used as-is."]

    n = len(p)
    lam = float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1))
    n_gw, n_sugg = int((p < GW).sum()), int((p < SUGG).sum())
    topn = int(spec["appearance"].get("top_n", 15))

    lv = sorted(pd.unique(chrs), key=chr_sort_key)
    index = {c: k for k, c in enumerate(lv)}
    ci = np.array([index[c] for c in chrs])
    idx = np.arange(n)
    o = np.lexsort((idx, p, ci))
    first = o[np.r_[True, ci[o][1:] != ci[o][:-1]]]
    leads = first[p[first] < GW]
    leads = leads[np.lexsort((leads, ci[leads], p[leads]))][:topn]

    lo = np.array([pos[ci == k].min() for k in range(len(lv))])
    hi = np.array([pos[ci == k].max() for k in range(len(lv))])
    span = hi - lo
    gap = GAP_FRAC * max(span.sum(), 1.0)
    offset = np.concatenate([[0.0], np.cumsum(span + gap)])[:len(lv)]
    cum = pos - lo[ci] + offset[ci]
    centre = offset + span / 2
    lab_ok = label_keep(centre, LABEL_SEP * max(cum.max() - cum.min(), 1.0))

    # A p of exactly 0 has no finite -log10 p of its own; it is drawn as a
    # triangle at a ceiling just above the largest other value.
    draw, k_step = thin_mask(p)
    ydisp = display_y(p, is_zero)
    ymax = max(float(ydisp.max()), -np.log10(GW)) * (1.14 if len(leads) else 1.06)
    small = int(draw.sum()) > 20000
    s_pt = 2.6 if small else 5.9

    ap = spec["appearance"]
    w, h = 10.5, 4.3
    fig, (ax, aq) = plt.subplots(1, 2, figsize=(w, h), gridspec_kw={"width_ratios": [2.7, 1]})
    sel = draw & ~is_zero
    ax.scatter(cum[sel], ydisp[sel], s=s_pt, linewidths=0,
               c=np.where(ci[sel] % 2 == 0, BLUES[0], BLUES[1]), zorder=1)
    ax.axhline(-np.log10(SUGG), ls="--", lw=0.8, color="#8C8C8C", zorder=2)
    ax.axhline(-np.log10(GW), ls="--", lw=0.9, color="#C1432B", zorder=2)
    if n_zero:
        ax.scatter(cum[is_zero], ydisp[is_zero], s=20, marker="^", c="#C1432B", linewidths=0, zorder=3)
    lead_round = leads[~is_zero[leads]]
    if len(lead_round):
        ax.scatter(cum[lead_round], ydisp[lead_round], s=18, c="#C1432B", linewidths=0, zorder=3)
    ax.set_xticks(centre[lab_ok])
    ax.set_xticklabels([lv[k] for k in np.where(lab_ok)[0]], fontsize=8.5)
    pad = 0.01 * (cum.max() - cum.min())
    ax.set_xlim(cum.min() - pad, cum.max() + pad)
    ax.set_ylim(0, ymax * 1.02)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10], integer=True))
    texts = [ax.text(cum[i], ydisp[i], short_label(snp[i]), fontsize=8, zorder=5) for i in leads]
    if texts:
        try:
            from adjustText import adjust_text
            adjust_text(texts, x=cum[leads], y=ydisp[leads], ax=ax, expand=(1.2, 1.6),
                        force_text=(0.5, 0.8), iter_lim=2000,
                        arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.5})
        except Exception:
            pass
    ax.set_xlabel(ap.get("x_label") or "Chromosome")
    ax.set_ylabel(ap.get("y_label") or "-log10 p")
    zero_txt = f"; {n_zero} with p = 0 (triangles)" if n_zero else ""
    ax.set_title(f"{fmt_count(n)} SNPs; {fmt_count(n_gw)} with p < 5e-8; "
                 f"{fmt_count(n_sugg)} with p < 1e-5{zero_txt}", fontsize=10, color="#4d4d4d")
    apply_pub_style(ax, ap.get("theme"))
    ax.tick_params(axis="x", labelsize=8.5)

    qq = qq_data(p, is_zero)
    qlim = max(qq["exp_max"], float(qq["obs"].max())) * 1.04
    aq.fill_between(qq["exp"], qq["lo"], qq["hi"], color="#D9D9D9", linewidth=0, zorder=1)
    aq.plot([0, qlim], [0, qlim], color="#C1432B", lw=0.9, zorder=2)
    nz = ~qq["zero"]
    aq.scatter(qq["exp"][nz], qq["obs"][nz], s=s_pt, c=BLUES[0], linewidths=0, zorder=3)
    if qq["zero"].any():
        aq.scatter(qq["exp"][~nz], qq["obs"][~nz], s=20, marker="^", c="#C1432B", linewidths=0,
                   zorder=4)
    aq.set_xlim(0, qq["exp_max"] * 1.04 * 1.02)
    aq.set_ylim(0, qlim * 1.02)
    aq.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10], integer=True))
    aq.yaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10], integer=True))
    aq.set_xlabel("Expected -log10 p")
    aq.set_ylabel("Observed -log10 p")
    aq.set_title(rf"$\lambda_{{GC}}$ = {lam:.3f}", fontsize=10, color="#4d4d4d")
    apply_pub_style(aq, ap.get("theme"))
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    fig.tight_layout()

    genome = {"row_type": "genome", "chr": None, "snp": None, "pos": np.nan,
              "p_value": float(p.min()), "n_snps": n, "lambda_gc": lam,
              "n_genome_wide": n_gw, "n_suggestive": n_sugg}
    rows = [genome] + [{"row_type": "lead", "chr": lv[ci[i]], "snp": snp[i], "pos": pos[i],
                        "p_value": p[i], "n_snps": None, "lambda_gc": np.nan,
                        "n_genome_wide": None, "n_suggestive": None} for i in leads]
    stats_df = pd.DataFrame(rows, columns=list(genome))
    test_meta = {"name": "Manhattan and QQ plot", "n_snps": n, "n_chromosomes": len(lv),
                 "lambda_gc": lam, "n_genome_wide": n_gw, "n_suggestive": n_sugg,
                 "min_p": float(p.min()), "genome_wide_threshold": GW,
                 "suggestive_threshold": SUGG, "n_zero_p": n_zero, "n_leads": int(len(leads)),
                 "thinning_step": int(k_step), "n_drawn": int(draw.sum())}
    df_used = pd.DataFrame({"chr": chrs, "pos": pos, "p": p})
    if lab_col:
        df_used["snp"] = snp
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": _methods(test_meta),
        "df_used": df_used, "clean_steps": steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]), "width": w, "height": h,
        "build_script": functools.partial(_script, chr_col, pos_col, p_col),
    }


def _methods(tm):
    """Write the methods paragraph, as the R engine words it.

    Args:
        tm: The test metadata.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    if tm["thinning_step"] > 1:
        thin = (f"For the size of the vector outputs, every SNP with p < 0.01 is drawn and one in "
                f"every {tm['thinning_step']} of the other SNPs in input order "
                f"({fmt_count(tm['n_drawn'])} of {fmt_count(tm['n_snps'])} points drawn); every "
                "statistic uses all SNPs. ")
    else:
        thin = f"All {fmt_count(tm['n_snps'])} SNPs are drawn in the Manhattan panel. "
    thin += ("The QQ panel draws every SNP with p < 0.01 and, for the others, one point per "
             "0.005 by 0.005 cell of -log10 p, because the sorted points lie on one curve. ")
    zero = (f"{tm['n_zero_p']} p value(s) of exactly 0 were set to the smallest positive double "
            "for the statistics and are drawn as triangles at the top of both panels. "
            if tm["n_zero_p"] else "")
    return (
        f"Association results for {fmt_count(tm['n_snps'])} SNPs on {tm['n_chromosomes']} "
        "chromosomes were shown as a Manhattan plot of -log10 p at the cumulative genomic "
        "position, with the genome-wide significance line at p = 5e-8 and the suggestive line at "
        f"p = 1e-5. {fmt_count(tm['n_genome_wide'])} SNPs reached genome-wide significance and "
        f"{fmt_count(tm['n_suggestive'])} passed the suggestive line (smallest p = "
        f"{tm['min_p']:.3g}); the lead SNP of each chromosome with p < 5e-8 is labelled "
        f"({tm['n_leads']} shown). {zero}The QQ plot compares the observed -log10 p with the "
        "expected uniform quantiles, (i - 0.5) / n, with a 95% band from the Beta(i, n - i + 1) "
        "distribution of the order statistics. The genomic inflation factor was lambda_GC = "
        f"{tm['lambda_gc']:.3f}, the median of the 1-df chi-square statistics over its expected "
        f"value of 0.455. {thin}The association test itself was performed upstream; this figure "
        f"summarises its output. Drawn in Python {platform.python_version()} with matplotlib "
        f"{v('matplotlib')}, scipy {v('scipy')} and adjustText {v('adjustText')}."
    )


def _script(chr_col, pos_col, p_col, in_name, fig_stub):
    """Write the standalone script that redraws both panels.

    Args:
        chr_col: The chromosome column.
        pos_col: The position column.
        p_col: The p-value column.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from scipy import stats

raw = pd.read_csv("{in_name}")
chrs = raw["{chr_col}"].astype(str).str.replace(r"^(?i:chr)", "", regex=True).to_numpy()
pos = pd.to_numeric(raw["{pos_col}"], errors="coerce").to_numpy(float)
p = pd.to_numeric(raw["{p_col}"], errors="coerce").to_numpy(float)
keep = np.isfinite(pos) & np.isfinite(p)
chrs, pos, p = chrs[keep], pos[keep], np.where(p[keep] == 0, np.finfo(float).tiny, p[keep])
special = ["X", "Y", "XY", "MT", "M"]
key = {{c: ((1, float(c), "") if c.isdigit() else (2, float(special.index(c.upper()) + 1), "")
           if c.upper() in special else (3, 0.0, c)) for c in set(chrs)}}
lv = sorted(set(chrs), key=key.get)
ci = np.array([lv.index(c) for c in chrs])
lo = np.array([pos[ci == k].min() for k in range(len(lv))])
span = np.array([pos[ci == k].max() for k in range(len(lv))]) - lo
offset = np.concatenate([[0.0], np.cumsum(span + 0.004 * span.sum())])[:len(lv)]
x = pos - lo[ci] + offset[ci]
lam = np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)
fig, (ax, aq) = plt.subplots(1, 2, figsize=(10.5, 4.3), gridspec_kw={{"width_ratios": [2.7, 1]}})
ax.scatter(x, -np.log10(p), s=4, c=np.where(ci % 2 == 0, "#3B6DB3", "#1F3A5F"), linewidths=0)
ax.axhline(-np.log10(1e-5), ls="--", color="grey"); ax.axhline(-np.log10(5e-8), ls="--", color="#C1432B")
ax.set_xticks(offset + span / 2); ax.set_xticklabels(lv, fontsize=8)
ax.set_xlabel("Chromosome"); ax.set_ylabel("-log10 p")
n = len(p); i = np.arange(1, n + 1)
e = -np.log10((i - 0.5) / n)
aq.fill_between(e, -np.log10(stats.beta.ppf(0.975, i, n - i + 1)),
                -np.log10(stats.beta.ppf(0.025, i, n - i + 1)), color="#D9D9D9")
aq.plot([0, e.max()], [0, e.max()], color="#C1432B")
aq.scatter(e, -np.log10(np.sort(p)), s=4, c="#3B6DB3", linewidths=0)
aq.set_xlabel("Expected -log10 p"); aq.set_ylabel("Observed -log10 p")
aq.set_title(f"lambda_GC = {{lam:.3f}}")
for a in (ax, aq):
    for s in ("top", "right"): a.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
