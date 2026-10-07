"""Recipe: PCA scatter with 95% confidence ellipses and a PERMANOVA p. Mirrors
the R recipe: rows are samples, numeric columns are features, --group colors the
points. Features are centered and scaled, PC1 vs PC2 is drawn with per-group
ellipses, and a seeded one-way PERMANOVA tests group separation."""
import itertools
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Ellipse
from scipy import stats as ss

from ..bundle import slugify
from ..theme import apply_pub_style, palette


def _scale(X):
    return (X - X.mean(0)) / X.std(0, ddof=1)


#: Designs with this many distinct labelings or fewer are tested exactly.
MAX_EXACT = 5000


def pm_fstat(g, xs, levels, ss_t):
    """Pseudo-F of a one-way PERMANOVA on Euclidean distances.

    Args:
        g: The group label of each row.
        xs: The scaled data.
        levels: The group levels.
        ss_t: The total sum of squares.

    Returns:
        float: The pseudo-F statistic.
    """
    ssw = sum(((xs[g == lv] - xs[g == lv].mean(0)) ** 2).sum() for lv in levels if (g == lv).any())
    return ((ss_t - ssw) / (len(levels) - 1)) / (ssw / (len(xs) - len(levels)))


def pm_assignments(idx, sizes):
    """Every distinct way to split the row indices into groups of the given sizes.

    Args:
        idx: The row indices still to assign.
        sizes: The size of each remaining group.

    Returns:
        list: Each item is a list of index tuples, one tuple per group.
    """
    if len(sizes) == 1:
        return [[tuple(idx)]]
    out = []
    for first in itertools.combinations(idx, sizes[0]):
        rest_idx = [i for i in idx if i not in first]
        for rest in pm_assignments(rest_idx, sizes[1:]):
            out.append([first] + rest)
    return out


def _permanova(Xs, groups, nperm=999, seed=1):
    """One-way PERMANOVA, exact for small designs, as the R engine runs it.

    A design with MAX_EXACT distinct labelings or fewer is tested by enumerating
    every labeling, so both engines report the same P. A larger design takes
    seeded random permutations, whose P differs between the engines.

    Args:
        Xs: The scaled data.
        groups: The group label of each row.
        nperm: The count of random permutations for a large design.
        seed: The random seed.

    Returns:
        tuple: F, P, the method ("exact" or "monte carlo") and the count of
        labelings or permutations.
    """
    levels = sorted(pd.unique(groups))
    sizes = [int((groups == lv).sum()) for lv in levels]
    ss_t = float(((Xs - Xs.mean(0)) ** 2).sum())
    obs = pm_fstat(groups, Xs, levels, ss_t)
    log_n = math.lgamma(len(groups) + 1) - sum(math.lgamma(s + 1) for s in sizes)
    if math.exp(log_n) <= MAX_EXACT:
        fs = []
        for assignment in pm_assignments(list(range(len(groups))), sizes):
            g = np.empty(len(groups), dtype=object)
            for lv, members in zip(levels, assignment):
                g[list(members)] = lv
            fs.append(pm_fstat(g, Xs, levels, ss_t))
        fs = np.asarray(fs)
        return float(obs), float(np.mean(fs >= obs - 1e-9 * abs(obs))), "exact", len(fs)
    rng = np.random.RandomState(seed)
    ge = sum(pm_fstat(rng.permutation(groups), Xs, levels, ss_t) >= obs for _ in range(nperm))
    return float(obs), float((ge + 1) / (nperm + 1)), "monte carlo", nperm


def _ellipse(ax, x, y, color, level=0.95):
    """Draw the normal-theory ellipse of ggplot2's stat_ellipse(type = "norm").

    The radius is sqrt(2 * F(level; 2, n - 1)) on the sample covariance, and a
    group of fewer than four points gets no ellipse, as in ggplot2.

    Args:
        ax: The axes.
        x: The PC1 scores of the group.
        y: The PC2 scores of the group.
        color: The edge colour.
        level: The confidence level.
    """
    if len(x) < 4:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    scale = np.sqrt(2 * ss.f.ppf(level, 2, len(x) - 1))
    w, h = 2 * scale * np.sqrt(np.maximum(vals, 0))
    ax.add_patch(Ellipse((x.mean(), y.mean()), w, h, angle=theta,
                         fill=False, edgecolor=color, lw=1.0))


def recipe_pca(df, spec):
    grp_col = spec["data"].get("group") or spec["data"].get("fill") or spec["data"].get("x")
    if not grp_col:
        raise ValueError("pca needs --group (the grouping column)")
    if grp_col not in df.columns:
        raise ValueError(f"column '{grp_col}' not found in data")
    lab_col = spec["data"].get("label")

    num = df.select_dtypes(include="number")
    feat = [c for c in num.columns if c not in (grp_col, lab_col)]
    if len(feat) < 2:
        raise ValueError("pca needs at least 2 numeric feature columns")
    X = df[feat].to_numpy(dtype=float)
    keepc = X.std(0, ddof=1) > 0
    X = X[:, keepc]
    keepr = ~np.isnan(X).any(1)
    X = X[keepr]
    groups = df[grp_col].astype(str).to_numpy()[keepr]
    if len(pd.unique(groups)) < 2:
        raise ValueError("pca needs at least 2 groups")
    labs = df[lab_col].astype(str).to_numpy()[keepr] if lab_col else None
    steps = [f"Used {X.shape[0]} samples and {X.shape[1]} numeric features."]

    Xs = _scale(X)
    U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
    scores = U * S
    ve = (S ** 2) / np.sum(S ** 2) * 100
    pc1, pc2 = scores[:, 0], scores[:, 1]

    fstat, pval, pm_method, pm_count = _permanova(Xs, groups, 999, 1)

    ap = spec["appearance"]
    levels = sorted(pd.unique(groups))
    pal = palette(len(levels), ap.get("palette"))
    w, h = 5.2, 4.6
    fig, ax = plt.subplots(figsize=(w, h))
    for j, lv in enumerate(levels):
        m = groups == lv
        ax.scatter(pc1[m], pc2[m], color=pal[j], s=40, alpha=0.9, label=lv)
        _ellipse(ax, pc1[m], pc2[m], pal[j])
    if labs is not None:
        for i in range(len(pc1)):
            ax.annotate(labs[i], (pc1[i], pc2[i]), fontsize=7, color="#4d4d4d")

    ax.set_xlabel(ap.get("x_label") or f"PC1 ({ve[0]:.1f}%)")
    ax.set_ylabel(ap.get("y_label") or f"PC2 ({ve[1]:.1f}%)")
    ax.set_title(f"PERMANOVA: F = {fstat:.2f}, p = {pval:.2g}", fontsize=10, color="#4d4d4d")
    apply_pub_style(ax, ap.get("theme"))
    ax.legend(title=grp_col, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)
    fig.tight_layout()

    k = min(5, len(ve))
    stats_df = pd.DataFrame({"PC": [f"PC{i + 1}" for i in range(k)],
                             "variance_explained_pct": np.round(ve[:k], 3)})
    test_meta = {"name": "PCA with PERMANOVA",
                 "pc1_variance_pct": float(ve[0]), "pc2_variance_pct": float(ve[1]),
                 "n_samples": int(X.shape[0]), "n_features": int(X.shape[1]),
                 "permanova": {"statistic_F": fstat, "p_value": pval, "method": pm_method,
                               "permutations": pm_count}}
    methods = _methods(X.shape[0], X.shape[1], ve, fstat, pval, pm_method, pm_count)
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df[keepr].reset_index(drop=True), "clean_steps": steps,
        "label": slugify(os.path.splitext(spec["data"]["file"])[0]),
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(grp_col, in_name, fig_stub),
    }


def pm_text(method, count):
    """Describe the PERMANOVA route for the methods paragraph.

    Args:
        method: "exact" or "monte carlo".
        count: The count of labelings or permutations.

    Returns:
        str: The phrase, as the R engine writes it.
    """
    if method == "exact":
        return f"exact test over all {count} distinct labelings"
    return f"{count} seeded random permutations"


def _methods(n_samp, n_feat, ve, fstat, pval, pm_method, pm_count):
    from importlib.metadata import version as v
    import platform
    return (
        f"Principal component analysis was run on {n_samp} samples by {n_feat} centered and scaled "
        f"features. PC1 and PC2 explained {ve[0]:.1f}% and {ve[1]:.1f}% of the variance and are "
        "shown with 95% normal-theory confidence ellipses per group (drawn for groups of four "
        "samples or more). Group separation was tested with a one-way PERMANOVA on Euclidean "
        f"distances ({pm_text(pm_method, pm_count)}): F = {fstat:.2f}, p = {pval:.2g}. "
        f"Drawn in Python {platform.python_version()} with matplotlib {v('matplotlib')}."
    )


def _script(grp_col, in_name, fig_stub):
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

raw = pd.read_csv("{in_name}")
grp = raw["{grp_col}"].astype(str).to_numpy()
num = raw.select_dtypes(include="number")
X = num.to_numpy(float); X = X[:, X.std(0, ddof=1) > 0]
Xs = (X - X.mean(0)) / X.std(0, ddof=1)
U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
sc = U * S; ve = S**2 / np.sum(S**2) * 100
fig, ax = plt.subplots(figsize=(5.2, 4.6))
for lv in pd.unique(grp):
    m = grp == lv
    ax.scatter(sc[m, 0], sc[m, 1], s=40, alpha=0.9, label=lv)
ax.set_xlabel(f"PC1 ({{ve[0]:.1f}}%)"); ax.set_ylabel(f"PC2 ({{ve[1]:.1f}}%)")
ax.legend(frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
