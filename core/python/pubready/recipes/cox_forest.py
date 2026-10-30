"""Recipe: Cox proportional-hazards model with a forest plot of hazard ratios
(lifelines). Mirrors the R cox_forest: HR with 95% CI and p per covariate, and a
proportional-hazards check.

The figure matches the R engine row for row: the label column on the left, one
row per model term (the first level of each factor gets its own "reference"
row), squares with 95% CI whiskers on a log axis, and the "HR (95% CI)" and "P"
columns on the right.
"""
import functools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from matplotlib import transforms
from matplotlib.ticker import NullLocator

from ..clean import coerce_numeric_columns
from ..theme import apply_pub_style, palette, plain_log_axis
from ..util import format_p, level_labels, level_order, p_phrase
from .survival_km import event_codes

# Tick sequences of the hazard-ratio axis, from coarse to fine (mantissas per decade).
TICK_SETS = [(1, 2, 5), (1, 1.5, 2, 3, 5, 7), (1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8)]


def forest_stars(p):
    """The significance stars of the forest plot, the same in both engines.

    The rule stops at three stars, as survminer::ggforest graded them, so a P
    below 1e-4 keeps "***".

    Args:
        p: The P value of one term.

    Returns:
        str: "***", "**", "*" or an empty string.
    """
    if p is None or np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def ticks_for(mant, lo, hi):
    """The values of one tick sequence in [lo, hi], extended across 1 when needed.

    On a side of 1 that the CIs reach but no tick covers, the next value of the
    sequence outward is added.

    Args:
        mant: The mantissas of the sequence in one decade, for example (1, 2, 5).
        lo: The smallest value on the axis, at most 1.
        hi: The largest value on the axis, at least 1.

    Returns:
        list: The tick positions in increasing order.
    """
    seq = sorted(m * 10.0 ** k for k in range(-4, 5) for m in mant)
    ticks = [v for v in seq if lo <= v <= hi]
    below = [v for v in seq if v < lo]
    above = [v for v in seq if v > hi]
    if lo < 1 and not any(v < 1 for v in ticks) and below:
        ticks = [below[-1]] + ticks
    if hi > 1 and not any(v > 1 for v in ticks) and above:
        ticks = ticks + [above[0]]
    return ticks


def forest_log_ticks(lo, hi):
    """Ticks of the hazard-ratio axis that cover the CIs on both sides of 1.

    The 1-2-5 sequence comes first. It is kept when it covers the CIs on both
    sides of 1 without widening the axis by more than 60%. A narrow range falls
    back to finer steps, and a range over many decades to the powers of ten. The
    R engine uses the same rule.

    Args:
        lo: The smallest lower CI bound.
        hi: The largest upper CI bound.

    Returns:
        list: The tick positions in increasing order, 1 included.
    """
    lo, hi = min(lo, 1.0), max(hi, 1.0)
    span = np.log(hi / lo)
    sets = TICK_SETS
    if len(ticks_for((1, 2, 5), lo, hi)) > 8:
        sets = [(1,)]
    ticks = [1.0]
    for mant in sets:
        ticks = ticks_for(mant, lo, hi)
        total = np.log(max([hi] + ticks) / min([lo] + ticks))
        if total <= 1.6 * span:
            return ticks
    return ticks


def label_inches(labels):
    """The width of the label column in inches, from the longest label.

    Args:
        labels: The row labels.

    Returns:
        float: The width, the same estimate as the R engine.
    """
    return 0.3 + 0.085 * max([len(s) for s in labels] + [4])


def left_km(times, events):
    """The Kaplan-Meier survival just before each time, S(t-).

    survival::cox.zph reads the left-continuous curve of survfitKM for its "km"
    transform: the product of (1 - d / n) over the event times below t.

    Args:
        times: The follow-up times.
        events: 1 for an event, 0 for a censored subject.

    Returns:
        numpy.ndarray: S(t-) for each subject.
    """
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=float)
    event_times = np.unique(times[events == 1])
    steps = np.ones(len(event_times) + 1)
    for i, u in enumerate(event_times):
        n_risk = np.sum(times >= u)
        d = np.sum((times == u) & (events == 1))
        steps[i + 1] = steps[i] * (1 - d / n_risk)
    return steps[np.searchsorted(event_times, times, side="left")]


def cox_zph_global(times, events, X, beta):
    """The global test of survival::cox.zph with the default "km" transform.

    The score test of beta(t) = beta + theta g(t) at theta = 0, with g(t) one
    minus the left-continuous Kaplan-Meier curve, centred on the event times.
    The sums follow the C routine zph1 of the survival package with Efron ties:
    the score u and the information of the 2p parameters, then u' I^-1 u on p
    degrees of freedom with the first p entries of u set to zero.

    Args:
        times: The follow-up times.
        events: 1 for an event, 0 for a censored subject.
        X: The design matrix of the model, one column per coefficient.
        beta: The fitted coefficients, in the column order of X.

    Returns:
        tuple: The chi-square statistic, its degrees of freedom and the P value.
    """
    from scipy import stats

    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=float)
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    eta = X @ np.asarray(beta, dtype=float)
    risk = np.exp(eta - eta.max())  # a common factor cancels in every ratio
    ttimes = 1 - left_km(times, events)
    gt = ttimes - ttimes[events == 1].mean()
    xc = X - X.mean(axis=0)
    u = np.zeros(2 * p)
    imat = np.zeros((2 * p, 2 * p))
    a = np.zeros(p)
    cmat = np.zeros((p, p))
    denom = 0.0
    for tval in np.unique(times)[::-1]:
        idx = np.where(times == tval)[0]
        dead = idx[events[idx] == 1]
        for i in idx[events[idx] == 0]:
            denom += risk[i]
            a += risk[i] * xc[i]
            cmat += risk[i] * np.outer(xc[i], xc[i])
        nd = len(dead)
        if nd == 0:
            continue
        timewt = gt[idx[0]]
        denom2 = risk[dead].sum()
        a2 = (risk[dead, None] * xc[dead]).sum(axis=0)
        cmat2 = (risk[dead, None, None] * np.einsum("ki,kj->kij", xc[dead], xc[dead])).sum(axis=0)
        u[:p] += xc[dead].sum(axis=0)
        u[p:] += timewt * xc[dead].sum(axis=0)
        for _ in range(nd):
            denom += denom2 / nd
            a += a2 / nd
            cmat += cmat2 / nd
            mean = a / denom
            u[:p] -= mean
            u[p:] -= timewt * mean
            v = (cmat - np.outer(mean, a)) / denom
            imat[:p, :p] += v
            imat[:p, p:] += timewt * v
            imat[p:, p:] += timewt * timewt * v
    imat[p:, :p] = imat[:p, p:].T
    score = np.concatenate([np.zeros(p), u[p:]])
    chisq = float(score @ np.linalg.solve(imat, score))
    return chisq, p, float(stats.chi2.sf(chisq, p))


def design_and_rows(df, t, e, covs):
    """Build the numeric design matrix and the display rows of the forest plot.

    A categorical covariate is one-hot encoded with its first level, by the
    level rule that both engines share (level_order), as the reference. It
    gives a "reference" row plus one row per other level. A numeric covariate
    gives one row.

    Args:
        df: The cleaned data.
        t: The time column.
        e: The event column.
        covs: The covariate columns.

    Returns:
        tuple: The design data frame and a list of row dicts with the keys
        ``label``, ``term`` (None for a reference row) and ``covariate``.
    """
    design = df[[t, e]].copy()
    rows = []
    for c in covs:
        s = df[c]
        if s.dtype == object or str(s.dtype).startswith("category") or not np.issubdtype(s.dtype, np.number):
            cats = level_order(s.astype(str))
            rows.append({"label": f"{c}: {cats[0]}", "term": None, "covariate": c})
            for lvl in cats[1:]:
                design[f"{c}{lvl}"] = (s.astype(str) == lvl).astype(int)
                rows.append({"label": f"{c}: {lvl}", "term": f"{c}{lvl}", "covariate": c})
        else:
            design[c] = pd.to_numeric(s, errors="coerce")
            rows.append({"label": c, "term": c, "covariate": c})
    return design, rows


def check_covariates(df, e, covs):
    """Check the data before the fit, with the rules and messages of the R engine.

    Args:
        df: The cleaned data.
        e: The event column.
        covs: The covariate columns.

    Raises:
        ValueError: If no row is left, if no event happened, or if a covariate
            holds one value.
    """
    if df.empty:
        raise ValueError("cox_forest found no row with a time, an event and a value of every "
                         "covariate")
    if df[e].sum() == 0:
        raise ValueError(f"cox_forest needs at least one event, and column '{e}' holds no 1 "
                         f"among the {len(df)} subject(s) kept")
    for c in covs:
        s = df[c]
        numeric = pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)
        u = list(pd.unique(level_labels(s) if numeric else s.astype(str)))
        if len(u) < 2:
            raise ValueError(f"covariate '{c}' holds one value ('{u[0]}') in the {len(df)} "
                             "subject(s) kept, so the Cox model cannot estimate its effect; "
                             "leave it out of --covariates")


def collinear_columns(x, tol=1e-8):
    """The model columns that are a linear combination of a constant and earlier columns.

    Each centred column is reduced against the orthonormal basis of the columns
    kept so far (modified Gram-Schmidt, two passes); a remainder below 1e-8 of
    its length marks it. The R engine runs the same arithmetic.

    Args:
        x: The model matrix, one column per term, without an intercept.
        tol: The largest remainder that marks a column.

    Returns:
        list: The 0-based indices of the marked columns.
    """
    x = np.asarray(x, dtype=float)
    hits, basis = [], []
    for j in range(x.shape[1]):
        v = x[:, j] - x[:, j].mean()
        nv = np.sqrt(np.sum(v ** 2))
        if nv == 0:
            hits.append(j)
            continue
        v = v / nv
        for _ in range(2):
            for q in basis:
                v = v - q * np.sum(q * v)
        r = np.sqrt(np.sum(v ** 2))
        if r < tol:
            hits.append(j)
        else:
            basis.append(v / r)
    return hits


def check_collinear(design, rows):
    """Stop when a model column is collinear with the columns before it.

    Args:
        design: The design data frame from design_and_rows.
        rows: The display rows from design_and_rows.

    Raises:
        ValueError: If a covariate is a linear combination of the ones before it.
    """
    terms = [r for r in rows if r["term"] is not None]
    hit = collinear_columns(design[[r["term"] for r in terms]].to_numpy(float))
    if hit:
        raise ValueError(f"covariate '{terms[hit[0]]['covariate']}' is collinear with the "
                         "covariates before it in --covariates (a linear combination of them), "
                         "so the Cox model cannot estimate its effect; leave it or one of them out")


def cox_separation(time, event, x):
    """The model columns that separate the events, as cox_separation() in R.

    The partial likelihood of such a column grows without bound, so its hazard
    ratio has no finite estimate. A column separates upward when, at every event
    time, each subject with an event holds the highest value of the risk set
    (the subjects whose time is that time or later), and the risk set holds a
    lower value at one event time at least; downward with the lowest value.

    Args:
        time: The follow-up times.
        event: 1 for an event, 0 for a censored subject.
        x: The model matrix, one column per term.

    Returns:
        numpy.ndarray: 1, -1 or 0 for each column.
    """
    x = np.asarray(x, dtype=float)
    event_times = np.unique(time[event == 1])
    out = np.zeros(x.shape[1], dtype=int)
    for j in range(x.shape[1]):
        at_max, at_min, spread = True, True, False
        for t in event_times:
            r = x[time >= t, j]
            d = x[(time == t) & (event == 1), j]
            if np.any(d < r.max()):
                at_max = False
            if np.any(d > r.min()):
                at_min = False
            if r.min() < r.max():
                spread = True
        if spread:
            out[j] = 1 if at_max else (-1 if at_min else 0)
    return out


def separation_text(terms, sep):
    """One sentence for each term that separates the events, as in the R engine.

    Args:
        terms: The model terms.
        sep: The separation direction of each term, 1, -1 or 0.

    Returns:
        str: The sentences, or an empty string.
    """
    return "".join(
        f" The term {tt} was not estimable because it separates the events: at every event "
        f"time, the subjects with an event held the {'highest' if s > 0 else 'lowest'} value of "
        f"{tt} in the risk set, so the hazard ratio has no finite estimate. The term stays in the "
        "model, and the figure and the stats table give no estimate for it."
        for tt, s in zip(terms, sep) if s != 0)


def row_texts(row, stats):
    """The "HR (95% CI)" and "P" texts of one display row.

    Args:
        row: The row dict.
        stats: The stats table indexed by term.

    Returns:
        tuple: The HR text and the P text.
    """
    if row["term"] is None:
        return "reference", ""
    r = stats.loc[row["term"]]
    if r.get("note", "estimated") != "estimated":
        return "not estimable, separation", ""
    star = forest_stars(r["p_value"])
    p_txt = format_p(r["p_value"]) + (" " + star if star else "")
    return f"{r['HR']:.2f} ({r['ci_low']:.2f}–{r['ci_high']:.2f})", p_txt


def cox_subtitle(n, events, ph_p):
    """Write the subtitle: the counts and, when the PH check fails, a warning.

    The warning keeps the figure from showing the hazard ratios of a model
    that breaks the proportional-hazards assumption as valid.

    Args:
        n: The number of subjects.
        events: The number of events.
        ph_p: The P value of the proportional-hazards check, or NaN.

    Returns:
        str: The subtitle.
    """
    txt = f"{n} subjects, {events} events"
    if ph_p is not None and np.isfinite(ph_p) and ph_p < 0.05:
        txt += f"; PH assumption violated ({p_phrase(ph_p)})"
    return txt


def forest_plot(rows, stats_df, spec, n, events, ph_p=np.nan):
    """Draw the forest plot with the label column and the two text columns.

    Args:
        rows: The display rows from ``design_and_rows``.
        stats_df: The stats table.
        spec: The figure specification.
        n: The number of subjects.
        events: The number of events.
        ph_p: The P value of the proportional-hazards check, or NaN.

    Returns:
        dict: The figure, its ticks, the label width, the width and the height.
    """
    ap = spec["appearance"]
    title = ap.get("title")
    pal = palette(1, ap.get("palette"))
    stats = stats_df.set_index("term")
    k = len(rows)
    yy = np.arange(k, 0, -1).astype(float)
    est = [(y_, stats.loc[r["term"]]) for y_, r in zip(yy, rows) if r["term"] is not None]
    est = [(y_, s) for y_, s in est if np.isfinite(s["HR"])]
    lows = [s["ci_low"] for _, s in est if np.isfinite(s["ci_low"]) and s["ci_low"] > 0]
    highs = [s["ci_high"] for _, s in est if np.isfinite(s["ci_high"])]
    lo = min(lows) if lows else 1.0
    hi = max(highs) if highs else 1.0
    ticks = forest_log_ticks(lo, hi)
    xa, xb = min([lo] + ticks), max([hi] + ticks)
    r = np.log10(xb / xa) if xb > xa else np.log10(2.0)
    ylim = (0.4, k + 1.6)

    lab_in = label_inches([row["label"] for row in rows])
    w = lab_in + 3.0 + 2.5
    h = max(3.0, 0.3 * (k + 1) + 1.6 + (0.3 if title else 0.0))
    top_in = 0.45 + (0.3 if title else 0.0)
    bottom_in = 0.85
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([lab_in / w, bottom_in / h, 2.85 / w, (h - top_in - bottom_in) / h])
    axt = fig.add_axes([(lab_in + 3.0) / w, bottom_in / h, 2.3 / w, (h - top_in - bottom_in) / h])

    ax.axvline(1.0, ls=":", color="0.4", lw=0.9, zorder=1)
    for y_, s in est:
        ax.hlines(y_, s["ci_low"], s["ci_high"], color="0.15", lw=1.0, zorder=2)
    ax.scatter([s["HR"] for _, s in est], [y_ for y_, _ in est], marker="s", s=48,
               color=pal[0], zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(10 ** (np.log10(xa) - 0.05 * r), 10 ** (np.log10(xb) + 0.05 * r))
    ax.set_xticks(ticks)
    plain_log_axis(ax.xaxis)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(*ylim)
    ax.set_yticks([])
    ax.set_xlabel(ap.get("x_label") or "Hazard ratio (95% CI, log scale)")
    ax.set_title(cox_subtitle(n, events, ph_p), fontsize=9.5, color="#4d4d4d", loc="left")
    apply_pub_style(ax, ap.get("theme"))
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=10.4)

    trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
    for y_, row in zip(yy, rows):
        fig.text(0.15 / w, y_, row["label"], transform=trans, ha="left", va="center",
                 fontsize=10.5)

    axt.set_xlim(0, 1)
    axt.set_ylim(*ylim)
    axt.axis("off")
    for y_, row in zip(yy, rows):
        hr_txt, p_txt = row_texts(row, stats)
        axt.text(0, y_, hr_txt, ha="left", va="center", fontsize=10)
        axt.text(0.66, y_, p_txt, ha="left", va="center", fontsize=10)
    for x_, head in ((0, "HR (95% CI)"), (0.66, "P")):
        axt.text(x_, k + 1, head, ha="left", va="center", fontsize=10, fontweight="bold")
    if title:
        fig.suptitle(title, fontweight="bold", fontsize=13)
    return {"fig": fig, "ticks": ticks, "label_in": lab_in, "width": w, "height": h}


def recipe_cox_forest(df, spec):
    t, e = spec["data"]["time"], spec["data"]["event"]
    covs = spec["data"].get("covariates")
    if not t or not e:
        raise ValueError("cox_forest needs --time and --event columns")
    if not covs:
        covs = [c for c in df.columns if c not in (t, e)]
    for c in [t, e] + covs:
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found")

    df = df.copy()
    df[t] = pd.to_numeric(df[t], errors="coerce")
    df[e] = event_codes(df[e], e)
    # A covariate whose cells read as numbers for 80% or more is numeric, and each
    # other cell, for example "45 y", becomes NaN and is logged.
    df, _, _, clean_steps = coerce_numeric_columns(df, covs)
    n0 = len(df)
    df = df.dropna(subset=[t, e] + covs)
    df = df[np.isfinite(df[t]) & (df[t] >= 0)].reset_index(drop=True)
    if len(df) < n0:
        clean_steps.append(f"Dropped {n0 - len(df)} row(s) with missing/invalid values "
                           f"({n0} -> {len(df)} rows).")
    if not clean_steps:
        clean_steps = ["No cleaning applied; input used as-is."]

    check_covariates(df, e, covs)
    # Build a numeric design matrix; one-hot encode categoricals (drop first = reference).
    design, rows = design_and_rows(df, t, e, covs)
    check_collinear(design, rows)

    cph = CoxPHFitter()
    # A tight precision makes lifelines converge to the estimate of survival::coxph;
    # the default stops about 2e-6 short of it. Terms are named as R names them.
    cph.fit(design, duration_col=t, event_col=e, fit_options={"precision": 1e-11})
    summ = cph.summary
    terms = list(summ.index)
    hr = summ["exp(coef)"].to_numpy()
    lo = summ["exp(coef) lower 95%"].to_numpy()
    hi = summ["exp(coef) upper 95%"].to_numpy()
    pv = summ["p"].to_numpy()

    # The global test of survival::cox.zph, as the R engine reports it.
    ph_p = np.nan
    try:
        coef = cph.params_
        ph_p = cox_zph_global(design[t].to_numpy(), design[e].to_numpy(),
                              design[list(coef.index)].to_numpy(), coef.to_numpy())[2]
    except (np.linalg.LinAlgError, ValueError):
        pass

    stats_df = pd.DataFrame({"term": terms, "HR": hr, "ci_low": lo, "ci_high": hi,
                             "z": summ["z"].to_numpy(), "p_value": pv})
    # A term that separates the events has no finite estimate. Its numbers are the
    # last step of the iteration, so the table, the figure and the text show none.
    sep = cox_separation(design[t].to_numpy(float), design[e].to_numpy(float),
                         design[terms].to_numpy(float))
    stats_df["note"] = np.where(sep != 0, "not estimable, separation", "estimated")
    stats_df.loc[sep != 0, ["HR", "ci_low", "ci_high", "z", "p_value"]] = np.nan
    hr, lo, hi, pv = (stats_df[k].to_numpy() for k in ("HR", "ci_low", "ci_high", "p_value"))
    n_subj, n_events = int(len(df)), int(df[e].sum())
    fp = forest_plot(rows, stats_df, spec, n_subj, n_events, ph_p)

    test_meta = {
        "name": "Cox proportional-hazards model",
        "n": n_subj, "events": n_events, "covariates": covs,
        "proportional_hazards": {"test": "cox.zph global", "p_value": ph_p},
        "hazard_ratios": [{"term": tt, "HR": float(h_), "ci_low": float(l_),
                           "ci_high": float(u_), "p_value": float(p_), "note": n_}
                          for tt, h_, l_, u_, p_, n_ in zip(terms, hr, lo, hi, pv,
                                                            stats_df["note"])],
    }
    methods = _methods(stats_df, ph_p, covs, n_subj, n_events, sep)
    return {
        "fig": fp["fig"], "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "label": "hazard_ratios",
        "width": fp["width"], "height": fp["height"],
        "build_script": functools.partial(_script, spec, t, e, covs, fp["ticks"],
                                          fp["label_in"], fp["width"], fp["height"]),
    }


def _methods(stats_df, ph_p, covs, n, events, sep):
    """Write the methods paragraph, with the same sentences as the R engine.

    Args:
        stats_df: The stats table.
        ph_p: The P value of the global cox.zph test, or NaN.
        covs: The covariates.
        n: The number of subjects.
        events: The number of events.
        sep: The separation direction of each term, 1, -1 or 0.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    sig = stats_df[stats_df["p_value"] < 0.05]
    if len(sig):
        parts = "; ".join(f"{r.term} (HR = {r.HR:.2f}, 95% CI {r.ci_low:.2f}-{r.ci_high:.2f}, {p_phrase(r.p_value)})"
                          for r in sig.itertuples())
        sig_txt = f" Significant predictors: {parts}."
    else:
        sig_txt = " No individual predictor reached significance."
    ph_txt = "" if np.isnan(ph_p) else (
        f" The proportional-hazards assumption was {'violated' if ph_p < 0.05 else 'not violated'}"
        f" (global cox.zph {p_phrase(ph_p)}).")
    return (
        f"Hazard ratios were estimated with a Cox proportional-hazards model ({n} subjects, {events} events) "
        f"including {', '.join(covs)}." + sig_txt + separation_text(list(stats_df["term"]), sep)
        + ph_txt +
        f" Analyses were performed in Python {platform.python_version()} with lifelines {v('lifelines')}"
        f", and the forest plot was drawn with matplotlib {v('matplotlib')}."
    )


def _script(spec, t, e, covs, ticks, lab_in, w, h, in_name, fig_stub):
    """Write the standalone script that refits the model and redraws the figure.

    Args:
        spec: The figure specification.
        t: The time column.
        e: The event column.
        covs: The covariate columns.
        ticks: The tick positions of the hazard-ratio axis.
        lab_in: The width of the label column in inches.
        w: The figure width in inches.
        h: The figure height in inches.
        in_name: The bundled input file name.
        fig_stub: The figure file stem.

    Returns:
        str: The script text.
    """
    ap = spec["appearance"]
    colour = palette(1, ap.get("palette"))[0]
    xlab = ap.get("x_label") or "Hazard ratio (95% CI, log scale)"
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from matplotlib import transforms
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter

df = pd.read_csv("{in_name}")
covs = {covs!r}
design = df[["{t}", "{e}"]].copy()
rows = []  # (label, term); term None marks the reference level of a factor
for c in covs:
    s = df[c]
    if s.dtype == object or not np.issubdtype(s.dtype, np.number):
        cats = sorted(s.astype(str).unique())
        rows.append(("%s: %s" % (c, cats[0]), None))
        for lvl in cats[1:]:
            design["%s%s" % (c, lvl)] = (s.astype(str) == lvl).astype(int)
            rows.append(("%s: %s" % (c, lvl), "%s%s" % (c, lvl)))
    else:
        design[c] = pd.to_numeric(s, errors="coerce")
        rows.append((c, c))
cph = CoxPHFitter().fit(design, "{t}", "{e}", fit_options={{"precision": 1e-11}})
s = cph.summary
print(s)

k = len(rows); yy = np.arange(k, 0, -1).astype(float)
ticks = {[float(x) for x in ticks]!r}
w, h, lab_in = {w:.4f}, {h:.4f}, {lab_in:.4f}
fig = plt.figure(figsize=(w, h))
ax = fig.add_axes([lab_in / w, 0.85 / h, 2.85 / w, (h - 1.3) / h])
axt = fig.add_axes([(lab_in + 3.0) / w, 0.85 / h, 2.3 / w, (h - 1.3) / h])
ax.axvline(1.0, ls=":", color="0.4", lw=0.9)
trans = transforms.blended_transform_factory(fig.transFigure, ax.transData)
for y_, (label, term) in zip(yy, rows):
    fig.text(0.15 / w, y_, label, transform=trans, ha="left", va="center", fontsize=10.5)
    if term is None:
        axt.text(0, y_, "reference", ha="left", va="center", fontsize=10)
        continue
    hr_, lo_, hi_, p_ = s.loc[term, ["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]]
    ax.hlines(y_, lo_, hi_, color="0.15", lw=1.0)
    ax.scatter([hr_], [y_], marker="s", s=48, color="{colour}", zorder=3)
    star = "***" if p_ < 0.001 else "**" if p_ < 0.01 else "*" if p_ < 0.05 else ""
    p_txt = ("< 0.001" if p_ < 0.001 else "%.2g" % p_) + (" " + star if star else "")
    axt.text(0, y_, "%.2f (%.2f\\u2013%.2f)" % (hr_, lo_, hi_), ha="left", va="center", fontsize=10)
    axt.text(0.66, y_, p_txt, ha="left", va="center", fontsize=10)
axt.text(0, k + 1, "HR (95% CI)", ha="left", va="center", fontsize=10, fontweight="bold")
axt.text(0.66, k + 1, "P", ha="left", va="center", fontsize=10, fontweight="bold")
ax.set_xscale("log")
ax.xaxis.set_major_locator(FixedLocator(ticks))
ax.xaxis.set_major_formatter(FuncFormatter("{{:g}}".format))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_ylim(0.4, k + 1.6); axt.set_ylim(0.4, k + 1.6); axt.set_xlim(0, 1); axt.axis("off")
ax.set_yticks([]); ax.set_xlabel({xlab!r}, fontsize=13)
ax.set_title("%d subjects, %d events" % (len(df), int(df["{e}"].sum())), fontsize=9.5,
             color="#4d4d4d", loc="left")
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
fig.savefig("{fig_stub}.pdf", bbox_inches="tight")
print("Reproduced figure.")
'''
