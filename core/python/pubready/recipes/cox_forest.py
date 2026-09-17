"""Recipe: Cox proportional-hazards model with a forest plot of hazard ratios
(lifelines). Mirrors the R cox_forest: HR with 95% CI and p per covariate, and a
proportional-hazards check."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from ..theme import apply_pub_style
from ..util import format_p


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
    ev = pd.to_numeric(df[e], errors="coerce")
    if not ev.dropna().isin([0, 1]).all():
        raise ValueError(f"event column '{e}' must be coded 1 = event, 0 = censored")
    df[e] = ev
    n0 = len(df)
    df = df.dropna(subset=[t, e] + covs)
    df = df[df[t] >= 0].reset_index(drop=True)
    clean_steps = (f"Dropped {n0 - len(df)} row(s) with missing/invalid values ({n0} -> {len(df)} rows)."
                   if len(df) < n0 else "No cleaning applied; input used as-is.")

    # Build a numeric design matrix; one-hot encode categoricals (drop first = reference).
    design = df[[t, e]].copy()
    for c in covs:
        s = df[c]
        if s.dtype == object or str(s.dtype).startswith("category") or not np.issubdtype(s.dtype, np.number):
            cats = sorted(s.astype(str).unique())
            for lvl in cats[1:]:
                design[f"{c}[{lvl}]"] = (s.astype(str) == lvl).astype(int)
        else:
            design[c] = pd.to_numeric(s, errors="coerce")

    cph = CoxPHFitter()
    cph.fit(design, duration_col=t, event_col=e)
    summ = cph.summary
    terms = list(summ.index)
    hr = summ["exp(coef)"].to_numpy()
    lo = summ["exp(coef) lower 95%"].to_numpy()
    hi = summ["exp(coef) upper 95%"].to_numpy()
    pv = summ["p"].to_numpy()

    ph_p = np.nan
    try:
        from lifelines.statistics import proportional_hazard_test
        ph = proportional_hazard_test(cph, design, time_transform="rank")
        ph_p = float(ph.summary["p"].min())
    except Exception:
        pass

    n = len(terms)
    w, h = 7.4, max(3.6, 0.5 * n + 2.0)
    fig, ax = plt.subplots(figsize=(w, h))
    yy = np.arange(n)[::-1]
    ax.errorbar(hr, yy, xerr=[hr - lo, hi - hr], fmt="s", color="black",
                capsize=3, ms=6, lw=1.0)
    ax.axvline(1.0, ls=":", color="grey")
    ax.set_xscale("log")
    ax.set_yticks(yy)
    ax.set_yticklabels(terms)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Hazard ratio (95% CI)")
    for y_, h_, l_, u_, p_ in zip(yy, hr, lo, hi, pv):
        star = " *" if p_ < 0.05 else ""
        ax.annotate(f"{h_:.2f} ({l_:.2f}–{u_:.2f})   P {format_p(p_) if p_ < 0.001 else '= ' + format(p_, '.2g')}{star}",
                    xy=(1.02, y_), xycoords=("axes fraction", "data"), va="center", fontsize=8)
    ax.set_title("Hazard ratio", fontweight="bold")
    apply_pub_style(ax)
    fig.tight_layout()

    stats_df = pd.DataFrame({"term": terms, "HR": hr, "ci_low": lo, "ci_high": hi,
                             "z": summ["z"].to_numpy(), "p_value": pv})
    test_meta = {
        "name": "Cox proportional-hazards model",
        "n": int(len(df)), "events": int(df[e].sum()), "covariates": covs,
        "proportional_hazards": {"test": "proportional_hazard_test (min p)", "p_value": ph_p},
        "hazard_ratios": [{"term": tt, "HR": float(h_), "ci_low": float(l_),
                           "ci_high": float(u_), "p_value": float(p_)}
                          for tt, h_, l_, u_, p_ in zip(terms, hr, lo, hi, pv)],
    }
    methods = _methods(stats_df, ph_p, covs, len(df), int(df[e].sum()))
    return {
        "fig": fig, "stats": stats_df, "test_meta": test_meta, "methods": methods,
        "df_used": df, "clean_steps": clean_steps, "label": "hazard_ratios",
        "width": w, "height": h,
        "build_script": lambda in_name, fig_stub: _script(spec, t, e, covs, in_name, fig_stub),
    }


def _methods(stats_df, ph_p, covs, n, events):
    from importlib.metadata import version as v
    import platform
    sig = stats_df[stats_df["p_value"] < 0.05]
    if len(sig):
        parts = "; ".join(f"{r.term} (HR = {r.HR:.2f}, 95% CI {r.ci_low:.2f}-{r.ci_high:.2f}, P = {format_p(r.p_value)})"
                          for r in sig.itertuples())
        sig_txt = f" Significant predictors: {parts}."
    else:
        sig_txt = " No individual predictor reached significance."
    ph_txt = "" if np.isnan(ph_p) else (
        f" The proportional-hazards assumption was {'violated' if ph_p < 0.05 else 'not violated'}"
        f" (minimum covariate cox.zph P = {format_p(ph_p)}).")
    return (
        f"Hazard ratios were estimated with a Cox proportional-hazards model ({n} subjects, {events} events) "
        f"including {', '.join(covs)}." + sig_txt + ph_txt +
        f" Analyses were performed in Python {platform.python_version()} with lifelines {v('lifelines')}."
    )


def _script(spec, t, e, covs, in_name, fig_stub):
    return f'''#!/usr/bin/env python3
# Standalone reproduction. Run inside the pinned container.
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from lifelines import CoxPHFitter

df = pd.read_csv("{in_name}")
covs = {covs}
design = df[["{t}", "{e}"]].copy()
for c in covs:
    s = df[c]
    if s.dtype == object or not np.issubdtype(s.dtype, np.number):
        for lvl in sorted(s.astype(str).unique())[1:]:
            design["%s[%s]" % (c, lvl)] = (s.astype(str) == lvl).astype(int)
    else:
        design[c] = pd.to_numeric(s, errors="coerce")
cph = CoxPHFitter().fit(design, "{t}", "{e}")
s = cph.summary
terms = list(s.index); hr = s["exp(coef)"].to_numpy()
lo = s["exp(coef) lower 95%"].to_numpy(); hi = s["exp(coef) upper 95%"].to_numpy(); pv = s["p"].to_numpy()
n = len(terms); yy = np.arange(n)[::-1]
fig, ax = plt.subplots(figsize=(7.4, max(3.6, 0.5 * n + 2.0)))
ax.errorbar(hr, yy, xerr=[hr - lo, hi - hr], fmt="s", color="black", capsize=3, ms=6, lw=1.0)
ax.axvline(1.0, ls=":", color="grey"); ax.set_xscale("log")
ax.set_yticks(yy); ax.set_yticklabels(terms); ax.set_ylim(-0.6, n - 0.4)
ax.set_xlabel("Hazard ratio (95% CI)"); ax.set_title("Hazard ratio", fontweight="bold")
for y_, h_, l_, u_, p_ in zip(yy, hr, lo, hi, pv):
    ax.annotate("%.2f (%.2f-%.2f)  P = %.2g" % (h_, l_, u_, p_), xy=(1.02, y_), xycoords=("axes fraction", "data"), va="center", fontsize=8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced figure.")
'''
