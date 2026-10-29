"""Recipe: Sankey (alluvial) diagram of patient state transitions across visits.

Mirrors the R recipe. Input is one row per patient visit: --id, --x (the visit)
and --y (the state). Visits are sorted numerically when every visit is a number,
else in code-point order; states numerically when every state is a number, else
by first appearance in the data. A patient missing at a visit takes no part in
the transitions into or out of that visit.

The layout (stacked nodes with small gaps, cubic ribbons coloured by the state
of origin, labels where they fit) restates the arithmetic of the R engine. With
exactly two visits that carry the same set of states, Bowker's test of symmetry
(McNemar's test for two states, no continuity correction) runs on the patients
present at both visits, with statsmodels' SquareTable.
"""
import functools
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy import stats
from statsmodels.stats.contingency_tables import SquareTable

from ..bundle import slugify
from ..theme import apply_pub_style, palette

#: Half width of a node, in visit units.
HW = 0.045
#: Gap between nodes, as a fraction of the largest visit.
GAP_FRAC = 0.03
#: Labels in one column stay at least this far apart, in inches.
LABEL_IN = 0.17


def fmt_value(v):
    """Format one value as R's as.character writes it.

    Args:
        v: The value.

    Returns:
        str or float: The text, or NaN for a missing value.
    """
    if pd.isna(v):
        return np.nan
    if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool):
        return f"{float(v):.15g}"
    return str(v).strip()


def labels_of(s):
    """Turn a column into strings that match the R engine.

    Args:
        s: The column.

    Returns:
        pandas.Series: The labels, with missing values kept as NaN.
    """
    return s.map(fmt_value)


def visit_order(v):
    """Visit order: numeric when every value is a number, else code-point order.

    Args:
        v: The visit labels.

    Returns:
        list: The visits in order.
    """
    u = list(pd.unique(v))
    num = pd.to_numeric(pd.Series(u), errors="coerce")
    if num.notna().all():
        return [u[i] for i in np.lexsort((np.array(u, dtype=object).astype(str), num.to_numpy()))]
    return sorted(u)


def state_order(s):
    """State order: numeric when every value is a number, else first appearance.

    Args:
        s: The state labels, in row order.

    Returns:
        list: The states in order.
    """
    u = list(pd.unique(s))
    num = pd.to_numeric(pd.Series(u), errors="coerce")
    if num.notna().all():
        return [u[i] for i in np.lexsort((np.array(u, dtype=object).astype(str), num.to_numpy()))]
    return u


def paired(d, a, b):
    """The patients present at two visits, with their state at each.

    Args:
        d: The cleaned table.
        a: The first visit.
        b: The second visit.

    Returns:
        pandas.DataFrame: One row per patient: id, state_from, state_to.
    """
    x = d.loc[d["visit"] == a, ["id", "state"]]
    y = d.loc[d["visit"] == b, ["id", "state"]]
    return x.merge(y, on="id", suffixes=("_from", "_to"))


def transitions(d, visits, states):
    """The transitions table between consecutive visits.

    Args:
        d: The cleaned table.
        visits: The visits in order.
        states: The states in order.

    Returns:
        pandas.DataFrame: from_visit, to_visit, from_state, to_state and n.
    """
    rows = []
    for a, b in zip(visits[:-1], visits[1:]):
        m = paired(d, a, b)
        for s in states:
            for t in states:
                n = int(((m["state_from"] == s) & (m["state_to"] == t)).sum())
                if n > 0:
                    rows.append({"from_visit": a, "to_visit": b, "from_state": s,
                                 "to_state": t, "n": n})
    return pd.DataFrame(rows, columns=["from_visit", "to_visit", "from_state", "to_state", "n"])


def symmetry(tab):
    """Bowker's test of the square table (McNemar's test for two states).

    Args:
        tab: The square table of counts, rows the first visit.

    Returns:
        dict or None: The statistic, df, p-value and count of empty pairs, or
        None when no patient changed state.
    """
    k = tab.shape[0]
    iu = np.triu_indices(k, 1)
    up, lo = tab[iu], tab.T[iu]
    keep = (up + lo) > 0
    if keep.sum() == 0:
        return None
    if keep.all():
        res = SquareTable(tab, shift_zeros=False).symmetry(method="bowker")
        return {"statistic": float(res.statistic), "df": float(res.df),
                "p_value": float(res.pvalue), "empty_pairs": 0}
    stat = float(np.sum((up[keep] - lo[keep]) ** 2 / (up[keep] + lo[keep])))
    df = float(keep.sum())
    return {"statistic": stat, "df": df, "p_value": float(stats.chi2.sf(stat, df)),
            "empty_pairs": int((~keep).sum())}


def run_test(d, visits, states):
    """Bowker's or McNemar's test when the design allows it.

    Args:
        d: The cleaned table.
        visits: The visits in order.
        states: The states in order.

    Returns:
        dict: The test, or the reason no test was run.
    """
    if len(visits) != 2:
        return {"reason": "more than two visits"}
    m = paired(d, visits[0], visits[1])
    s1 = [s for s in states if s in set(m["state_from"])]
    s2 = [s for s in states if s in set(m["state_to"])]
    if set(s1) != set(s2):
        return {"reason": "the two visits do not carry the same states"}
    if len(s1) < 2:
        return {"reason": "fewer than two states"}
    tab = np.array([[int(((m["state_from"] == s) & (m["state_to"] == t)).sum()) for t in s1]
                    for s in s1], dtype=float)
    res = symmetry(tab)
    if res is None:
        return {"reason": "no patient changed state"}
    res["name"] = "McNemar's test" if len(s1) == 2 else "Bowker's test of symmetry"
    res["n_paired"] = len(m)
    res["states"] = s1
    res["table"] = tab.astype(int)
    return res


def layout_nodes(d, visits, states):
    """Node positions: states stacked from the top, each column centred.

    Args:
        d: The cleaned table.
        visits: The visits in order.
        states: The states in order.

    Returns:
        tuple: The node table, the tallest column height and the gap.
    """
    counts = np.array([[int(((d["visit"] == v) & (d["state"] == s)).sum()) for s in states]
                       for v in visits])
    nmax = counts.sum(axis=1).max()
    gap = GAP_FRAC * nmax
    heights = counts.sum(axis=1) + gap * ((counts > 0).sum(axis=1) - 1)
    hmax = heights.max()
    rows = []
    for vi, v in enumerate(visits):
        top = hmax - (hmax - heights[vi]) / 2
        for si, s in enumerate(states):
            n = counts[vi, si]
            if n == 0:
                continue
            rows.append({"visit": v, "vi": vi, "state": s, "n": int(n), "ytop": top,
                         "ybot": top - n})
            top = top - n - gap
    return pd.DataFrame(rows), float(hmax), float(gap)


def ribbon_polygons(trans, nodes, visits, states):
    """Ribbon polygons between consecutive visits.

    Outflows leave a node from its top in the order of the target state, and
    inflows enter a node from its top in the order of the source state.

    Args:
        trans: The transitions table.
        nodes: The node table.
        visits: The visits in order.
        states: The states in order.

    Returns:
        list: One (state of origin, x array, y array) per transition.
    """
    key = {(r.visit, r.state): r for r in nodes.itertuples()}
    out_off = {k: 0.0 for k in key}
    in_off = {k: 0.0 for k in key}
    tr = trans.reset_index(drop=True)
    vrank = {v: i for i, v in enumerate(visits)}
    srank = {s: i for i, s in enumerate(states)}
    y0 = np.zeros(len(tr))
    y1 = np.zeros(len(tr))
    for i in sorted(range(len(tr)), key=functools.partial(_rank_out, tr, vrank, srank)):
        k = (tr.at[i, "from_visit"], tr.at[i, "from_state"])
        y0[i] = key[k].ytop - out_off[k]
        out_off[k] += tr.at[i, "n"]
    for i in sorted(range(len(tr)), key=functools.partial(_rank_in, tr, vrank, srank)):
        k = (tr.at[i, "to_visit"], tr.at[i, "to_state"])
        y1[i] = key[k].ytop - in_off[k]
        in_off[k] += tr.at[i, "n"]
    u = np.linspace(0, 1, 40)
    s = 3 * u ** 2 - 2 * u ** 3
    polys = []
    for i in range(len(tr)):
        x0 = key[(tr.at[i, "from_visit"], tr.at[i, "from_state"])].vi + HW
        x1 = key[(tr.at[i, "to_visit"], tr.at[i, "to_state"])].vi - HW
        n = tr.at[i, "n"]
        xs = x0 + (x1 - x0) * u
        top = y0[i] + (y1[i] - y0[i]) * s
        bot = (y0[i] - n) + ((y1[i] - n) - (y0[i] - n)) * s
        polys.append((tr.at[i, "from_state"], np.concatenate([xs, xs[::-1]]),
                      np.concatenate([top, bot[::-1]])))
    return polys


def _rank_out(tr, vrank, srank, i):
    """Sort key of a transition among the outflows of its node.

    Args:
        tr: The transitions table.
        vrank: The rank of each visit.
        srank: The rank of each state.
        i: The row.

    Returns:
        tuple: The key.
    """
    return (vrank[tr.at[i, "from_visit"]], srank[tr.at[i, "from_state"]],
            srank[tr.at[i, "to_state"]])


def _rank_in(tr, vrank, srank, i):
    """Sort key of a transition among the inflows of its node.

    Args:
        tr: The transitions table.
        vrank: The rank of each visit.
        srank: The rank of each state.
        i: The row.

    Returns:
        tuple: The key.
    """
    return (vrank[tr.at[i, "to_visit"]], srank[tr.at[i, "to_state"]],
            srank[tr.at[i, "from_state"]])


def figure_dims(nodes, visits, states, ncol_leg):
    """The figure size and the label room, matched to the R engine.

    Args:
        nodes: The node table.
        visits: The visits in order.
        states: The states in order.
        ncol_leg: The legend columns.

    Returns:
        dict: The width, the height, the panel height, the x limits and the
        inches per visit.
    """
    chars = max(len(f"{s} ({n})") for s, n in zip(nodes["state"], nodes["n"]))
    label_in = 0.075 * chars + 0.25
    nv = len(visits)
    # Visits sit far enough apart for their axis labels, and for the node labels
    # of a middle visit, which sit between two columns.
    unit_in = max(2.0, 0.09 * max(len(v) for v in visits) + 0.4, label_in + 0.6 if nv > 2 else 0)
    legend_rows = math.ceil(len(states) / ncol_leg)
    legend_h = 0.3 + 0.25 * legend_rows
    height = max(4.2, 2.4 + 0.32 * len(states)) + legend_h
    pad = HW + 0.03 + label_in / unit_in
    return {"width": unit_in * (nv - 1) + 2 * (HW * unit_in + 0.03 * unit_in + label_in) + 0.6,
            "height": height, "panel_h": height - 1.5 - legend_h,
            "xlim": (-pad, nv - 1 + pad), "unit_in": unit_in}


def legend_ncol(states, width):
    """The legend columns, matched to the R engine.

    Args:
        states: The states.
        width: The figure width in inches.

    Returns:
        int: The count of columns.
    """
    per = 0.085 * max(len(s) for s in states) + 0.6
    ncol = min(len(states), max(1, math.floor((width - 0.5) / per)))
    return int(math.ceil(len(states) / math.ceil(len(states) / ncol)))


def subtitle(test, n_patients, n_visits, width):
    """The subtitle with the test result.

    A subtitle wider than the figure breaks after the patient count.

    Args:
        test: The test, or the reason none was run.
        n_patients: The count of patients.
        n_visits: The count of visits.
        width: The figure width in inches.

    Returns:
        str: The subtitle.
    """
    base = f"{n_patients} patients, {n_visits} visits"
    if "statistic" not in test:
        return base
    p = "P < 0.001" if test["p_value"] < 0.001 else f"P = {test['p_value']:.2g}"
    tst = (f"{test['name']} (n = {test['n_paired']}): chi-square({int(test['df'])}) = "
           f"{test['statistic']:.2f}, {p}")
    one = f"{base}; {tst}"
    return f"{base}\n{tst}" if 0.07 * len(one) + 0.4 > width else one


def recipe_sankey(df, spec):
    """Draw the state transitions of patients across visits.

    Args:
        df: The tidy table, one row per patient visit.
        spec: The figure specification.

    Returns:
        dict: The bundle payload.

    Raises:
        ValueError: If a column is missing, if a patient has two rows at one
            visit, or if there are fewer than two visits.
    """
    dat = spec["data"]
    idc, x, y = dat.get("id"), dat.get("x"), dat.get("y")
    if not idc:
        raise ValueError("sankey needs --id (the patient column)")
    if not x or not y:
        raise ValueError("sankey needs --x (the visit) and --y (the state)")
    for c in (idc, x, y):
        if c not in df.columns:
            raise ValueError(f"column '{c}' not found in data")

    n0 = len(df)
    d = pd.DataFrame({"id": labels_of(df[idc]), "visit": labels_of(df[x]),
                      "state": labels_of(df[y])})
    keep = d.notna().all(axis=1) & (d["state"].astype(str) != "")
    d = d[keep].reset_index(drop=True)
    if d.duplicated(["id", "visit"]).any():
        raise ValueError("sankey needs one row per patient and visit; some patients have two "
                         "rows at the same visit")
    visits = visit_order(d["visit"])
    if len(visits) < 2:
        raise ValueError("sankey needs at least 2 visits")
    states = state_order(d["state"])
    trans = transitions(d, visits, states)
    test = run_test(d, visits, states)
    n_pat = int(d["id"].nunique())
    n_partial = int((d.groupby("id").size() < len(visits)).sum())
    steps = []
    if int((~keep).sum()):
        steps.append(f"Dropped {int((~keep).sum())} row(s) with a missing patient, visit or state.")
    steps.append(f"Kept {len(d)} row(s) of {n0}: {n_pat} patient(s), {len(visits)} visit(s) "
                 f"({', '.join(visits)}), {len(states)} state(s) ({', '.join(states)}).")
    steps.append(f"{n_partial} patient(s) are missing at one visit or more and take no part in "
                 "the transitions into or out of that visit.")

    nodes, hmax, _ = layout_nodes(d, visits, states)
    width0 = figure_dims(nodes, visits, states, len(states))["width"]
    ncol_leg = legend_ncol(states, width0)
    dims = figure_dims(nodes, visits, states, ncol_leg)
    sub = subtitle(test, n_pat, len(visits), dims["width"])
    if "\n" in sub:
        dims["height"] += 0.2
    ap = spec["appearance"]
    pal = dict(zip(states, palette(len(states), ap.get("palette"))))
    fig = draw_figure(trans, nodes, hmax, visits, states, pal, dims, ap, x, y, ncol_leg,
                      sub)

    has_test = "statistic" in test
    test_meta = {"name": test["name"] if has_test else "Sankey diagram (descriptive)",
                 "n_patients": n_pat, "n_visits": len(visits), "visits": visits,
                 "states": states, "n_transitions": int(trans["n"].sum())}
    if has_test:
        test_meta.update({
            "statistic": test["statistic"], "df": test["df"], "p_value": test["p_value"],
            "n_paired": test["n_paired"], "empty_pairs": test["empty_pairs"],
            "table": [dict(zip(test["states"], map(int, row))) for row in test["table"]]})
    else:
        test_meta["no_test_reason"] = test["reason"]
    return {
        "fig": fig, "stats": trans, "test_meta": test_meta,
        "methods": _methods(test, trans, visits, n_pat, idc, x, y),
        "df_used": d, "clean_steps": steps, "label": slugify(y),
        "width": dims["width"], "height": dims["height"],
        "build_script": functools.partial(_script, spec, visits, dims),
    }


def draw_figure(trans, nodes, hmax, visits, states, pal, dims, ap, x, y, ncol_leg, sub):
    """Draw the nodes, the ribbons and the labels.

    Args:
        trans: The transitions table.
        nodes: The node table.
        hmax: The height of the tallest column.
        visits: The visits in order.
        states: The states in order.
        pal: The colour of each state.
        dims: The figure size.
        ap: The appearance options.
        x: The visit column.
        y: The state column.
        ncol_leg: The legend columns.
        sub: The subtitle.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    fig = plt.figure(figsize=(dims["width"], dims["height"]), layout="constrained")
    ax = fig.add_subplot()
    for state, xs, ys in ribbon_polygons(trans, nodes, visits, states):
        ax.fill(xs, ys, color=pal[state], alpha=0.45, lw=0)
    show = label_fit(nodes, dims["panel_h"] / (1.04 * hmax))
    box = {"boxstyle": "square,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.75}
    for r, shown in zip(nodes.itertuples(), show):
        ax.fill([r.vi - HW, r.vi + HW, r.vi + HW, r.vi - HW], [r.ybot, r.ybot, r.ytop, r.ytop],
                color=pal[r.state], lw=0)
        if not shown:
            continue
        ymid = (r.ytop + r.ybot) / 2
        if r.vi == 0:
            ax.text(r.vi - HW - 0.03, ymid, f"{r.state} ({r.n})", ha="right", va="center",
                    fontsize=9, bbox=box, clip_on=False)
        else:
            ax.text(r.vi + HW + 0.03, ymid, f"{r.state} ({r.n})", ha="left", va="center",
                    fontsize=9, bbox=box, clip_on=False)
    ax.set_xlim(*dims["xlim"])
    ax.set_ylim(-0.02 * hmax, 1.02 * hmax)
    ax.set_xticks(range(len(visits)))
    ax.set_xticklabels(visits)
    ax.set_yticks([])
    apply_pub_style(ax, ap.get("theme"))
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="x", length=0, labelsize=11, pad=4)
    ax.set_xlabel(ap.get("x_label") or x)
    ax.set_title(sub, fontsize=10, color="#4d4d4d")
    handles = [Patch(facecolor=pal[s], edgecolor="none") for s in states]
    fig.legend(handles, states, title=y, loc="outside lower center", ncol=ncol_leg,
               frameon=False, fontsize=10, title_fontsize=11, handlelength=1.1,
               handleheight=1.1)
    if ap.get("title"):
        fig.suptitle(ap["title"], fontweight="bold")
    return fig


def _methods(test, trans, visits, n_pat, idc, x, y):
    """Write the methods paragraph.

    Args:
        test: The test, or the reason none was run.
        trans: The transitions table.
        visits: The visits in order.
        n_pat: The count of patients.
        idc: The patient column.
        x: The visit column.
        y: The state column.

    Returns:
        str: The paragraph.
    """
    from importlib.metadata import version as v
    import platform
    txt = (f"The state {y} of {n_pat} patients ({idc}) was followed across {len(visits)} visits "
           f"of {x} ({', '.join(visits)}). Visits are ordered numerically when every visit is a "
           "number, else alphabetically, and states numerically when every state is a number, "
           "else by their first appearance in the data. The diagram shows "
           f"{int(trans['n'].sum())} transitions between consecutive visits; a patient missing at "
           "a visit takes no part in the transitions into or out of it, and the transitions table "
           "gives every count. ")
    if "statistic" in test:
        p = "P < 0.001" if test["p_value"] < 0.001 else f"P = {test['p_value']:.2g}"
        txt += (f"{test['name']}, without continuity correction, compared the {test['n_paired']} "
                f"patients present at both visits (chi-square({int(test['df'])}) = "
                f"{test['statistic']:.2f}, {p}).")
        if test["empty_pairs"] > 0:
            txt += (f" {test['empty_pairs']} pair(s) of cells held no patient and were left out "
                    "of the statistic and its degrees of freedom.")
    else:
        txt += f"The diagram is descriptive, and no test was applied because {test['reason']}."
    return (f"{txt} Drawn in Python {platform.python_version()} with matplotlib "
            f"{v('matplotlib')}; the test uses statsmodels {v('statsmodels')}.")


def _script(spec, visits, dims, in_name, fig_stub):
    """Write the standalone script that recomputes the transitions and the test.

    Args:
        spec: The figure specification.
        visits: The visits in order.
        dims: The figure size.
        in_name: The name of the input copy inside the bundle.
        fig_stub: The name of the figure file, without its extension.

    Returns:
        str: The script.
    """
    dat = spec["data"]
    idc, x, y = dat["id"], dat["x"], dat["y"]
    return f'''#!/usr/bin/env python3
# Standalone reproduction of the transition counts and the test. Run inside the
# pinned pubready container.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.contingency_tables import SquareTable

raw = pd.read_csv("{in_name}")
d = raw[["{idc}", "{x}", "{y}"]].dropna().astype(str)
d.columns = ["id", "visit", "state"]
d["visit"] = d["visit"].str.replace(r"\\.0$", "", regex=True)
visits = {visits!r}
fig, ax = plt.subplots(figsize=({dims["width"]:.2f}, {dims["height"]:.2f}))
for k, (a, b) in enumerate(zip(visits[:-1], visits[1:])):
    m = d[d["visit"] == a].merge(d[d["visit"] == b], on="id")
    tab = pd.crosstab(m["state_x"], m["state_y"])
    print(tab)
    if len(visits) == 2 and list(tab.index) == list(tab.columns):
        print(SquareTable(tab.to_numpy(), shift_zeros=False).symmetry(method="bowker"))
    tab.sum(axis=1).plot.bar(ax=ax, position=k, width=0.3)
ax.set_xlabel("{y}"); ax.set_ylabel("Patients")
fig.tight_layout(); fig.savefig("{fig_stub}.pdf")
print("Reproduced the transition counts.")
'''


def label_fit(nodes, in_per_unit):
    """Which nodes show their label, matched to the R engine.

    From the top down, a node shows its label when the label centre is at least
    LABEL_IN below the last label shown in the same column, so labels never
    overlap.

    Args:
        nodes: The node table, in column and stacking order.
        in_per_unit: Inches per patient on the count axis.

    Returns:
        list: One flag per node.
    """
    show = [False] * len(nodes)
    last = {}
    for i, r in enumerate(nodes.itertuples()):
        y_in = (r.ytop + r.ybot) / 2 * in_per_unit
        if last.get(r.vi, math.inf) - y_in >= LABEL_IN:
            show[i] = True
            last[r.vi] = y_in
    return show
