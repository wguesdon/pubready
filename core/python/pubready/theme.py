"""The pubready house style for matplotlib, matched to the R theme."""
import matplotlib
matplotlib.use("Agg")  # headless container
import matplotlib.pyplot as plt

PALETTE = ["#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A"]


def palette(n, custom=None):
    cols = custom if custom else PALETTE
    return [cols[i % len(cols)] for i in range(n)]


def apply_pub_style(ax, theme="pubready_house"):
    if theme == "prism":
        return apply_prism_style(ax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.8)
        ax.spines[s].set_color("black")
    ax.tick_params(colors="black", labelsize=11)
    for lbl in (ax.xaxis.label, ax.yaxis.label):
        lbl.set_color("black")
        lbl.set_fontsize(13)
    return ax


def apply_prism_style(ax):
    """A GraphPad Prism look: thick black axes, outward bold ticks, larger
    labels. Matches the R ggprism::theme_prism variant."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.3)
        ax.spines[s].set_color("black")
    ax.tick_params(colors="black", labelsize=12, width=1.3, length=6,
                   direction="out")
    for lbl in (ax.xaxis.label, ax.yaxis.label):
        lbl.set_color("black")
        lbl.set_fontsize(14)
        lbl.set_fontweight("bold")
    for lbl in ax.get_xticklabels():
        lbl.set_fontweight("bold")
    return ax


def new_fig(w=3.8, h=4.0):
    fig, ax = plt.subplots(figsize=(w, h))
    return fig, ax


def _raincloud(ax, data, x, y, order, pal, show_points=True):
    import numpy as np
    for i, g in enumerate(order):
        vals = data.loc[data[x] == g, y].dropna().to_numpy()
        if len(vals) < 2:
            continue
        v = ax.violinplot([vals], positions=[i], showextrema=False, widths=0.6)
        for body in v["bodies"]:
            verts = body.get_paths()[0].vertices
            verts[:, 0] = np.clip(verts[:, 0], i, np.inf)  # keep the right half
            body.set_facecolor(pal[i]); body.set_edgecolor("none"); body.set_alpha(0.85)
        bp = ax.boxplot([vals], positions=[i], widths=0.1, patch_artist=True,
                        showfliers=False, manage_ticks=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(pal[i]); patch.set_alpha(0.9)
        for med in bp["medians"]:
            med.set_color("black")
        if show_points:
            jit = np.random.RandomState(i).uniform(-0.18, -0.06, size=len(vals))
            ax.scatter(np.full(len(vals), i) + jit, vals, s=10, color=pal[i],
                       alpha=0.6, edgecolors="black", linewidths=0.2, zorder=3)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order)


def add_dist_geom(ax, data, x, y, order, pal, geom="box", show_points=True):
    """Draw the distribution geom for the comparison recipes: box | violin | bar |
    raincloud, plus optional jittered points."""
    import numpy as np
    import seaborn as sns
    if geom == "violin":
        sns.violinplot(data=data, x=x, y=y, hue=x, order=order, hue_order=order,
                       palette=pal, ax=ax,
                       legend=False, cut=0, inner=None, alpha=0.9)
    elif geom == "bar":
        means = [data.loc[data[x] == g, y].mean() for g in order]
        sems = [data.loc[data[x] == g, y].sem() for g in order]
        ax.bar(range(len(order)), means, yerr=sems, width=0.6, color=pal,
               edgecolor="black", linewidth=0.5, alpha=0.9, capsize=4)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
    elif geom == "raincloud":
        _raincloud(ax, data, x, y, order, pal, show_points)
        return
    else:
        sns.boxplot(data=data, x=x, y=y, hue=x, order=order, hue_order=order,
                    palette=pal, ax=ax,
                    legend=False, width=0.6, fliersize=0)
    if show_points:
        sns.stripplot(data=data, x=x, y=y, hue=x, order=order, hue_order=order,
                      palette=pal, ax=ax,
                      legend=False, size=4, alpha=0.7, edgecolor="black",
                      linewidth=0.3, jitter=0.12)
