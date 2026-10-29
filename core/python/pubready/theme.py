"""The pubready house style for matplotlib, matched to the R theme."""
import matplotlib
matplotlib.use("Agg")  # headless container
import matplotlib.pyplot as plt

# The first six colours are the house palette. The fourteen after them keep a
# figure with up to twenty groups distinct; palette() cycles past twenty.
PALETTE = ["#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A",
           "#56B4E9", "#E78AC3", "#8C6D31", "#66C2A5", "#F28E2B", "#A6D854",
           "#1F3A5F", "#B07AA1", "#9C755F", "#BAB0AC", "#17BECF", "#FF9DA7",
           "#4E7A27", "#D4A6C8"]


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


def diverging_cmap():
    """The blue-white-red map of the R heatmaps (#3B6DB3, white, #C1432B).

    Returns:
        matplotlib.colors.LinearSegmentedColormap: The colour map.
    """
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("pubready_diverging", ["#3B6DB3", "#FFFFFF", "#C1432B"])


def significance_cmap():
    """The blue-to-red map of R's scale_colour_gradient for a significance scale.

    Returns:
        matplotlib.colors.LinearSegmentedColormap: Blue for low, red for high values.
    """
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("pubready_significance", ["#3B6DB3", "#C1432B"])


def plain_tick(value, pos=None):
    """Format a log-axis tick as a plain number (0.1, 1, 10), as ggplot2 does.

    Args:
        value: The tick value.
        pos: The tick position, unused; matplotlib passes it.

    Returns:
        str: The label.
    """
    return f"{value:g}"


def log_ticks(lo, hi):
    """Ticks of the 1-2-5 sequence that cover a positive range, for a log axis.

    Args:
        lo: The smallest value on the axis.
        hi: The largest value on the axis.

    Returns:
        list: The tick positions, at least the two that bracket the range.
    """
    seq = [m * 10.0 ** k for k in range(-4, 5) for m in (1, 2, 5)]
    inside = [v for v in seq if lo <= v <= hi]
    below = [v for v in seq if v < lo]
    above = [v for v in seq if v > hi]
    if len(inside) < 2:
        inside = ([below[-1]] if below else []) + inside + ([above[0]] if above else [])
    return inside


def plain_log_axis(axis):
    """Label a log-scaled matplotlib axis with plain numbers.

    Args:
        axis: The axis, for example ``ax.xaxis``.
    """
    from matplotlib.ticker import FuncFormatter, NullFormatter
    axis.set_major_formatter(FuncFormatter(plain_tick))
    axis.set_minor_formatter(NullFormatter())


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
