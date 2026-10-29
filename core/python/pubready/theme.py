"""The pubready house style for matplotlib, matched to the R theme."""
import matplotlib
matplotlib.use("Agg")  # headless container
import matplotlib.pyplot as plt
from matplotlib.ticker import Formatter

# The first six colours are the house palette. The fourteen after them keep a
# figure with up to twenty groups distinct; palette() cycles past twenty.
PALETTE = ["#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A",
           "#56B4E9", "#E78AC3", "#8C6D31", "#66C2A5", "#F28E2B", "#A6D854",
           "#1F3A5F", "#B07AA1", "#9C755F", "#BAB0AC", "#17BECF", "#FF9DA7",
           "#4E7A27", "#D4A6C8"]


# The diameter in points of a filled ggplot2 point (shape 16, stroke 0) of size
# 1. ggplot2 sets the symbol size to size * .pt points (.pt = 72.27 / 25.4), and
# R draws shape 16 as a circle whose radius is 0.375 of the symbol size.
R_POINT_PT = 0.75 * 72.27 / 25.4


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


def log_ticks_ref(lo, hi, ref=1.0):
    """Log-axis ticks of the 1-2-5 sequence for a plot with a reference line.

    The ticks inside the data range are kept. When the data reach past ``ref``
    on one side but no tick lies on that side, the next value of the sequence on
    that side is added, so the axis reads on both sides of the line. More than 7
    ticks are thinned to the powers of ten. The R engine has the same rule in
    log_ticks_ref().

    Args:
        lo: The smallest value to show.
        hi: The largest value to show.
        ref: The reference value, for example 1 for a ratio.

    Returns:
        list: The sorted tick positions.
    """
    import numpy as np
    seq = sorted(m * 10.0 ** k for k in range(-4, 5) for m in (1, 2, 5))
    tol = 1e-9
    lo2, hi2 = min(lo, ref), max(hi, ref)
    ticks = [v for v in seq if lo2 * (1 - tol) <= v <= hi2 * (1 + tol)]
    if hi > ref and not any(t > ref * (1 + tol) for t in ticks):
        ticks.append(min(v for v in seq if v > ref * (1 + tol)))
    if lo < ref and not any(t < ref * (1 - tol) for t in ticks):
        ticks.insert(0, max(v for v in seq if v < ref * (1 - tol)))
    if len(ticks) > 7:
        ticks = [t for t in ticks if abs(np.log10(t) - round(np.log10(t))) < 1e-9]
    return sorted(set(ticks))


def wrap_label(s, width=12):
    """Wrap one label at the spaces so that no line is longer than ``width``.

    A word longer than ``width`` stays whole on its own line. The R engine has
    the same greedy rule (wrap_label), so both engines break a label at the same
    places. ``textwrap`` is not used because it treats whitespace differently.

    Args:
        s: The label.
        width: The largest count of characters on one line.

    Returns:
        str: The label with a newline at each break.
    """
    s = str(s)
    if len(s) <= width:
        return s
    words = [w for w in s.split(" ") if w]
    if not words:
        return s
    lines, cur = [], words[0]
    for w in words[1:]:
        if len(cur) + 1 + len(w) <= width:
            cur = f"{cur} {w}"
        else:
            lines.append(cur)
            cur = w
    return "\n".join(lines + [cur])


def wrap_labels(labels, width=12):
    """Wrap each label of a list with wrap_label.

    Args:
        labels: The labels.
        width: The largest count of characters on one line.

    Returns:
        list: The wrapped labels.
    """
    return [wrap_label(s, width) for s in labels]


def longest_line(s):
    """The length in characters of the longest line of a wrapped label.

    Args:
        s: The label.

    Returns:
        int: The length.
    """
    return max((len(part) for part in str(s).split("\n")), default=0)


def category_width(n, labels=None, min_width=4.2):
    """The figure width of a category axis, as the R engine computes it.

    A wrapped label whose longest line (a long word) is wider than the slot
    widens every slot to 0.085 in per character.

    Args:
        n: The count of categories.
        labels: The wrapped labels, or None.
        min_width: The smallest width in inches.

    Returns:
        float: 1.1 in per category plus 1 in, and never less than ``min_width``.
    """
    chars = max((longest_line(s) for s in labels), default=0) if labels else 0
    return max(min_width, max(1.1, 0.085 * chars) * n + 1.0)


def legend_width(labels, title=""):
    """The width of a legend right of the panel, as legend_width() of the R engine.

    Args:
        labels: The legend labels.
        title: The legend title.

    Returns:
        float: The key plus 0.075 in per character of the longest line or title.
    """
    chars = max([longest_line(s) for s in labels] + [len(str(title)), 0])
    return 0.55 + 0.075 * chars


def plain_number_labels(values):
    """Format axis breaks as plain numbers with a comma between the thousands.

    The count of decimals is the smallest one that shows every break exactly (at
    most 6), never in scientific notation, as plain_number_labels() of the R
    engine does.

    Args:
        values: The break values.

    Returns:
        list: The labels.
    """
    import numpy as np
    vals = np.asarray(values, dtype=float)
    fin = vals[np.isfinite(vals)]
    digits = 0
    while digits < 6 and len(fin) and np.any(
            np.abs(np.round(fin, digits) - fin) > 1e-9 * np.maximum(1, np.abs(fin))):
        digits += 1
    return [f"{v + 0.0:,.{digits}f}" if np.isfinite(v) else "" for v in vals]


class PlainNumberFormatter(Formatter):
    """A matplotlib tick formatter that applies plain_number_labels to all ticks."""

    def __init__(self):
        self.locs = []

    def set_locs(self, locs):
        """Store the tick positions, so every label uses the same decimals.

        Args:
            locs: The tick positions.
        """
        self.locs = list(locs)

    def __call__(self, x, pos=None):
        """Format one tick.

        Args:
            x: The tick value.
            pos: The tick index, unused.

        Returns:
            str: The label.
        """
        return plain_number_labels(list(self.locs) + [x])[-1]


def plain_number_axis(axis):
    """Label a linear matplotlib axis with plain numbers and thousands separators.

    Args:
        axis: The axis, for example ``ax.yaxis``.
    """
    axis.set_major_formatter(PlainNumberFormatter())


def geom_top(values, geom="box"):
    """The top of the drawn geom of one group, as geom_top() of the R engine.

    The seaborn violin is cut at the data (cut=0), so it ends at the largest
    value; a bar ends at the mean plus its standard error when no point is
    higher.

    Args:
        values: The values of the group.
        geom: The geometry: box, violin, bar or raincloud.

    Returns:
        float: The top, or NaN for an empty group.
    """
    import numpy as np
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return float("nan")
    top = float(v.max())
    if geom == "bar" and len(v) >= 2:
        top = max(top, float(v.mean() + v.std(ddof=1) / np.sqrt(len(v))))
    return top


def rain_offsets(n, lo=-0.18, hi=-0.06):
    """Deterministic horizontal offsets for the points of a raincloud.

    The n values in [lo, hi] come from the golden-ratio sequence, so the points
    spread to the left of the box without a random seed, as in the R engine.

    Args:
        n: The count of points.
        lo: The leftmost offset.
        hi: The rightmost offset.

    Returns:
        numpy.ndarray: The offsets.
    """
    import numpy as np
    k = np.arange(1, n + 1)
    return lo + (hi - lo) * np.mod(k * 0.6180339887498949, 1.0)


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
            jit = rain_offsets(len(vals))
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
