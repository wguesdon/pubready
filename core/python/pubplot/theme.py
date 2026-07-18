"""The pubplot house style for matplotlib, matched to the R theme."""
import matplotlib
matplotlib.use("Agg")  # headless container
import matplotlib.pyplot as plt

PALETTE = ["#3B6DB3", "#C1432B", "#2E8B57", "#7A5195", "#E0A100", "#5A5A5A"]


def palette(n, custom=None):
    cols = custom if custom else PALETTE
    return [cols[i % len(cols)] for i in range(n)]


def apply_pub_style(ax):
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


def new_fig(w=3.8, h=4.0):
    fig, ax = plt.subplots(figsize=(w, h))
    return fig, ax
