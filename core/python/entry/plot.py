#!/usr/bin/env python3
"""Entry: figkit plot (Python engine)."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBPLOT_CORE", "/opt/pubplot/core"), "python"))
from pubplot.run import run_recipe          # noqa: E402
from pubplot.spec import spec_from_args     # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--recipe", required=True)
    p.add_argument("--data", required=True)
    p.add_argument("--x"); p.add_argument("--y")
    p.add_argument("--fill"); p.add_argument("--facet")
    p.add_argument("--time"); p.add_argument("--event"); p.add_argument("--covariates")
    p.add_argument("--annotation"); p.add_argument("--scale"); p.add_argument("--cluster")
    p.add_argument("--id"); p.add_argument("--group"); p.add_argument("--label")
    p.add_argument("--fc_cutoff", type=float); p.add_argument("--p_cutoff", type=float)
    p.add_argument("--top_n", type=int)
    p.add_argument("--test", default="auto")
    p.add_argument("--paired", action="store_true")
    p.add_argument("--p_adjust", default="none")
    p.add_argument("--geom"); p.add_argument("--theme")
    p.add_argument("--xlab"); p.add_argument("--ylab"); p.add_argument("--title")
    p.add_argument("--palette"); p.add_argument("--sheet"); p.add_argument("--stamp")
    p.add_argument("--engine")  # consumed by the figkit dispatcher; ignored here
    p.add_argument("--out", default="pubplot_output")
    a, _ = p.parse_known_args()
    spec = spec_from_args(a)
    run_recipe(spec, a.data, a.out, sheet=a.sheet, stamp=a.stamp)


if __name__ == "__main__":
    main()
