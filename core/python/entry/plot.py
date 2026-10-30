#!/usr/bin/env python3
"""Entry: figkit plot (Python engine)."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBREADY_CORE", "/opt/pubready/core"), "python"))
from pubready.run import run_entry, run_recipe            # noqa: E402
from pubready.spec import add_plot_arguments, spec_from_args   # noqa: E402


def main():
    p = add_plot_arguments(argparse.ArgumentParser())
    a, _ = p.parse_known_args()
    spec = spec_from_args(a)
    run_recipe(spec, a.data, a.out, sheet=a.sheet, stamp=a.stamp)


if __name__ == "__main__":
    run_entry(main)
