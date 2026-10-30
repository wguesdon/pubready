#!/usr/bin/env python3
"""Entry: figkit plot (Python engine)."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBREADY_CORE", "/opt/pubready/core"), "python"))
from pubready.run import run_entry, run_recipe            # noqa: E402
from pubready.spec import (add_plot_arguments, check_recipe_options,  # noqa: E402
                           given_options, spec_from_args)


def main():
    p = add_plot_arguments(argparse.ArgumentParser(allow_abbrev=False))
    a = p.parse_args()
    check_recipe_options(a.recipe, given_options(sys.argv[1:]))
    spec = spec_from_args(a)
    run_recipe(spec, a.data, a.out, sheet=a.sheet, stamp=a.stamp)


if __name__ == "__main__":
    run_entry(main)
