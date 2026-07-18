#!/usr/bin/env python3
"""Entry: figkit inspect (Python engine)."""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBPLOT_CORE", "/opt/pubplot/core"), "python"))
from pubplot.io import read_tidy                      # noqa: E402
from pubplot.run import inspect_data, print_inspect   # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", required=True)
    p.add_argument("--sheet"); p.add_argument("--format", default="text")
    p.add_argument("--engine")
    a, _ = p.parse_known_args()
    df = read_tidy(a.data, a.sheet)
    info = inspect_data(df)
    if a.format == "json":
        print(json.dumps(info, indent=2))
    else:
        print_inspect(info)


if __name__ == "__main__":
    main()
