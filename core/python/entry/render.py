#!/usr/bin/env python3
"""Entry: figkit render (Python engine)."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBPLOT_CORE", "/opt/pubplot/core"), "python"))
from pubplot.run import run_recipe          # noqa: E402
from pubplot.spec import spec_from_config    # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--data"); p.add_argument("--sheet"); p.add_argument("--stamp")
    p.add_argument("--engine")
    p.add_argument("--out", default="pubplot_output")
    a, _ = p.parse_known_args()
    spec = spec_from_config(a.config)
    raw = a.data or spec.get("data", {}).get("file")
    if not raw:
        raise SystemExit("no input data: pass --data or set data.file in the config")
    run_recipe(spec, raw, a.out, sheet=a.sheet, stamp=a.stamp)


if __name__ == "__main__":
    main()
