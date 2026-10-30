#!/usr/bin/env python3
"""Entry point for the Python-only advice walker."""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBREADY_CORE", "/opt/pubready/core"), "python"))
from pubready.advise import advise, render_result  # noqa: E402
from pubready.spec import PLOT_TEXT_OPTIONS, PLOT_VALUE_OPTIONS  # noqa: E402


def main():
    """Parse explicit flags, print one result and return its public exit code."""
    parser = argparse.ArgumentParser(allow_abbrev=False, argument_default=argparse.SUPPRESS)
    for key in ("data", "engine", "out") + PLOT_TEXT_OPTIONS:
        parser.add_argument(f"--{key}")
    for key, kind in PLOT_VALUE_OPTIONS:
        parser.add_argument(f"--{key}", type=kind)
    parser.add_argument("--paired", action="store_true", default=argparse.SUPPRESS)
    parser.add_argument("--answer", action="append", default=[])
    parser.add_argument("--answers")
    parser.add_argument("--format", choices=("text", "json"), default="text")
    args = vars(parser.parse_args())
    entries = args.pop("answer")
    filename = args.pop("answers", None)
    output_format = args.pop("format")
    try:
        result = advise(args, entries, filename)
    except (ValueError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(render_result(result, output_format), end="")
    return result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
