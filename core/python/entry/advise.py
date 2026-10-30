#!/usr/bin/env python3
"""Entry point for the Python-only advice walker."""

import argparse
import os
import sys
import yaml

sys.path.insert(0, os.path.join(os.environ.get("PUBREADY_CORE", "/opt/pubready/core"), "python"))
from pubready.advise import AdviceError, advise, empty_result, load_tree, render_result  # noqa: E402
from pubready.spec import PLOT_TEXT_OPTIONS, PLOT_VALUE_OPTIONS  # noqa: E402


class AdviceParser(argparse.ArgumentParser):
    """An argument parser that returns public usage errors."""

    def error(self, message):
        """Raise a usage error instead of printing an unstructured response.

        Args:
            message: Parser error text.
        """
        raise AdviceError(message)


class UniqueOption(argparse.Action):
    """Accept identical repetitions and reject contradictory plot flags."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Store an explicit value once.

        Args:
            parser: The argument parser.
            namespace: Parsed arguments.
            values: The supplied value.
            option_string: The option spelling.
        """
        if hasattr(namespace, self.dest) and getattr(namespace, self.dest) != values:
            raise AdviceError(f"Option {option_string} has two different values.")
        setattr(namespace, self.dest, values)


def main():
    """Parse explicit flags, print one result and return its public exit code."""
    parser = AdviceParser(allow_abbrev=False, argument_default=argparse.SUPPRESS)
    for key in ("data", "engine", "out") + PLOT_TEXT_OPTIONS:
        parser.add_argument(f"--{key}", action=UniqueOption)
    for key, kind in PLOT_VALUE_OPTIONS:
        parser.add_argument(f"--{key}", type=kind, action=UniqueOption)
    parser.add_argument("--paired", action="store_true", default=argparse.SUPPRESS)
    parser.add_argument("--answer", action="append", default=[])
    parser.add_argument("--answers")
    parser.add_argument("--format", choices=("text", "json"), default="text")
    output_format = "json" if "--format=json" in sys.argv or any(
        a == "--format" and b == "json" for a, b in zip(sys.argv, sys.argv[1:])) else "text"
    try:
        args = vars(parser.parse_args())
        entries = args.pop("answer")
        filename = args.pop("answers", None)
        output_format = args.pop("format")
        result = advise(args, entries, filename)
    except AdviceError as error:
        result = empty_result(None)
        result.update(status="usage_error", exit_code=2, conflicts=[{"reason": str(error)}])
    except (ValueError, OSError, yaml.YAMLError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    print(render_result(result, output_format), end="")
    return result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
