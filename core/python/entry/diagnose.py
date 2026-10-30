#!/usr/bin/env python3
"""Entry: figkit diagnose (Python engine). Show the assumption evidence behind
the test choice so a scientist can confirm or override it. Writes a normality QQ
panel per checked quantity plus a machine-readable report, and prints the
recommended test. It does not draw the final figure."""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get("PUBREADY_CORE", "/opt/pubready/core"), "python"))
from pubready.bundle import _ts, slugify, write_json   # noqa: E402
from pubready.io import read_tidy          # noqa: E402
from pubready.qc import qc_normality       # noqa: E402
from pubready.recipes import get_recipe    # noqa: E402
from pubready.run import run_entry         # noqa: E402
from pubready.spec import (add_plot_arguments, check_recipe_options,  # noqa: E402
                           given_options, spec_from_args)


def _fmt(x, g=False):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return f"{x:.3g}" if g else f"{x:.3f}"


def main():
    # diagnose accepts every option of figkit plot, so the recipe sees the same spec.
    p = add_plot_arguments(argparse.ArgumentParser(allow_abbrev=False))
    a = p.parse_args()
    check_recipe_options(a.recipe, given_options(sys.argv[1:]))

    spec = spec_from_args(a)
    df = read_tidy(a.data, a.sheet)
    result = get_recipe(spec["recipe"])(df, spec)

    rec = (result.get("test_meta") or {}).get("name") or "(chosen at plot time)"
    ts = _ts(a.stamp)
    base = slugify(result.get("label") or spec["data"].get("y") or
                   spec["data"].get("x") or spec["recipe"])
    ddir = os.path.join(a.out, f"diagnose_{spec['recipe']}_{base}_{ts}")
    os.makedirs(ddir, exist_ok=True)

    print(f"Recipe:            {spec['recipe']}")
    print(f"Recommended test:  {rec}")

    qc = result.get("qc")
    if qc:
        fig, panels = qc_normality(qc)
        png = os.path.join(ddir, "qc_normality.png")
        fig.savefig(png, dpi=200, facecolor="white")
        fig.savefig(os.path.join(ddir, "qc_normality.pdf"), facecolor="white")

        print(f"Normality check on: {qc['quantity']}\n")
        for s in panels:
            print(f"  {s['panel']:<16} n={s['n']:<4} Shapiro p={_fmt(s['shapiro_p'], True):<8} "
                  f"skew={_fmt(s['skewness']):<7} kurtosis={_fmt(s['kurtosis'])}")
        print("\nShapiro-Wilk is a default, not a verdict: it is underpowered at small n")
        print("and over-rejects at large n. Read the QQ plot; mild deviation is usually fine.")
        print(f"\nQQ plot: {png}")
        report = {"recipe": spec["recipe"], "engine": "python", "recommended_test": rec,
                  "normality_quantity": qc["quantity"], "panels": panels,
                  "qq_plot": "qc_normality.png"}
    else:
        print("\nThis recipe does not choose its test from a normality check, so there is")
        print("no QQ diagnostic. Run figkit advise to check the recipe and its test choices.")
        report = {"recipe": spec["recipe"], "engine": "python", "recommended_test": rec,
                  "normality_quantity": None}

    write_json(report, os.path.join(ddir, "diagnostics.json"))

    print(f"\nWrote diagnostics: {ddir}")
    print("Confirm or override the test, then run `figkit plot` to draw the figure.")


if __name__ == "__main__":
    run_entry(main)
