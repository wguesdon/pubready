#!/usr/bin/env bash
# Execute one returned argument array without shell evaluation.
set -eu
out="$1"
shift
"$@" --out "$out"
