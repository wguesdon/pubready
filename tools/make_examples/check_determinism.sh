#!/usr/bin/env bash
# Determinism check: regenerate every reference bundle into a temp tree and
# confirm each stats_*.csv matches the committed copy under example/expected.
# Stats hold the scientific numbers with no environment stamps (no git commit,
# image digest, or timestamp), so a clean match proves same input + same recipe
# gives the same result. Slow: one container run per figure. Run on demand.
set -euo pipefail
cd "$(dirname "$0")/../.."

# The output dir must sit under the repo root. figkit mounts the working dir into
# the container, so an --out path outside that mount is written inside the
# container and never lands on the host. Use a repo-relative temp dir, not /tmp.
tmp="$(mktemp -d ".determinism.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

echo "Regenerating reference bundles into $tmp (this is slow)..."
PUBREADY_EXPECTED_OUT="$tmp" ./tools/make_examples/generate_expected.sh >/dev/null

fail=0
for want in example/expected/*/stats_*.csv; do
  bundle="$(basename "$(dirname "$want")")"
  got="$tmp/$bundle/$(basename "$want")"
  if [ ! -f "$got" ]; then
    echo "MISSING regenerated stats for $bundle"
    fail=1
    continue
  fi
  if ! diff -q "$want" "$got" >/dev/null; then
    echo "DIFF in $bundle/$(basename "$want")"
    fail=1
  fi
done

if [ "$fail" -eq 0 ]; then
  echo "determinism OK: all stats tables reproduced"
else
  echo "determinism FAILED: see diffs above" >&2
  exit 1
fi
