#!/usr/bin/env bash
# End-to-end test of the uv install. It installs the committed HEAD from a
# git+file URL into a temporary uv tool folder, as `uv tool install` from GitHub
# does, and checks that the installed figkit finds its files, records the
# installed commit and writes a complete bundle. It needs uv, Podman and the
# built image, and uv fetches hatchling from the network. Commit the changes
# first: the test installs HEAD, not the working tree.
cd "$(dirname "$0")/.."
if ! command -v uv > /dev/null; then
  echo "SKIP: uv is not installed"; exit 0
fi

repo="$PWD"
commit="$(git rev-parse HEAD)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
export UV_TOOL_DIR="$tmp/tools" UV_TOOL_BIN_DIR="$tmp/bin"
fail=0

echo "[1] uv tool install from git+file at $commit"
if ! uv tool install --quiet "pubready @ git+file://$repo@$commit"; then
  echo "FAIL: uv tool install"; exit 1
fi
fk="$tmp/bin/figkit"
if [ -x "$fk" ]; then echo "PASS: figkit is on the tool PATH"; else echo "FAIL: no figkit"; exit 1; fi

echo "[2] figkit version names the installed commit"
if "$fk" version | grep -q "commit: ${commit:0:7}"; then
  echo "PASS: the version names ${commit:0:7}"
else
  echo "FAIL: expected commit ${commit:0:7}"; "$fk" version; fail=1
fi

echo "[3] figkit where holds the adapters, the reference notes and the example"
where="$("$fk" where)"
for f in skills/claude-code/SKILL.md skills/codex/AGENTS.md skills/opencode/AGENTS.md \
         reference/decision_tree.md reference/decision_tree.yaml docs/how_it_works.md container/Containerfile \
         core/r/bootstrap.R core/python/pubready/run.py core/recipe_options.txt \
         example/tumor_volume.csv; do
  if [ ! -f "$where/$f" ]; then echo "  MISSING: $f"; fail=1; fi
done
[ "$fail" -eq 0 ] && echo "PASS: the installed tree is complete -> $where"

echo "[4] the installed figkit writes a bundle that records the installed commit"
work="$tmp/work"
mkdir -p "$work"
cp "$where/example/tumor_volume.csv" "$work/"
for engine in r python; do
  if (cd "$work" && "$fk" plot --engine "$engine" --recipe two_group_compare \
        --data tumor_volume.csv --x group --y volume --out "out_$engine" > /dev/null 2>&1); then
    if grep -q "$commit" "$work"/out_"$engine"/*/manifest_*.json && \
       grep -q '"git_dirty": false' "$work"/out_"$engine"/*/report_provenance.json; then
      echo "PASS: $engine bundle records $commit, not dirty"
    else
      echo "FAIL: $engine bundle does not record the installed commit"; fail=1
    fi
  else
    echo "FAIL: $engine plot from the installed figkit"; fail=1
  fi
done

if (cd "$work" && "$fk" advise --data tumor_volume.csv --x group --y volume \
      --answer aim=compare --answer design=independent --format json > advise.json) && \
      grep -q '"recipe": "two_group_compare"' "$work/advise.json"; then
  echo "PASS: installed advise reads the tree and selects the recipe"
else
  echo "FAIL: installed advise"; fail=1
fi

echo
if [ "$fail" -eq 0 ]; then echo "INSTALL TEST PASSED"; else echo "INSTALL TEST FAILED"; exit 1; fi
