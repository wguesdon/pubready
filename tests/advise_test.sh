#!/usr/bin/env bash
# Exercise the public command, including its exit codes and target engines.
set -eu
cd "$(dirname "$0")/.."
out=pubready_output/advise_tests
mkdir -p "$out"

check() {
  local code="$1" recipe="$2" label="$3" rc=0
  shift 3
  ./cli/figkit advise --format json "$@" > "$out/$label.json" || rc=$?
  if [ "$rc" -ne "$code" ]; then
    cat "$out/$label.json"
    echo "FAIL: $label returned $rc, expected $code"
    exit 1
  fi
  if [ -n "$recipe" ]; then
    grep -q '"recipe": "'"$recipe"'"' "$out/$label.json"
  fi
  echo "PASS: $label"
}

for engine in r python; do
  check 0 paired_compare "paired_$engine" --engine "$engine" \
    --data benchmarks/data/sleep.csv --x group --y extra --answer aim=compare --answer design=paired
  check 0 multi_group_compare "multi_$engine" --engine "$engine" \
    --data benchmarks/data/plant_growth.csv --x group --y weight --answer aim=compare --answer design=independent
  check 0 timecourse "repeated_$engine" --engine "$engine" \
    --data benchmarks/data/chick_weight.csv --x Time --y weight --id Chick --group Diet \
    --answer aim=compare --answer design=repeated
  check 0 survival_km "survival_$engine" --engine "$engine" \
    --data benchmarks/data/lung.csv --time time --event event --x sex \
    --answer aim=time_to_event --answer tte_aim=compare_groups
  check 0 cumulative_incidence "competing_$engine" --engine "$engine" \
    --data benchmarks/data/mgus2.csv --time etime --event event --x sex --cause 1 \
    --answer aim=time_to_event --answer competing_events=present
  for label in paired multi repeated survival competing; do
    jq -j '.command[] + "\u0000"' "$out/${label}_$engine.json" |
      xargs -0 ./tests/run_advised_command.sh "$out/plots_$engine/$label" > "$out/plot_${label}_$engine.log" 2>&1
  done
  grep -q 'Wilcoxon' "$out/plots_$engine/paired"/*/stats_*.csv
  grep -q 'one-way ANOVA' "$out/plots_$engine/multi"/*/stats_*.csv
  echo "PASS: $engine plot subset"
done
check 3 "" missing_design --data benchmarks/data/sleep.csv --x group --y extra --answer aim=compare
check 2 "" duplicate --data benchmarks/data/sleep.csv --answer aim=compare --answer aim=curve
check 4 "" conflicting_codes --data benchmarks/data/mgus2.csv --event event --time etime \
  --answer aim=time_to_event --answer competing_events=absent
check 5 "" unsupported --data benchmarks/data/sleep.csv --x group --y extra \
  --answer aim=compare --answer design=stratified
check 0 two_group_compare deterministic_a --data example/tumor_volume.csv --x group --y volume \
  --answer aim=compare --answer design=independent
check 0 two_group_compare deterministic_b --data example/tumor_volume.csv --x group --y volume \
  --answer aim=compare --answer design=independent
cmp "$out/deterministic_a.json" "$out/deterministic_b.json"

# Keep every current engine note in the operating guide for each client.
check 0 factorial_anova note_factorial --engine python \
  --data benchmarks/data/tooth_growth.csv --x dose --y len --fill supp \
  --answer aim=compare --answer design=independent
check 0 proportions note_proportions --engine python \
  --data benchmarks/data/titanic.csv --x Class --y Survived \
  --answer aim=compare --answer category_design=independent
check 0 pca note_pca --engine python \
  --data example/pca_samples.csv --group group \
  --answer aim=describe --answer describe_kind=pca
for adapter in skills/claude-code/SKILL.md skills/codex/AGENTS.md skills/opencode/AGENTS.md; do
  sed -n '/^## Engine limits/,/^## Read more/p' "$adapter" |
    tr '\n' ' ' | tr -s '[:space:]' ' ' > "$out/engine_limits.txt"
  for result in note_factorial note_proportions note_pca survival_python; do
    note=$(jq -er '.notes[0] | select(type == "string" and length > 0)' "$out/$result.json")
    grep -Fq -- "$note" "$out/engine_limits.txt"
  done
done
echo "PASS: engine notes match all three adapters"

# The tree contract already equates its leaves with recipe_options.txt.
./cli/figkit --help | sed -n '/^Recipes:/,$p' | sed 's/^Recipes://' |
  tr ', ' '\n' | sed '/^$/d' | sort > "$out/help_recipes.txt"
awk '!/^#/ && NF && $1 != "always" {print $1}' core/recipe_options.txt |
  sort > "$out/option_recipes.txt"
cmp "$out/help_recipes.txt" "$out/option_recipes.txt"
echo "PASS: help lists every tree recipe"
echo "All advise checks passed."
