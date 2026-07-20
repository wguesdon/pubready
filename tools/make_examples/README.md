# make_examples

Maintainer tooling that rebuilds the example assets. None of it is part of a
figure bundle, and an end user never runs it. It lives here, out of `example/`,
so that folder shows only inputs and reference bundles.

Two stages.

**1. Generators (`make_*.R`).** Deterministic, seeded scripts. Each writes one
input CSV into `example/`, the file a scientist would drop in. Run one inside the
pinned container from the repo root:

```bash
podman run --rm -v "$PWD":/work -w /work localhost/pubplot:0.4.1 \
  Rscript tools/make_examples/make_survival_data.R
```

**2. Reference bundles (`generate_expected.sh`).** Runs `figkit plot` on those
inputs for both engines and writes the committed bundles under `example/expected/`.
It uses fixed stamps and a fixed created date, so reruns diff cleanly:

```bash
./tools/make_examples/generate_expected.sh
```

## The volcano fixture

`make_volcano_data.R` writes the small synthetic DE table (1200 genes) that the
committed volcano example uses, keeping the repo light. `make_volcano_airway.R` is
optional: it rebuilds the real airway RNA-seq benchmark with DESeq2 into
`example/de_results_airway.csv`, which is gitignored.

## Determinism check

`check_determinism.sh` regenerates every bundle into a temp tree and confirms each
`stats_*.csv` matches the committed copy. Stats carry the scientific numbers with
no environment stamps, so a clean match proves same input plus same recipe gives
the same result. It runs one container per figure, so it is slow. Run it on
demand, not in fast CI:

```bash
./tools/make_examples/check_determinism.sh
```
