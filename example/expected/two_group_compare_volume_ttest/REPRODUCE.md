# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `0fdf1a82f66f110d35bc2f0fe8b60a9bacc9d3c2d190d6c71e44034c0e4417e1`
- Image digest: `sha256:994931c28c7be349e19f358de1261a4886f092519ef389779e2d0852bde97edf`
- pubplot commit: `d8ea256`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 Rscript script_ttest.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
