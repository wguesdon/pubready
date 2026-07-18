# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `7113ab1438c72d6287b277e29ae50e1612e26d679e505d5cdac6f799ff52ff79`
- Image digest: `sha256:193a667e0981eb819e289eff0b06cbb512622d94fdb508a6e26fd81b8f8a7c74`
- pubplot commit: `3b5b248`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 Rscript script_demo.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
