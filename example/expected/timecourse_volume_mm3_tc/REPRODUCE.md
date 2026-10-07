# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:0.4.1`
- Image id: `be5dbb7703261c36a854ae80ed11cfc2a61dc669053a2c0cd63fffc3aa0a1bf8`
- Image digest: `sha256:339bf57d4f33cf35d6ffa492d0e2faccd3ab5e0a435e20512b201d2c5fffd9e4`
- pubready commit: `f82215f4b833c9d70432f64150af9bc357a8d5e4`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:0.4.1 Rscript script_tc.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
