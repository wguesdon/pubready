# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `ee9e40c9eb1d90c899f403012aabd98793acbeb39521f197dc3e532e26749cfb`
- Image digest: `sha256:4e3fbeae25d9038f7ddc12c0bb5713b87f5669a976b573faf2abfd4c8c30fcbc`
- pubplot commit: `90f0285`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 Rscript script_km.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
