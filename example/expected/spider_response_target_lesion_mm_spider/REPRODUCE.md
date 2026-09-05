# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.4.1`
- Image id: `4db9054289ed8c3ad23c6e14c7328d47cfba2ecaf29e993caaa92fc6f11a5c2c`
- Image digest: `sha256:4466918382e0c859745bafb62bbd098cc39c7a5f4ff607541e9cd90388e03db6`
- pubplot commit: `c55152b`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.4.1 Rscript script_spider.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
