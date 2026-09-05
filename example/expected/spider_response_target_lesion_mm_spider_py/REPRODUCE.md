# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.4.1`
- Image id: `4db9054289ed8c3ad23c6e14c7328d47cfba2ecaf29e993caaa92fc6f11a5c2c`
- pubplot commit: `c55152b`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.4.1 python3 script_spider_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
