# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:0.4.1`
- Image id: `c6b5111edf89b3fa4919d18a52fe9f5c45ce02c5fba5672bbb6d97e988739a0e`
- pubready commit: `b35bc2f`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:0.4.1 python3 script_km_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
