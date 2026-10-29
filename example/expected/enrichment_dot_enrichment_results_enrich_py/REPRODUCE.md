# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:2.0.0`
- Image id: `1e752da378ea70f07f72cd70305231f4f079f1b4e2d99c87cdf2e3d35e1f44b2`
- pubready commit: `1b8bb2fd31d9315d60137be3bc664c2a3e33f671`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:2.0.0 python3 script_enrich_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
