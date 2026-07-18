# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `6192160892c9cd60dc3c82d9f1fc3a7bc849f2da45750786d034ed45e9a6b073`
- pubplot commit: `cf6e272`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 python3 script_demo_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
