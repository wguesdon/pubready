# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:0.4.1`
- Image id: `be5dbb7703261c36a854ae80ed11cfc2a61dc669053a2c0cd63fffc3aa0a1bf8`
- pubready commit: `f82215f4b833c9d70432f64150af9bc357a8d5e4`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:0.4.1 python3 script_paired_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
