# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:1.0.0`
- Image id: `b3c6d3cc597f0d6a85f7685cd0840776ca2773f64c4acdb958cae700bc785f83`
- pubready commit: `31d6bbe65e2043a4e8204e63902cf8f8a135acf6`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:1.0.0 python3 script_roc_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
