# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:3.0.0`
- Image id: `0511fa40d3aa87e86b6e5e4b4852cf88897b4ef8dabb1d2646f211c7c72b8cf6`
- pubready commit: `4657df8148b43cbd66ec3a2bc04a756a38e7843c`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:3.0.0 python3 script_embedding_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
