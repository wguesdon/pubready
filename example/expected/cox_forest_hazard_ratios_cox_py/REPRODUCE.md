# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.3.0`
- Image id: `2d43d34f07864b4572a53398e14ad3d1631527725b7665e04f5b88bdbf622240`
- pubplot commit: `1007a90`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.3.0 python3 script_cox_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
