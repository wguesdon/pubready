# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:0.4.1`
- Image id: `b893dd6f1d130c5b6dd3a3e82750620f264f759546b7cb415512e8231a8ae02b`
- pubready commit: `e8f923bef415cd1719ad862c56f1043770302fc8`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:0.4.1 python3 script_anova_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
