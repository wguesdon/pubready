# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:0.4.1`
- Image id: `b893dd6f1d130c5b6dd3a3e82750620f264f759546b7cb415512e8231a8ae02b`
- Image digest: `sha256:097c136b375dd40cd790ef1528d70d748ab67d06fc44e73c70de82b4a4782d17`
- pubready commit: `e8f923bef415cd1719ad862c56f1043770302fc8`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:0.4.1 Rscript script_prop.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
