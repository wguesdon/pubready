# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.3.0`
- Image id: `2d43d34f07864b4572a53398e14ad3d1631527725b7665e04f5b88bdbf622240`
- Image digest: `sha256:f0d2e1c018a637f042a69afaca8805e828f081f54ce88287af9643aa227417a5`
- pubplot commit: `1007a90`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.3.0 Rscript script_demo.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
