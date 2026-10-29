# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:3.0.0`
- Image id: `0511fa40d3aa87e86b6e5e4b4852cf88897b4ef8dabb1d2646f211c7c72b8cf6`
- Image digest: `sha256:6af675049bd5433b2d94050941a2f5aaa8644e1656f74e21f139a93e36308290`
- pubready commit: `4657df8148b43cbd66ec3a2bc04a756a38e7843c`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:3.0.0 Rscript script_bar.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
