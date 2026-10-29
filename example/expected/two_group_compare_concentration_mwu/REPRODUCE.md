# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:2.0.0`
- Image id: `1e752da378ea70f07f72cd70305231f4f079f1b4e2d99c87cdf2e3d35e1f44b2`
- Image digest: `sha256:f22d119d4121702c30e3632421cdf83f9bc1b190c13c3547876874aed2f0419e`
- pubready commit: `1b8bb2fd31d9315d60137be3bc664c2a3e33f671`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:2.0.0 Rscript script_mwu.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
