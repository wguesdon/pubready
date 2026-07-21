# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.4.1`
- Image id: `de93d6bacf2e827e806f040319b0b20ea7dc1f65c20388406bf4343fb26e2667`
- Image digest: `sha256:e47e117619ee83ef8194b4a007960cce130937a5827590cff65f0330f817ba25`
- pubplot commit: `15c0014`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.4.1 Rscript script_enrich.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
