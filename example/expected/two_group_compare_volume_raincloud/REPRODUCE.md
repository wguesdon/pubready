# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.4.0`
- Image id: `4cf3e55927bb05afb9692649c012295c7599388c7be89ffaf893c590e80bd530`
- Image digest: `sha256:feccf0bcc3ec0b49b272c9f5910ab221b2dbb649016ee2ee654df7e503e3bead`
- pubplot commit: `f26dea0`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.4.0 Rscript script_raincloud.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
