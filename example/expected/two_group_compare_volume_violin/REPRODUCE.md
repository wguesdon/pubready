# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `6192160892c9cd60dc3c82d9f1fc3a7bc849f2da45750786d034ed45e9a6b073`
- Image digest: `sha256:0b27b73f590b4a5c3eb79ce87d8ed965bbf52e94918d7fd3fabfe303c547d856`
- pubplot commit: `cf6e272`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 Rscript script_violin.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
