# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubplot:0.1.0`
- Image id: `8e759d0550177f4eec47a0d1ce44e8184a08b4703ade063539cdf2f9b2ed2a59`
- Image digest: `sha256:abd4a9e761aaeb5f2694699bd76d1800abf656cdedd37d8183c1f8e3ebf07766`
- pubplot commit: `5040627`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubplot:0.1.0 Rscript script_anova.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
