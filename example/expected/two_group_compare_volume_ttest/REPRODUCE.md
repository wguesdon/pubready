# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:1.0.0`
- Image id: `adb9d0cb25c7c7836a02626b9e9da08d7912abc1ea7a07b3cc0689ac498bfa1a`
- Image digest: `sha256:4aad0f46de6f65d943f66f6dd38001a64036986b3ec44c955a52a594a496fb6e`
- pubready commit: `f66e11ade8775d5c33ce81c38926927708b194cb`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm -v "$PWD":/work -w /work \
  localhost/pubready:1.0.0 Rscript script_ttest.R
```

The script reads the bundled input copy, reruns the same test, and redraws
the figure. The stats are recomputed, so the figure and the numbers cannot
drift apart.

See `manifest_*.json` for the full package versions and the test result,
and `session_info.txt` for the complete environment.
