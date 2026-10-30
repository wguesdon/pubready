# Reproduce this figure

This bundle is self-contained. To regenerate the figure in the exact
environment it was made in:

## 1. Get the pinned environment

- Container image: `localhost/pubready:4.0.0`
- Image id: `e76db7ec4edffd2951b7bef6a38a1e4abf91499f706c57fa73d8b8efcdac01b5`
- Image digest: `sha256:0f57027a3fbca5ec9d41d74ca91226fbdf45d94b1a85f271637d12733db5d483`
- pubready commit: `567ce8900bdb5db88369d131c8d787b042b3af57`

## 2. Run the standalone script

From inside this folder:

```bash
podman run --rm --security-opt label=disable -v "$PWD":/work -w /work \
  localhost/pubready:4.0.0 python3 script_twoway_py.py
```

The script reads the bundled input copy, reruns the analysis, and redraws
the figure. See `manifest_*.json` for package versions and
`session_info.txt` for the full environment.
