"""Reading tidy data and checksumming inputs (mirrors the R io helpers)."""
import hashlib
import os

import pandas as pd


def read_tidy(path, sheet=None):
    if not os.path.exists(path):
        raise FileNotFoundError(f"input file not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(path, sheet_name=0 if sheet is None else sheet)
    return pd.read_csv(path)


def file_checksums(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            md5.update(chunk)
            sha.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha.hexdigest()}
