"""Minimal, logged cleaning for x/y recipes (mirrors the R clean_tidy)."""
import numpy as np
import pandas as pd


def clean_xy(df, x, y):
    df = df.copy()
    steps = []
    n0 = len(df)
    if x not in df.columns:
        raise ValueError(f"column '{x}' not found in data")
    if y not in df.columns:
        raise ValueError(f"column '{y}' not found in data")

    yv = pd.to_numeric(df[y], errors="coerce")
    coerced = int((yv.isna() & df[y].notna()).sum())
    if coerced:
        steps.append(f"Coerced column '{y}' to numeric; {coerced} value(s) became NA.")
    df[y] = yv

    mask = df[x].notna() & df[y].notna() & np.isfinite(df[y])
    dropped = int((~mask).sum())
    if dropped:
        steps.append(f"Dropped {dropped} row(s) with missing '{x}' or '{y}' "
                     f"({n0} -> {int(mask.sum())} rows).")
    df = df[mask].reset_index(drop=True)

    if not steps:
        steps = ["No cleaning applied; input used as-is."]
    return df, steps
