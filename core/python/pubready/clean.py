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


def coerce_numeric_columns(df, cols=None, threshold=0.8):
    """Classify and coerce the numeric columns of a table.

    A column is numeric when it is already numeric, or when 80% or more of its
    cells that are not missing read as finite numbers. A cell of such a column
    that does not read as a number, for example "<0.1" or "45 y", becomes NaN,
    and a step records the count and the first such value. Every other column
    is text. The R engine has the same rule.

    Args:
        df: The table.
        cols: The columns to check, or None for every column.
        threshold: The share of readable cells that makes a column numeric.

    Returns:
        tuple: The table, the list of numeric columns, the list of text columns
        and the list of steps.
    """
    df = df.copy()
    steps, numeric, text = [], [], []
    for cn in (list(df.columns) if cols is None else list(cols)):
        v = df[cn]
        if pd.api.types.is_numeric_dtype(v) and not pd.api.types.is_bool_dtype(v):
            numeric.append(cn)
            continue
        present = v.notna() & (v.astype(str).str.strip() != "")
        num = pd.to_numeric(v.where(present), errors="coerce")
        ok = present & num.notna() & np.isfinite(num.astype(float))
        n_present = int(present.sum())
        if n_present > 0 and int(ok.sum()) / n_present >= threshold:
            bad = present & ~ok
            if bad.any():
                steps.append(f"Coerced column '{cn}' to numeric; {int(bad.sum())} value(s) "
                             f"became NA, for example '{v[bad].iloc[0]}'.")
            df[cn] = num.where(ok).astype(float)
            numeric.append(cn)
        else:
            text.append(cn)
    return df, numeric, text, steps
