"""Reading tidy data and checksumming inputs (mirrors the R io helpers).

One read rule for both engines. Every cell of a CSV file is read as text, the
spaces and tabs around each value are removed, and only an empty cell or "NA"
reads as missing. A cell that is a number with a comma ("85,45" or "4,380")
stops the run, because the comma can be a decimal mark or a thousands mark. A
column becomes numeric when each present cell is a number with "." as the
decimal mark, and logical when each present cell is TRUE or FALSE. Every other
column stays text. The R engine (core/r/R/io.R) applies the same rule.
"""
import hashlib
import os

import numpy as np
import pandas as pd

#: The cells that read as missing.
READ_NA_TOKENS = ("", "NA")
#: A number with "." as the decimal mark and no thousands mark.
READ_NUMBER_RE = r"^[+-]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][+-]?[0-9]+)?$"
#: A number with a comma: a decimal comma or a thousands mark.
READ_COMMA_NUMBER_RE = r"^[+-]?[0-9]*(,[0-9]+)+([.][0-9]*)?([eE][+-]?[0-9]+)?$"
#: A whole number, which pandas stores as an integer when no cell is missing.
READ_INTEGER_RE = r"^[+-]?[0-9]+$"
#: The cells that read as logical.
READ_TRUE_TOKENS = ("TRUE", "True", "true")
READ_FALSE_TOKENS = ("FALSE", "False", "false")
#: The characters that the read rule removes around a value.
READ_SPACE = " \t"


def read_counts_text(counts):
    """Join "'name' (count)" pairs for a read step.

    Args:
        counts: A dict of column name to count, in column order.

    Returns:
        str: The joined pairs, as the R engine writes them.
    """
    return ", ".join(f"'{k}' ({v})" for k, v in counts.items())


def _as_text(v):
    """Give one cell of an Excel text column as text, or None when missing."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return v if isinstance(v, str) else str(v)


def _numeric_column(tv, present):
    """Convert a text column whose present cells are numbers.

    Args:
        tv: The trimmed text column, with None for a missing cell.
        present: The mask of the present cells.

    Returns:
        pandas.Series: int64 when each cell is a whole number and none is
        missing, as pandas reads such a column, else float64.
    """
    pv = tv[present]
    if bool(present.all()) and len(pv) and bool(pv.str.fullmatch(READ_INTEGER_RE).all()):
        try:
            return pd.Series([int(v) for v in tv], index=tv.index, dtype="int64")
        except OverflowError:
            pass
    return pd.Series([float(v) if p else np.nan for v, p in zip(tv, present)],
                     index=tv.index, dtype="float64")


def tidy_text_columns(df):
    """Apply the read rule to the text columns of a table.

    Args:
        df: The table as read, with text cells as str.

    Returns:
        tuple: The table and the list of read steps for the data log.

    Raises:
        ValueError: If a cell is a number with a comma.
    """
    steps = []
    old_names = [str(c) for c in df.columns]
    new_names = [c.strip(READ_SPACE) for c in old_names]
    n_renamed = sum(a != b for a, b in zip(old_names, new_names))
    if n_renamed:
        steps.append(f"Removed the spaces around {n_renamed} column name(s).")
    trimmed, missing = {}, {}
    columns = []
    for i, cn in enumerate(new_names):
        v = df.iloc[:, i]
        if v.dtype != object:
            columns.append(v)
            continue
        tv = pd.Series([_as_text(x) for x in v], index=v.index, dtype=object)
        present = tv.notna()
        stripped = tv.where(~present, tv[present].str.strip(READ_SPACE))
        n_trim = int((present & (stripped != tv)).sum())
        if n_trim:
            trimmed[cn] = n_trim
        is_na_token = present & stripped.isin(READ_NA_TOKENS)
        if is_na_token.any():
            missing[cn] = int(is_na_token.sum())
        tv = stripped.where(~is_na_token, None)
        present = tv.notna()
        comma = present & tv.where(present, "").str.fullmatch(READ_COMMA_NUMBER_RE)
        if comma.any():
            raise ValueError(
                f"column '{cn}' holds '{tv[comma].iloc[0]}', a number with a comma. "
                "pubready reads '.' as the decimal mark and no thousands mark. Replace "
                "each decimal comma with '.' and remove each thousands mark, then run again.")
        pv = tv[present]
        if bool(pv.str.fullmatch(READ_NUMBER_RE).all()):
            columns.append(_numeric_column(tv, present))
        elif bool(pv.isin(READ_TRUE_TOKENS + READ_FALSE_TOKENS).all()):
            flags = tv.isin(READ_TRUE_TOKENS)
            columns.append(flags.astype(bool) if bool(present.all())
                           else flags.astype(object).where(present, np.nan))
        else:
            columns.append(tv.where(present, np.nan))
    if columns:
        df = pd.concat(columns, axis=1, ignore_index=True)
    df.columns = new_names
    if trimmed:
        steps.append(f"Removed the spaces around {sum(trimmed.values())} value(s): "
                     f"{read_counts_text(trimmed)}.")
    if missing:
        steps.append(f"Read {sum(missing.values())} cell(s) that hold \"NA\" or nothing as "
                     f"missing: {read_counts_text(missing)}.")
    return df, steps


def read_tidy_logged(path, sheet=None):
    """Read a tidy CSV or Excel file with the read rule.

    Args:
        path: The input file.
        sheet: The Excel sheet, or None for the first sheet.

    Returns:
        tuple: The table and the list of read steps for the data log.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a cell is a number with a comma.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"input file not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls"):
        # An Excel cell keeps the type that the sheet gives it. The text cells
        # get the read rule.
        raw = pd.read_excel(path, sheet_name=0 if sheet is None else sheet,
                            keep_default_na=False, na_values=[""])
    else:
        raw = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    return tidy_text_columns(raw)


def read_tidy(path, sheet=None):
    """Read a tidy CSV or Excel file into a data frame, with the read rule.

    Args:
        path: The input file.
        sheet: The Excel sheet, or None for the first sheet.

    Returns:
        pandas.DataFrame: The table.
    """
    return read_tidy_logged(path, sheet)[0]


def file_checksums(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            md5.update(chunk)
            sha.update(chunk)
    return {"md5": md5.hexdigest(), "sha256": sha.hexdigest()}
