import pandas as pd
import pytest

from pubready.clean import clean_xy


def test_clean_xy_passes_clean_data_through():
    df = pd.DataFrame({"group": ["A", "A", "B", "B"], "value": [1, 2, 3, 4]})
    out, steps = clean_xy(df, "group", "value")
    assert len(out) == 4
    assert steps == ["No cleaning applied; input used as-is."]


def test_clean_xy_drops_missing_and_logs():
    df = pd.DataFrame({"group": ["A", "A", None, "B"], "value": [1, None, 3, 4]})
    out, steps = clean_xy(df, "group", "value")
    assert len(out) == 2
    assert any("Dropped 2 row" in s for s in steps)


def test_clean_xy_coerces_non_numeric_outcome():
    df = pd.DataFrame({"group": ["A", "B", "A"], "value": ["1.5", "oops", "2.0"]})
    out, steps = clean_xy(df, "group", "value")
    assert len(out) == 2
    assert any("Coerced column 'value'" in s for s in steps)


def test_clean_xy_errors_on_missing_column():
    df = pd.DataFrame({"group": ["A", "B"], "value": [1, 2]})
    with pytest.raises(ValueError, match="not found"):
        clean_xy(df, "group", "missing")
