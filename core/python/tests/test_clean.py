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


def test_coerce_numeric_columns_keeps_a_mostly_numeric_column():
    import pandas as pd
    from pubready.clean import coerce_numeric_columns
    df = pd.DataFrame({"tnf": ["1.5", "2.0", "<0.1", "3.2", "4.1"],
                       "group": ["a", "a", "b", "b", "b"],
                       "il6": [1.0, 2.0, 3.0, 4.0, 5.0]})
    out, numeric, text, steps = coerce_numeric_columns(df)
    assert numeric == ["tnf", "il6"]
    assert text == ["group"]
    assert out["tnf"].isna().sum() == 1
    assert steps == ["Coerced column 'tnf' to numeric; 1 value(s) became NA, for example '<0.1'."]


def test_coerce_numeric_columns_reads_a_logical_column_as_text():
    import pandas as pd
    from pubready.clean import coerce_numeric_columns
    df = pd.DataFrame({"flag": [True, False, True], "obj": pd.Series([True, None, False], dtype=object),
                       "v": [1.0, 2.0, 3.0]})
    _, numeric, text, _ = coerce_numeric_columns(df)
    assert numeric == ["v"]
    assert text == ["flag", "obj"]
