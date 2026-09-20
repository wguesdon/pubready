import pandas as pd

from pubready.bundle import dataframe_markdown


def test_dataframe_markdown_escapes_table_values():
    frame = pd.DataFrame({"name": ["a|b"], "value": [2]})

    table = dataframe_markdown(frame)

    assert "a\\|b" in table
    assert "| name | value |" in table
