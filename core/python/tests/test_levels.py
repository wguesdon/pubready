import numpy as np
import pandas as pd
import pytest

from pubready.util import level_labels, level_order


def test_level_labels_writes_whole_floats_as_integers():
    s = pd.Series([80.0, 160.0, np.nan, 0.5, 2.5e-9])
    out = level_labels(s)
    assert list(out[:2]) == ["80", "160"]
    assert np.isnan(out[2])
    assert list(out[3:]) == ["0.5", "2.5e-09"]


def test_level_order_sorts_numbers_as_numbers():
    assert level_order(["0", "160", "235", "310", "80"]) == ["0", "80", "160", "235", "310"]
    assert level_order(pd.Series([80.0, 0.0, np.nan, 160.0])) == ["0", "80", "160"]


def test_level_order_sorts_text_by_code_point():
    assert level_order(["control", "Drug", "non-USA", "USA"]) == ["Drug", "USA", "control", "non-USA"]


def test_level_order_follows_levels_and_drops_unused_levels():
    assert level_order(["wt", "ko", "het"], ["wt", "het", "ko", "rescue"]) == ["wt", "het", "ko"]


def test_level_order_stops_on_a_value_that_levels_does_not_list():
    with pytest.raises(ValueError, match="does not list the group\\(s\\) 'ko'"):
        level_order(["wt", "ko"], ["wt"])
