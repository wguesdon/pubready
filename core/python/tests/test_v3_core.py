"""Tests of the area B fixes of the v3 findings: the read rule, the bundle, the
spec and the inspect report. The R suite (tests/testthat/test-v3-core.R) holds
the same cases with the same expected values."""
import sys

import numpy as np
import pytest

from pubready.io import read_tidy, read_tidy_logged


def write_tmp_csv(tmp_path, lines):
    """Write lines to a CSV file in the test folder and return its path."""
    path = tmp_path / "data.csv"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_decimal_comma_stops_and_names_the_column(tmp_path):
    path = write_tmp_csv(tmp_path, ["group,value", 'control,"46,2"', "treated,51.0"])
    with pytest.raises(ValueError, match="column 'value' holds '46,2', a number with a comma"):
        read_tidy(path)


def test_thousands_mark_stops_and_names_the_column(tmp_path):
    path = write_tmp_csv(tmp_path, ["group,value", 'control,"4,380"', "treated,4570"])
    with pytest.raises(ValueError, match="column 'value' holds '4,380', a number with a comma"):
        read_tidy(path)


def test_text_value_with_a_comma_stays_text(tmp_path):
    path = write_tmp_csv(tmp_path, ["id,state", '1,"edema, no diuretics"', "2,no edema"])
    df = read_tidy(path)
    assert list(df["state"]) == ["edema, no diuretics", "no edema"]


def test_only_an_empty_cell_or_na_reads_as_missing(tmp_path):
    path = write_tmp_csv(tmp_path, ["arm,response,score", "a,None,1", "a,N/A,", "b,NULL,NA",
                                    "b,nan,4", "c,,5"])
    df, steps = read_tidy_logged(path)
    assert list(df["response"][:4]) == ["None", "N/A", "NULL", "nan"]
    assert df["response"].isna().tolist() == [False, False, False, False, True]
    assert df["score"].isna().tolist() == [False, True, True, False, False]
    assert list(df["score"].dropna()) == [1.0, 4.0, 5.0]
    assert steps == ["Read 3 cell(s) that hold \"NA\" or nothing as missing: "
                     "'response' (1), 'score' (2)."]


def test_spaces_around_a_value_and_a_column_name_are_removed(tmp_path):
    path = write_tmp_csv(tmp_path, ["group , value", "control ,52.2", " control, 49.9",
                                    "treated,48.6"])
    df, steps = read_tidy_logged(path)
    assert list(df.columns) == ["group", "value"]
    assert list(df["group"]) == ["control", "control", "treated"]
    assert list(df["value"]) == [52.2, 49.9, 48.6]
    assert steps == ["Removed the spaces around 2 column name(s).",
                     "Removed the spaces around 3 value(s): 'group' (2), 'value' (1)."]


def test_one_text_cell_keeps_text_and_logical_tokens_read_as_logical(tmp_path):
    path = write_tmp_csv(tmp_path, ["a,b,c,d", "1,TRUE,T,1", "nan,false,F,2", "3,,T,3"])
    df = read_tidy(path)
    assert list(df["a"]) == ["1", "nan", "3"]
    assert list(df["b"][:2]) == [True, False] and df["b"].isna().tolist() == [False, False, True]
    assert list(df["c"]) == ["T", "F", "T"]
    assert list(df["d"]) == [1, 2, 3] and str(df["d"].dtype) == "int64"


def test_subnormal_p_value_keeps_its_value(tmp_path):
    path = write_tmp_csv(tmp_path, ["snp,p", "rs1,1e-320", "rs2,0.5", "rs3,2.5E-10"])
    df = read_tidy(path)
    assert df["p"][0] == 1e-320
    assert 0 < df["p"][0] < sys.float_info.min
    assert df["p"][2] == 2.5e-10
    assert np.isfinite(df["p"]).all()
