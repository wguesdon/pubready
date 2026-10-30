"""Tests of the area B fixes of the v3 findings: the read rule, the bundle, the
spec and the inspect report. The R suite (tests/testthat/test-v3-core.R) holds
the same cases with the same expected values."""
import json
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.bundle import (_markdown_value, dataframe_markdown, manifest_arguments,
                             reproduce_md, write_json,
                             write_stats_csv)
from pubready.io import read_tidy, read_tidy_logged
from pubready.qc import qc_layout, qc_normality
from pubready.run import inspect_data, print_inspect, run_entry
from pubready.spec import spec_from_args


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


def test_manifest_writes_null_for_nan_and_infinity(tmp_path):
    path = str(tmp_path / "manifest.json")
    write_json({"a": float("nan"), "b": float("inf"), "c": None,
                "d": [np.float64(1.5), np.float64("-inf")], "e": np.int64(3)}, path)
    txt = open(path).read()
    assert "NaN" not in txt and "Infinity" not in txt
    assert json.loads(txt) == {"a": None, "b": None, "c": None, "d": [1.5, None], "e": 3}


def test_manifest_keeps_the_full_precision_of_each_number(tmp_path):
    path = str(tmp_path / "manifest.json")
    write_json({"intercept": 0.00508612345678912, "p": 1.5e-12}, path)
    back = json.loads(open(path).read())
    assert back == {"intercept": 0.00508612345678912, "p": 1.5e-12}


#: The arguments block of the manifest, the same list in both engines (F101).
MANIFEST_ARG_NAMES = [
    "x", "y", "fill", "facet", "time", "event", "covariates", "annotation", "id",
    "group", "label", "positive", "se", "ci_low", "ci_high", "measure", "markers",
    "ref", "cause", "levels", "visit_levels", "features", "test", "paired",
    "p_adjust", "scale", "cluster", "fc_cutoff", "p_cutoff", "top_n",
    "pd_threshold", "pr_threshold", "geom", "theme", "xlab", "ylab", "title",
    "palette"]


def test_manifest_records_every_plot_option(tmp_path):
    a = SimpleNamespace(recipe="sankey", data="d.csv", x="week", y="response", id="patient",
                        levels="CR|PR", covariates="age", xlab="Week", test="auto")
    args = manifest_arguments(spec_from_args(a))
    assert list(args) == MANIFEST_ARG_NAMES
    path = str(tmp_path / "manifest.json")
    write_json({"arguments": args}, path)
    back = json.loads(open(path).read())["arguments"]
    assert back["id"] == "patient"
    assert back["levels"] == ["CR", "PR"]
    assert back["covariates"] == ["age"]
    assert back["fill"] is None
    assert back["xlab"] == "Week"
    assert back["paired"] is False
    assert back["top_n"] == 15


def test_python_reproduce_md_names_the_image_digest():
    txt = reproduce_md({"image": "img", "image_id": "abc", "digest": ""}, "c0ffee", "script.py")
    assert "- Image digest: `n/a (local build)`" in txt
    txt = reproduce_md({"image": "img", "image_id": "abc", "digest": "sha256:1"}, "c0ffee",
                       "script.py")
    assert "- Image digest: `sha256:1`" in txt


#: One CSV form for a stats table in both engines (F126).
STATS_CSV_EXPECTED = [
    "name,n,p,ok,d,flag",
    "a,3,0.5,TRUE,-Inf,TRUE",
    "b,NA,NA,FALSE,2,FALSE",
    "c,5,Inf,NA,NA,TRUE"]


def test_stats_csv_writes_na_booleans_infinity_and_integers_in_one_form(tmp_path):
    df = pd.DataFrame({"name": ["a", "b", "c"], "n": [3.0, np.nan, 5.0],
                       "p": [0.5, np.nan, np.inf], "ok": [True, False, None],
                       "d": [-np.inf, 2.0, np.nan], "flag": [True, False, True]})
    path = str(tmp_path / "stats.csv")
    write_stats_csv(df, path)
    assert open(path).read().splitlines() == STATS_CSV_EXPECTED


def test_report_table_writes_the_r_forms():
    assert _markdown_value(float("nan")) == ""
    assert _markdown_value(np.bool_(True)) == "TRUE"
    assert _markdown_value(float("-inf")) == "-Inf"


def test_one_panel_qc_figure_wraps_its_title_and_subtitle():
    lay = qc_layout("differences between the methods (new_assay - reference_assay)", 1)
    assert lay["subtitle"] == ["Points near the line support normality;",
                               "mild, judged deviation is usually",
                               "acceptable."]
    assert lay["title"] == ["Normality QC:", "differences between the", "methods (new_assay -",
                            "reference_assay)"]
    assert all(len(s) <= 23 for s in lay["title"])
    assert lay["width"] == pytest.approx(3.3)
    assert lay["height"] == pytest.approx(3.9 + 0.25 * 3 + 0.18 * 2)
    three = qc_layout("group values", 3)
    assert len(three["subtitle"]) == 1
    assert three["title"] == ["Normality QC: group values"]
    fig, _ = qc_normality({"quantity": "x", "panels": {"a": np.arange(10.0)}})
    assert tuple(fig.get_size_inches()) == pytest.approx((3.3, 3.9 + 0.18 * 2))


def test_report_table_writes_a_category_with_na():
    df = pd.DataFrame({"group": pd.Categorical(["a", None, "b|c"])})
    assert dataframe_markdown(df).splitlines() == [
        "| group |", "| --- |", "| a |", "|  |", "| b\\|c |"]


def test_report_table_writes_whole_numbers_as_r_does():
    df = pd.DataFrame({"n": [1.0, np.nan, 3.0], "p": [0.5, np.inf, np.nan]})
    assert dataframe_markdown(df).splitlines() == [
        "| n | p |", "| --- | --- |", "| 1 | 0.5 |", "|  | Inf |", "| 3 |  |"]


def _raise_data_error():
    """Stand in for an entry main that meets a data error."""
    raise ValueError("column 'value' not found in data")


def _raise_defect():
    """Stand in for an entry main that meets a defect of pubready."""
    raise KeyError("value")


def test_entry_prints_a_data_error_without_a_traceback(capsys):
    with pytest.raises(SystemExit) as stop:
        run_entry(_raise_data_error)
    assert stop.value.code == 1
    err = capsys.readouterr().err
    assert err == "Error: column 'value' not found in data\n"
    with pytest.raises(KeyError):
        run_entry(_raise_defect)


#: The inspect report, the same words in both engines (F72).
INSPECT_CSV = ["extra,group,ID,sex", "0.7,1,1,f", "-1.6,1,2,m", "-0.2,1,3,f", "-1.2,1,4,m",
               "-0.1,1,5,f", "1.9,2,1,f", "0.8,2,2,m", "1.1,2,3,f", "0.1,2,4,m", "-0.1,2,5,f"]
INSPECT_EXPECTED = [
    "Rows: 10   Columns: 4",
    "",
    "  extra                numeric   missing=0  unique=9  [-1.6 .. 1.9, mean 0.14]",
    "  group                numeric   missing=0  unique=2  levels: \"1\" (5), \"2\" (5)"
    "  [1 .. 2, mean 1.5]",
    "  ID                   numeric   missing=0  unique=5  levels: \"1\" (2), \"2\" (2), "
    "\"3\" (2), \"4\" (2), \"5\" (2)  [1 .. 5, mean 3]",
    "  sex                  text      missing=0  unique=2  levels: \"f\" (6), \"m\" (4)",
    "",
    "Candidate grouping columns (x):  group (numeric code), ID (numeric code), sex",
    "Candidate outcome columns (y):   extra",
    "Possible subject id columns:     ID (repeats across the levels of group)"]


def test_inspect_gives_level_counts_numeric_codes_and_a_subject_id(tmp_path, capsys):
    info = inspect_data(read_tidy(write_tmp_csv(tmp_path, INSPECT_CSV)))
    print_inspect(info)
    assert capsys.readouterr().out.splitlines() == INSPECT_EXPECTED
    assert info["id_candidates"] == [{"column": "ID", "across": "group", "n_subjects": 5,
                                      "complete": True}]
    assert info["group_candidates"] == ["group", "ID", "sex"]
    assert info["outcome_candidates"] == ["extra"]
