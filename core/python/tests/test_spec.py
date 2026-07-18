from types import SimpleNamespace

from pubplot.spec import default_appearance, spec_from_args


def test_spec_from_args_builds_python_spec_with_defaults():
    a = SimpleNamespace(recipe="two_group_compare", data="/p/tumor_volume.csv",
                        x="group", y="volume", test="auto")
    spec = spec_from_args(a)
    assert spec["engine"] == "python"
    assert spec["recipe"] == "two_group_compare"
    assert spec["data"]["file"] == "tumor_volume.csv"
    assert spec["data"]["x"] == "group"
    assert spec["test"]["method"] == "auto"
    assert spec["test"]["paired"] is False
    assert spec["appearance"]["geom"] == "box"


def test_spec_from_args_splits_palette_and_covariates():
    a = SimpleNamespace(recipe="cox_forest", data="trial.csv", x=None, y=None,
                        palette="#3B6DB3, #C1432B", covariates="age, sex, stage")
    spec = spec_from_args(a)
    assert spec["appearance"]["palette"] == ["#3B6DB3", "#C1432B"]
    assert spec["data"]["covariates"] == ["age", "sex", "stage"]


def test_default_appearance_shape():
    ap = default_appearance()
    assert ap["bracket"] == {"show": True, "label": "p.signif"}
    assert ap["y_limits"] == [None, None]
