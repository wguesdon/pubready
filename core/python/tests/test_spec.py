from types import SimpleNamespace

from pubready.spec import default_appearance, spec_from_args


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


def test_spec_from_args_carries_theme():
    a = SimpleNamespace(recipe="two_group_compare", data="d.csv",
                        x="g", y="v", theme="prism")
    assert spec_from_args(a)["appearance"]["theme"] == "prism"
    b = SimpleNamespace(recipe="two_group_compare", data="d.csv", x="g", y="v")
    assert spec_from_args(b)["appearance"]["theme"] == "pubready_house"


def test_spec_from_args_carries_new_data_and_cutoff_fields():
    a = SimpleNamespace(recipe="pca", data="d.csv", x=None, y=None, id="subj",
                        group="grp", label="gene", fc_cutoff=2.0, p_cutoff=0.01, top_n=20)
    s = spec_from_args(a)
    assert s["data"]["id"] == "subj"
    assert s["data"]["group"] == "grp"
    assert s["data"]["label"] == "gene"
    assert s["appearance"]["fc_cutoff"] == 2.0
    assert s["appearance"]["p_cutoff"] == 0.01
    assert s["appearance"]["top_n"] == 20
