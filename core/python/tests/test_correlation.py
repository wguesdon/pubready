from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pubready.recipes.correlation import _resolve, recipe_correlation
from pubready.spec import spec_from_args


def test_resolve_picks_method_from_normality():
    assert _resolve("auto", True, True)["method"] == "pearson"
    assert _resolve("auto", True, False)["method"] == "spearman"
    assert _resolve("auto", False, True)["method"] == "spearman"
    assert _resolve("spearman", True, True)["method"] == "spearman"
    with pytest.raises(ValueError, match="unknown correlation method"):
        _resolve("nonsense", True, True)


def test_recipe_correlation_runs_end_to_end():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 10, 40)
    y = 2 * x + rng.normal(0, 1, 40)
    df = pd.DataFrame({"x": x, "y": y})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="y", test="pearson")
    spec = spec_from_args(a)

    res = recipe_correlation(df, spec)

    row = res["stats"].iloc[0]
    assert row["method"] == "Pearson correlation"
    assert row["n"] == 40
    assert row["estimate"] > 0.9
    assert row["p_value"] < 1e-6


def test_recipe_correlation_needs_present_columns():
    df = pd.DataFrame({"x": [1, 2, 3]})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="missing", test="auto")
    spec = spec_from_args(a)
    with pytest.raises(ValueError, match="not found"):
        recipe_correlation(df, spec)


def test_the_correlation_subtitle_names_the_method_the_coefficient_and_p():
    from pubready.recipes.correlation import cor_label
    assert (cor_label("Pearson correlation", "r", 0.7123, 2.1e-10)
            == "Pearson correlation, r = 0.71, P = 2.1e-10")
    assert (cor_label("Spearman rank correlation", "rho", -0.456, 0.0312)
            == "Spearman rank correlation, rho = -0.46, P = 0.031")
    x = np.arange(1, 21, dtype=float)
    df = pd.DataFrame({"x": x, "y": x + np.sin(x)})
    a = SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="y", test="kendall")
    res = recipe_correlation(df, spec_from_args(a))
    assert res["fig"].axes[0].get_title().startswith("Kendall rank correlation, tau = ")


def test_the_band_is_the_confidence_interval_of_the_fitted_mean():
    """lm_band gives the 95% CI of the fitted mean, the band of geom_smooth(method = "lm")."""
    import statsmodels.api as sm

    from pubready.recipes.correlation import lm_band
    rng = np.random.default_rng(4)
    x = rng.normal(size=30)
    y = 2 * x + rng.normal(size=30)
    grid, fit, lo, hi = lm_band(x, y)
    pred = (sm.OLS(y, sm.add_constant(x)).fit()
            .get_prediction(sm.add_constant(grid)).summary_frame(alpha=0.05))
    assert len(grid) == 80 and grid[0] == x.min() and grid[-1] == x.max()
    assert fit == pytest.approx(pred["mean"].to_numpy(), rel=1e-10)
    assert lo == pytest.approx(pred["mean_ci_lower"].to_numpy(), rel=1e-10)
    assert hi == pytest.approx(pred["mean_ci_upper"].to_numpy(), rel=1e-10)


def test_the_emitted_script_draws_the_same_band_and_runs(tmp_path):
    import subprocess
    import sys
    x = np.linspace(0, 10, 20)
    df = pd.DataFrame({"x": x, "y": 2 * x + np.sin(x)})
    spec = spec_from_args(SimpleNamespace(recipe="correlation", data="toy.csv", x="x", y="y",
                                          test="pearson"))
    script = recipe_correlation(df, spec)["build_script"]("input.csv", "figure")
    assert "def lm_band(" in script and "regplot" not in script
    df.to_csv(tmp_path / "input.csv", index=False)
    (tmp_path / "script.py").write_text(script)
    run = subprocess.run([sys.executable, "script.py"], cwd=tmp_path, capture_output=True,
                         text=True, check=False)
    assert run.returncode == 0, run.stderr
    assert (tmp_path / "figure.pdf").exists()


def test_correlation_statistic_and_ci_match_stats_cor_test():
    from pubready.recipes.correlation import fisher_ci, kendall_statistic, spearman_s
    # mtcars wt and mpg: the conf.int of stats::cor.test(method = "pearson").
    low, high = fisher_ci(-0.8676593765172278, 32)
    assert low == pytest.approx(-0.9338264132849942, rel=1e-9)
    assert high == pytest.approx(-0.7440871964601128, rel=1e-9)
    # airquality Temp and Ozone: the S of stats::cor.test(method = "spearman").
    assert spearman_s(0.7740429554613013, 116) == pytest.approx(58778.20599585169, rel=1e-9)
    # Without a tie and below 50 pairs, R reports T, the count of concordant pairs.
    assert kendall_statistic([1, 2, 3, 4, 5], [1, 3, 2, 5, 4], 0.6) == 8


def test_a_rank_correlation_scatter_gets_no_fitted_line():
    df = pd.DataFrame({"a": np.arange(1, 31), "b": np.arange(1, 31) ** 2.0})
    base = dict(recipe="correlation", data="c.csv", x="a", y="b")
    sp = recipe_correlation(df, spec_from_args(SimpleNamespace(**base, test="spearman")))
    assert len(sp["fig"].axes[0].lines) == 0
    assert "No line is drawn" in sp["methods"]
    pe = recipe_correlation(df, spec_from_args(SimpleNamespace(**base, test="pearson")))
    assert len(pe["fig"].axes[0].lines) == 1
