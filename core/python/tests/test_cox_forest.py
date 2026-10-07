import pytest

from pubready.recipes.cox_forest import forest_stars


@pytest.mark.parametrize("p, stars", [
    (1e-9, "***"), (5e-4, "***"), (5e-3, "**"), (0.02, "*"), (0.05, ""), (0.3, ""),
])
def test_forest_stars_follow_ggforest(p, stars):
    # survminer::ggforest, which draws the R figure, stops at three stars.
    assert forest_stars(p) == stars
