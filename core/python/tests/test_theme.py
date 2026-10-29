"""The shared label and axis helpers. The R suite (test-theme.R) expects the
same strings."""
import pytest

from pubready.theme import (category_width, longest_line, plain_number_labels, wrap_label,
                            wrap_labels)


def test_wrap_label_breaks_a_long_label_at_the_spaces_as_the_r_engine_does():
    assert wrap_label("Vehicle control (DMSO 0.1%)") == "Vehicle\ncontrol\n(DMSO 0.1%)"
    assert wrap_label("Compound 7b, 10 uM (72 h)") == "Compound 7b,\n10 uM (72 h)"
    assert wrap_label("Overexpression") == "Overexpression"
    assert wrap_label("WT") == "WT"
    assert wrap_labels(["KO + rescue", "High dose + adjuvant"]) == ["KO + rescue",
                                                                    "High dose +\nadjuvant"]
    assert longest_line("High dose +\nadjuvant") == 11


def test_plain_number_labels_prints_plain_numbers_with_a_thousands_separator():
    assert plain_number_labels([-2e5, 0, 2e5, 4e5]) == ["-200,000", "0", "200,000", "400,000"]
    assert plain_number_labels([0, 0.25, 0.5]) == ["0.00", "0.25", "0.50"]
    assert plain_number_labels([1000, 1500, 2000]) == ["1,000", "1,500", "2,000"]
    assert plain_number_labels([10, 20, float("nan")]) == ["10", "20", ""]


def test_category_width_grows_with_the_number_of_groups():
    assert category_width(4) == pytest.approx(5.4)
    assert category_width(8) == pytest.approx(9.8)
    assert category_width(2) == pytest.approx(4.2)


def test_category_width_widens_the_slots_for_a_long_unbreakable_word():
    assert category_width(8, wrap_labels(["WT", "Overexpression"])) == pytest.approx(
        0.085 * 14 * 8 + 1)


def test_rain_offsets_puts_the_raincloud_points_left_of_the_box_as_r_does():
    from pubready.theme import rain_offsets
    assert list(rain_offsets(3)) == pytest.approx([-0.1058359, -0.1516718, -0.0775078], abs=1e-6)
    off = rain_offsets(500)
    assert ((off >= -0.18) & (off <= -0.06)).all()
