"""Check the container's fallback glyphs through the actual renderer."""

import warnings

import matplotlib.pyplot as plt

from pubready import theme


def test_cjk_emoji_and_superscript_plus_render_without_missing_glyphs():
    """The default figure font stack must cover the reported missing glyphs."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig, ax = plt.subplots()
        ax.set_title("患者 治疗 🧬 😀 CD4⁺")
        fig.canvas.draw()
    missing = [str(w.message) for w in caught if "Glyph" in str(w.message) and "missing" in str(w.message)]
    assert missing == []


def test_standalone_script_uses_the_same_fonts(tmp_path):
    """A fresh script must render Unicode without importing the engine theme."""
    import subprocess
    import sys

    script = tmp_path / "font_check.py"
    script.write_text("import matplotlib.pyplot as plt\n"
                      "fig, ax = plt.subplots()\n"
                      "ax.set_title('患者 治疗 🧬 😀 CD4⁺')\n"
                      "fig.canvas.draw()\n")
    result = subprocess.run([sys.executable, str(script)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    assert "missing from font" not in result.stderr
