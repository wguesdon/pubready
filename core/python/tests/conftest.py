"""Use a non-interactive matplotlib backend so recipe tests never open a window."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
    """Close every figure a test opened, so a long suite does not hold them all."""
    yield
    plt.close("all")
