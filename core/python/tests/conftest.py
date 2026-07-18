"""Use a non-interactive matplotlib backend so recipe tests never open a window."""
import matplotlib

matplotlib.use("Agg")
