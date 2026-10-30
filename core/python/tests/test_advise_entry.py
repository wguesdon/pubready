"""Exercise the installed entry boundary on each approved benchmark route."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from test_advise_walk import CASES


@pytest.mark.parametrize("dataset,flags,answers,recipe,code", CASES)
def test_entry_benchmark_contract(dataset, flags, answers, recipe, code):
    """Verify exit codes and parseable JSON from the actual entry script."""
    entry = Path(__file__).resolve().parents[1] / "entry/advise.py"
    command = [sys.executable, str(entry), "--format", "json", "--data",
               str(Path(os.environ["PUBREADY_BENCHMARKS"]) / "data" / f"{dataset}.csv")]
    for key, value in flags.items():
        command.extend([f"--{key}", value])
    for answer in answers:
        command.extend(["--answer", answer])
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    assert completed.returncode == code, completed.stderr + completed.stdout
    result = json.loads(completed.stdout)
    assert result["recipe"] == recipe
    assert result["exit_code"] == code
