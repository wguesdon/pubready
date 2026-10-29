"""The figkit command of an installed pubready.

``uv add`` or ``uv tool install`` puts this console script on the PATH. It runs
the bundled ``cli/figkit`` with bash, so an installed figkit and a cloned
``./cli/figkit`` are the same program and drive the same pinned Podman image.

An installed copy is not a Git checkout, so figkit cannot ask Git for the commit
of the code it mounts. uv records that commit when it installs from a Git URL
(the PEP 610 ``direct_url.json`` file), and the launcher passes it on in
``PUBREADY_SOURCE_COMMIT`` so every bundle still names the exact code.
"""
import json
import os
import shutil
import sys
from importlib import metadata
from pathlib import Path

DIST_NAME = "pubready"


def tree_root():
    """Find the folder that holds cli/figkit and the files it mounts.

    A wheel carries the repository files under ``pubready_cli/tree``. An
    editable install has no copy, so the source checkout two levels above this
    file holds them.

    Returns:
        pathlib.Path: The folder.

    Raises:
        SystemExit: If neither location holds cli/figkit.
    """
    here = Path(__file__).resolve().parent
    for root in (here / "tree", here.parents[1]):
        if (root / "cli" / "figkit").is_file():
            return root
    raise SystemExit("figkit: the pubready files are missing from this install; reinstall pubready.")


def source_provenance():
    """Read the commit that uv installed from the direct_url.json record.

    Returns:
        tuple: The commit and the dirty state, as the strings that figkit
        records. A Git URL install gives its commit and ``false``; any other
        install gives ``unknown`` and ``not_recorded``.
    """
    try:
        raw = metadata.distribution(DIST_NAME).read_text("direct_url.json")
    except metadata.PackageNotFoundError:
        raw = None
    if not raw:
        return "unknown", "not_recorded"
    commit = (json.loads(raw).get("vcs_info") or {}).get("commit_id")
    if commit:
        return commit, "false"
    return "unknown", "not_recorded"


def build_env(base_env, commit, dirty):
    """Add the source provenance to the environment of figkit.

    A value that the caller already set wins, so a maintainer can override it.

    Args:
        base_env: The current environment.
        commit: The installed commit.
        dirty: The dirty state of the installed code.

    Returns:
        dict: The environment for figkit.
    """
    env = dict(base_env)
    env.setdefault("PUBREADY_SOURCE_COMMIT", commit)
    env.setdefault("PUBREADY_SOURCE_DIRTY", dirty)
    return env


def main(argv=None):
    """Run the bundled cli/figkit with the given arguments.

    Args:
        argv: The arguments; the command line when None.

    Raises:
        SystemExit: If bash is not available.
    """
    args = sys.argv[1:] if argv is None else list(argv)
    bash = shutil.which("bash")
    if bash is None:
        raise SystemExit("figkit needs bash. On Windows, run it inside WSL2.")
    script = tree_root() / "cli" / "figkit"
    commit, dirty = source_provenance()
    os.execve(bash, [bash, str(script), *args], build_env(os.environ, commit, dirty))


if __name__ == "__main__":
    main()
