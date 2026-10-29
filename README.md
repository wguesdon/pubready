# pubready

An agentic skill to create reproducible, publication-ready figures from tidy
data, designed to work across several agent harnesses and LLMs.

## Setup

Everything runs inside a pinned Podman container, so the only things you install
on the host are **Podman** and **git**. The `./cli/figkit` command is a bash
script, so you also need a POSIX shell: native on Linux and macOS, and WSL2 on
Windows.

You can get the `figkit` command in three ways. Each way runs the same pinned
image, and each bundle records the commit of the code that drew the figure.

| Way | Command | Use it when |
|---|---|---|
| uv tool | `uv tool install "pubready @ git+https://github.com/wguesdon/pubready@main"` | You want `figkit` in each folder of your machine. |
| uv project | `uv add "pubready @ git+https://github.com/wguesdon/pubready@main"` | You want pubready pinned in the `pyproject.toml` of an analysis project. |
| git clone | `git clone https://github.com/wguesdon/pubready.git` | You want to change pubready, or to run its tests and examples. |

Step 2 below gives the details of each way.

### AWS SageMaker notebooks

Use the notebook terminal to run the command line workflow. The container supports
both `amd64` instances and `arm64` Graviton instances.

Run this script if you also want Quarto in the notebook environment:

```bash
./scripts/install_quarto.sh
export PATH="$HOME/.local/bin:$PATH"
quarto --version
```

The script installs Quarto 1.11.4 under `$HOME/.local`. It does not require root
access. The script verifies the downloaded archive before it installs Quarto.

Build the image with rootless Podman, and then run `figkit` from the terminal:

```bash
./cli/figkit build
./cli/figkit plot --recipe two_group_compare \
  --data example/tumor_volume.csv --x group --y volume
```

Each generated bundle contains `report.qmd`, `report.html`, and
`report_provenance.json`. The HTML report contains its images and other resources.
The report command stops when the repository contains an untracked file.

### 1. Install Podman

**Ubuntu / Debian**

```bash
sudo apt update && sudo apt install -y podman git
```

Rootless Podman works out of the box. On an older release you may want a newer
Podman; see the [official install guide](https://podman.io/docs/installation).

**macOS** (Intel or Apple Silicon)

```bash
brew install podman            # or install Podman Desktop from podman.io
podman machine init
podman machine start
```

Containers are Linux, so Podman runs a small Linux VM that `podman machine`
manages. Clone pubready somewhere under your home directory so the container can
mount it.

**Windows**

Use WSL2 with Ubuntu and treat it as Linux. In an elevated PowerShell:

```powershell
wsl --install -d Ubuntu
```

Reboot, open the Ubuntu shell, and install Podman with the Ubuntu commands above.
Clone and run pubready from inside the WSL filesystem (your Linux home), not from a
`/mnt/c/...` Windows path.

### 2. Get figkit and build the image

The first build downloads the R and Python stacks. Later builds use the cache.

#### With uv

The uv install needs `uv` and Python 3.9 or later. It does not install the R
engine or the Python engine on the host. It installs the files that `figkit`
mounts into the image, and `figkit build` builds the image from them.

Install `figkit` as a uv tool to use it in each folder of your machine:

```bash
uv tool install "pubready @ git+https://github.com/wguesdon/pubready@main"
figkit build
```

Add pubready to a uv project to pin it with the other dependencies of the
analysis. In the project, run `figkit` through `uv run`:

```bash
uv add "pubready @ git+https://github.com/wguesdon/pubready@main"
uv run figkit build
```

To pin a release, replace `main` with the tag of that release. Release 1.0.0 and
the releases before it do not have the uv install.

uv records the installed commit, and `figkit` writes that commit into each
bundle. An install from a local folder records no commit, so `figkit plot` stops
and asks for an install from the Git URL.

#### With git clone

Clone the repository, then build the image from its root:

```bash
git clone https://github.com/wguesdon/pubready.git
cd pubready
./cli/figkit build
```

#### Check the setup

Run `figkit doctor` after the build, or when a command reports a Podman problem.
It checks Podman, rootless mode, the macOS Podman machine, the image and the
folders that `figkit` mounts. For each failed check, it prints the fix. It does
not install anything, and it does not run `sudo`.

```bash
figkit doctor        # from a clone: ./cli/figkit doctor
```

### 3. Install the skill in your agent

Run the commands of a clone from the root of the repository. Run the commands of
a uv install from the folder where you start the agent. `figkit where` prints the
folder that holds the adapters of a uv install.

| Host | From a clone | From a uv install |
|---|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/pubready` | `mkdir -p ~/.claude/skills && ln -s "$(figkit where)/skills/claude-code" ~/.claude/skills/pubready` |
| Codex | `ln -s skills/codex/AGENTS.md AGENTS.md` (append if a root `AGENTS.md` exists) | `ln -s "$(figkit where)/skills/codex/AGENTS.md" AGENTS.md` |
| opencode | add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`; for the `/pubready` command also `mkdir -p .opencode/command && ln -s ../../skills/opencode/command/pubready.md .opencode/command/pubready.md` | add the path that `echo "$(figkit where)/skills/opencode/AGENTS.md"` prints to `instructions` in `opencode.json`; for the `/pubready` command also `mkdir -p .opencode/command && ln -s "$(figkit where)/skills/opencode/command/pubready.md" .opencode/command/pubready.md` |

If `figkit where` prints a different folder after `uv tool upgrade pubready`,
make the links again.

Confirm the install by running `figkit` directly once. From a clone:

```bash
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume
```

From a uv install, copy the example data into the current folder first.
`figkit` mounts only the current folder into the container, so the data file
must be in that folder or below it.

```bash
cp "$(figkit where)/example/tumor_volume.csv" .
figkit plot --recipe two_group_compare --data tumor_volume.csv --x group --y volume
```

### 4. Ask for a figure

Open the agent in this folder and describe the figure in plain language. The agent
picks the recipe and runs `figkit`; you do not call `figkit` yourself. How you
invoke pubready differs per host, because each one loads the skill differently.

**Claude Code** — the skill triggers from your request. Just ask:

> Use pubready to compare tumor volume between the two groups in `example/tumor_volume.csv`, with the right statistical test on top.

**Codex** — it reads the root `AGENTS.md` and follows it, so ask the same way. Run
Codex outside its sandbox so Podman can use a user namespace.

> Use pubready to compare tumor volume between the two groups in `example/tumor_volume.csv`.

**opencode** — use the `/pubready` slash command with a data file, or a plain prompt:

> /pubready example/tumor_volume.csv

Or drive it headless with any model. Pass `--auto` so opencode approves the
`figkit` calls without a prompt:

```bash
opencode run --auto -m openrouter/openai/gpt-oss-120b \
  "Using pubready, compare tumor volume between the two groups in example/tumor_volume.csv."
```

Each run writes a self-contained bundle under `pubready_output/`. The bundle
contains the figure, statistics, a runnable script, a manifest, and a Quarto report.

## Tested models

pubready is host and model agnostic: any agent that reads the adapter and can run
`./cli/figkit` will drive it. The combinations below have produced a correct
figure end to end.

- **Claude Code** with its default model.
- **opencode** with the models below, run headless as `opencode run --auto -m <id>`:

| Model (`opencode -m` id) | Result |
|---|---|
| `opencode/claude-sonnet-4-5` | works, including the `/pubready` command |
| `opencode/glm-5.2` | works |
| `opencode/kimi-k2.7-code` | works |
| `openrouter/qwen/qwen3-coder` | works |
| `opencode/deepseek-v4-pro` | works |
| `opencode/minimax-m3` | works |
| `opencode/grok-4.5` | works |
| `openrouter/openai/gpt-oss-120b` | works, but wants a directive prompt |

**Codex** reads the adapter and runs the correct `figkit` command, but producing
a figure needs Codex run outside its sandbox, since Podman needs a user namespace.

## Documentation

- [docs/how_it_works.md](docs/how_it_works.md) — the workflow, the recipes, the bundle, and how to add one.
- [example/README.md](example/README.md) — worked examples for every recipe.
- [benchmarks/README.md](benchmarks/README.md) lists the public datasets that check every recipe, in both engines, against an independent reference value.
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release.
