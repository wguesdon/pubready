# pubplot

An agentic skill to create reproducible, publication-ready figures from tidy
data, designed to work across several agent harnesses and LLMs.

## Setup

Everything runs inside a pinned Podman container, so the only things you install
on the host are **Podman** and **git**. The `./cli/figkit` command is a bash
script, so you also need a POSIX shell: native on Linux and macOS, and WSL2 on
Windows.

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
manages. Clone pubplot somewhere under your home directory so the container can
mount it.

**Windows**

Use WSL2 with Ubuntu and treat it as Linux. In an elevated PowerShell:

```powershell
wsl --install -d Ubuntu
```

Reboot, open the Ubuntu shell, and install Podman with the Ubuntu commands above.
Clone and run pubplot from inside the WSL filesystem (your Linux home), not from a
`/mnt/c/...` Windows path.

### 2. Build the image

Clone this repo, then from its root build the image. The first build downloads
the R and Python stacks and takes a few minutes; after that it is cached.

```bash
git clone https://github.com/wguesdon/pubplot.git
cd pubplot
./cli/figkit build
```

### 3. Install the skill in your agent

| Host | Install (run from the repo root) |
|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/pubplot` |
| Codex | `ln -s skills/codex/AGENTS.md AGENTS.md` (append if a root `AGENTS.md` exists) |
| opencode | add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`; for the `/pubplot` command also `mkdir -p .opencode/command && ln -s ../../skills/opencode/command/pubplot.md .opencode/command/pubplot.md` |

Open the agent in this folder and it can make figures. Confirm the install with a
quick render:

```bash
./cli/figkit plot --recipe two_group_compare --data example/tumor_volume.csv \
  --x group --y volume
```

## Tested models

pubplot is host and model agnostic: any agent that reads the adapter and can run
`./cli/figkit` will drive it. The combinations below have produced a correct
figure end to end.

- **Claude Code** with its default model.
- **opencode** with the models below, run as `opencode run -m <id>`:

| Model (`opencode -m` id) | Result |
|---|---|
| `opencode/claude-sonnet-4-5` | works, including the `/pubplot` command |
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
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release.
