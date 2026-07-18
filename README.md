# pubplot

An agentic skill to create reproducible, publication-ready figures from tidy
data, designed to work across several agent harnesses and LLMs.

## Setup

Two steps, once.

1. Build the container image. Needs Podman.

```bash
./cli/figkit build
```

2. Install the skill in your agent.

| Host | Install (run from the repo root) |
|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/pubplot` |
| Codex | `ln -s skills/codex/AGENTS.md AGENTS.md` (append if a root `AGENTS.md` exists) |
| opencode | add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`; for the `/pubplot` command also `mkdir -p .opencode/command && ln -s ../../skills/opencode/command/pubplot.md .opencode/command/pubplot.md` |

Open the agent in this folder and it can make figures.

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
