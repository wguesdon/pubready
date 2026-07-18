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

## Documentation

- [docs/how_it_works.md](docs/how_it_works.md) — the workflow, the recipes, the bundle, and how to add one.
- [example/README.md](example/README.md) — worked examples for every recipe.
- [CHANGELOG.md](CHANGELOG.md) — what changed in each release.
