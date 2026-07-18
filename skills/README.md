# skills

Thin adapters that teach each host tool how to call `figkit` and where the
reference docs live. The reasoning and the code stay in the shared core, so the
three adapters do not drift.

- `claude-code/SKILL.md` — Claude Code Agent Skill (`name` + `description`
  frontmatter, triggers on figure requests).
- `codex/AGENTS.md` — Codex instructions, read from the nearest `AGENTS.md`.
- `opencode/` — opencode `AGENTS.md` instructions plus a `/pubplot` custom
  command in `command/pubplot.md`.

Each adapter carries the same short operating guide: the `figkit` contract, the
workflow, the recipe table, and how the `auto` test is chosen. The detail lives in
`docs/how_it_works.md`, `example/README.md`, and `reference/`. Update the adapters
together when the `figkit` surface changes.

## Install

| Host | Activate |
|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/pubplot` |
| Codex | `ln -s skills/codex/AGENTS.md AGENTS.md` (append if a root `AGENTS.md` exists) |
| opencode | add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`, and link `command/pubplot.md` into `.opencode/command/` |

See each adapter's own file for the full install note.
