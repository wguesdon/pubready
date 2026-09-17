# opencode adapter

Two ways to wire pubready into opencode. Use either or both.

## Always-on instructions

opencode reads `AGENTS.md` as project instructions. Point it at the adapter with
the `instructions` array in `opencode.json` at the repo root:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "instructions": ["skills/opencode/AGENTS.md"]
}
```

Or link the adapter as the repo-root `AGENTS.md`:

```bash
ln -s skills/opencode/AGENTS.md AGENTS.md
```

If a root `AGENTS.md` already exists, add a line that references
`skills/opencode/AGENTS.md` rather than overwriting it.

## The /pubready command

`command/pubready.md` is an opencode custom command. Install it so opencode finds
it, then invoke it with a data file:

```bash
mkdir -p .opencode/command
ln -s ../../skills/opencode/command/pubready.md .opencode/command/pubready.md
```

```
/pubready example/tumor_volume.csv
```

The command inspects the columns first, then walks the design with you before
drawing. The reasoning and the code stay in the shared core, so this adapter does
not drift from the Claude Code and Codex ones.
