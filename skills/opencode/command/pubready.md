---
description: Make a publication-ready figure from tidy data with the stats on top (pubready).
agent: build
---

The user wants a publication figure from tidy data. Use pubready. Follow the
adapter in `skills/opencode/AGENTS.md`: choose the test with the user, then call
`./cli/figkit`. Do not hand-write ggplot or matplotlib.
After a uv install there is no `./cli/figkit`. Call `figkit` in its place, in
the commands below too.

Data file: $ARGUMENTS

Start by inspecting the columns, then walk the design with the user before drawing:

`!./cli/figkit inspect --data $ARGUMENTS`

Run `./cli/figkit advise --data FILE --format json` with the known plot flags
and `--answer KEY=VALUE` answers. Ask the scientist each unanswered question.
Never supply a scientist answer yourself. Repeat until the result is complete,
then follow `next_step`. Resolve conflicts before a plot. Exit code 5 means
that no recipe supports the design. Point the scientist to the final bundle.
