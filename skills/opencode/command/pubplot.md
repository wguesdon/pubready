---
description: Make a publication-ready figure from tidy data with the stats on top (pubplot).
agent: build
---

The user wants a publication figure from tidy data. Use pubplot. Follow the
adapter in `skills/opencode/AGENTS.md`: choose the test with the user, then call
`./cli/figkit`. Do not hand-write ggplot or matplotlib.

Data file: $ARGUMENTS

Start by inspecting the columns, then walk the design with the user before drawing:

`!./cli/figkit inspect --data $ARGUMENTS`

Then pick the recipe (`two_group_compare`, `multi_group_compare`, `factorial_anova`,
`survival_km`, `cox_forest`, `heatmap`, `correlation`, `correlation_heatmap`,
`upset`, `volcano`, `enrichment_dot`, `paired_compare`, `proportions`, `pca`),
confirm the columns and the test with the user, and run `./cli/figkit plot`. Leave
`--test auto` unless the user forces a test; add `--theme prism` for the GraphPad
Prism look, or `--geom raincloud|bar` on the comparison recipes. When it finishes,
point them at the bundle folder under `pubplot_output/`.
