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

Then pick the recipe (`two_group_compare`, `multi_group_compare`, `factorial_anova`,
`survival_km`, `cox_forest`, `heatmap`, `correlation`, `correlation_heatmap`,
`upset`, `volcano`, `enrichment_dot`, `paired_compare`, `proportions`, `pca`,
`spider_response`, `waterfall`, `dose_response`, `roc_curve`, `bland_altman`,
`forest_meta`, `timecourse`, `swimmer`, `estimation_plot`, `enzyme_kinetics`,
`standard_curve`, `cumulative_incidence`, `calibration`, `sankey`, `ma_plot`,
`manhattan`, `feature_panel`, `embedding`, `distribution`, `composition`),
confirm the columns and the test with the user, and run `./cli/figkit plot`. Leave
`--test auto` unless the user forces a test; add `--theme prism` for the GraphPad
Prism look, or `--geom raincloud|bar` on the comparison recipes. When it finishes,
point them at the bundle folder under `pubready_output/`.
