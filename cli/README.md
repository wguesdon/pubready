# cli

The stable command surface the agent calls: `figkit`.

The agent does not render figures directly. It calls `figkit`, which runs the
right recipe inside the container. Examples:

    figkit inspect data.csv
    figkit plot --engine r --recipe two_group_compare --x group --y value --test auto

Keeping this surface stable is what lets the LLM stay thin.
