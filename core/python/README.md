# pubplot (Python engine)

The Python engine for [pubplot](../../README.md): recipe functions that turn
tidy data into a publication-ready figure with the chosen statistical test drawn
on top, mirroring the R engine one to one and writing the same reproducibility
bundle.

Not meant to be installed standalone. The figkit CLI loads this package from a
runtime mount (`PYTHONPATH`) inside the pinned container. This `pyproject.toml`
exists for versioning, dependency declaration, and running the test suite.

```bash
figkit test          # runs this suite (pytest) and the R suite in-container
```

Stack: pandas, numpy, scipy, matplotlib, seaborn, pingouin, statannotations,
statsmodels, scikit-posthocs, lifelines, PyComplexHeatmap.
