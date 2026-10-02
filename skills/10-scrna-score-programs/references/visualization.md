# Program-score visualization

The built-in figures preserve the useful plot families from the sc06 pathway, metabolism, and Hallmark scripts while removing project-specific metadata names and endothelial labels.

## Group heatmaps

Set `summarize_by` to include every grouping field needed for scientific summaries. Configure `visualization.group_heatmap.x` as the comparison field (usually condition) and `facet` as the population field (usually cell type). The renderer averages over any additional summary fields, such as sample, so including sample yields an equal-sample descriptive mean rather than weighting samples by cell count. Keep `scale: none` when score magnitudes within a method are meaningful; use `row_zscore` only to compare each signature's relative pattern, and label it accordingly.

Set `focus_levels` to selected values of the facet field to produce an additional focused heatmap. This replaces hard-coded `EC`, `Endocardial`, and `LEC` logic. The all-population heatmap is always retained.

## UMAP activity maps

UMAP figures are opt-in because many signatures and split levels can create very large outputs. Set `visualization.umap.enabled: true`, name an existing reduction, optionally set `split_by`, and choose `features_per_page`. The renderer uses one shared scale across split panels and paginates every retained signature; it never silently shows only a top subset.

UMAP localization and group heatmaps are descriptive. They do not test condition effects. Use the sample-aware guidance in `interpretation-and-inference.md` for formal inference.

## External figure references

When ScientificFigureLibrary is available, search it by plot family and data shape (for example `heatmap`, `single-cell`, `pathway activity`, or `faceted tile`) and inspect the selected template's code, license, and provenance before adaptation. Keep scoring, signature selection, thresholds, and complete-result coverage unchanged. Treat the library as an optional visual reference: the executor remains reproducible without it, and third-party templates are not copied into this skill.

Default exports are white-background PNG at 300 dpi. Preserve complete scientific tables independently of any display filtering or pagination.
