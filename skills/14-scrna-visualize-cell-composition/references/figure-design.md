# Figure design

The default figure family follows the repository's `paper_v1` contract: white background, quiet axes, deterministic label colours, readable legends, and 300-dpi PNG (or matching PDF when requested).

- Use 100% stacked bars only for showing complete sample composition; retain one bar per sample.
- Use jittered sample points for comparisons. Overlay mean and standard error only as a summary layer.
- Use a blue-white-red heatmap for proportions, with the scale stated in the subtitle.
- Pair every proportion view with a count view.
- Facet by condition or batch only when the grouping has at least two observed levels.
- Paginate high-cardinality category axes rather than selecting a hidden top-N subset.
- UMAP is a diagnostic view of representation and metadata distribution, not a batch-effect test by itself.
