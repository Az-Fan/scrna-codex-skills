# Single-page marker overview

One primary figure: `top_marker_dotplot.pdf` by default, or a single PNG instead when explicitly selected. Never paginate this overview. Rows show marker genes in source-cluster blocks; columns show all clusters in numeric order when cluster IDs are numeric. Each block prints its source-cluster ID and cell count. Duplicate genes across source blocks are retained and reuse the same expression and detection values. Missing-marker clusters get an explicit no-marker block.

`reporting.dotplot_top_n` defaults to 4, independently of `reporting.top_n=20`. Select adjusted P<`reporting.adjusted_p_threshold` (default 0.05), sorting by adjusted P, descending effect and gene symbol. `marker_overview_selection.tsv` records every selected row, source cluster and original statistics. Do not shorten the full or top20 marker tables to match the figure.

Color: gene-wise scaled mean RNA expression across all clusters, using the Seurat DotPlot calculation once on unique selected genes. Dot size: detected fraction, limits 0–100%. Plot annotations do not imply cell-type or subtype identities. No sample-level inference is performed.

The default page size adapts to marker rows and cluster count; override `plots.overview_width` and `plots.overview_height` in inches if necessary. Larger figures remain one page; never silently paginate or drop clusters to fit a target paper size.

Use `workflow.action=plot_existing`, `input.object`, `input.markers`, `metadata.cluster` and `analysis.assay` for drawing from a previous complete marker table. This branch does not rerun FindAllMarkers, normalize or rewrite marker tables or Seurat objects. Normal mode continues to produce full markers, ranked top markers and cluster summaries. Both modes require the existing project environment.

Illustrative EC layout: 13 clusters, top4 each, 52 displayed marker rows; 8.8 × 13.4 inches. Numeric source-cluster blocks make shared or conflicting marker expression visible without splitting groups across pages.
