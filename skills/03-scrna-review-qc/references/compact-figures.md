# Compact scientific QC figures

Use `display.sample_labels` to map sample IDs to unique short labels and `display.sample_colors` to assign reproducible colors by display label. Do not merge samples sharing a time point. Counts use log10 axes; metric tables retain original units. Set `display.mito_unit` to `fraction` or `percent`; the default infers percent only for the source column `percent.mt` and fraction for `mito_frac`. Mitochondrial figures show percent.

`approved_rules` has optional lower/upper bounds in original metric units. Populate it only from an already approved decision, recording `approved_rules_source`. Red dashed overlays describe those decisions, never new candidate thresholds. The default is no approved overlay. `threshold_review.tsv` remains a distinct candidate table with blank approval fields; candidate retention is hypothetical and never shown as actual filtering.

Page 1: sample sizes, median summary, four core distributions (genes, UMI, mitochondrial fraction, doublet score when available), and shared-axis UMI/gene binned density. Boxes show median and interquartile range. There is no downsampling.

Page 2: numeric cluster order, cell counts, median heatmap and sample fractions. Heatmap colors standardize each metric across clusters, applying log10 to count medians first; tile text is the raw displayed median. Do not compare colors across metrics as absolute severity. `display.focus_clusters` optionally selects clusters for the lower row. Low complexity and sample imbalance alone do not justify exclusion.

Page 3: cluster, sample and core metric UMAP panels with identical coordinate limits and seeded point order. Set `display.umap_reduction` explicitly or use auto-detection of umap/standard_umap. Counts use log10 color values; other metrics use displayed units. Legends span all finite values without percentile clipping. Dense point layers are 1100-pixel PNG rasters; PDF axes, titles, legends and labels remain vector. Doublet scores are diagnostic, not probabilities.

Supplemental diagnostics use raw source units in a separate supplementary PDF. No plots infer cell types, determine formal batch effects or change filtering rules. Never save or alter the source Seurat object.

Dependencies: patchwork, png, viridisLite and scales in the existing selected environment. Never install them automatically.

The atlas has three sections, normally one page each. Sample overview pages contain at most eight samples; cluster review pages at most 24 clusters. Every group remains in tables and plots. Preview names add `_pageN` when a section spans multiple pages. Missing cluster/UMAP data produce an availability panel. Coordinates also fall back to UMAP_1/UMAP_2 metadata when no reduction is configured. Non-finite metrics are unavailable values; zero counts cannot appear on log axes. Required sample IDs cannot be missing.

`sample_display_summary.tsv` has sample IDs, display labels, cell counts and medians in display units. `qc_summary_by_sample.tsv`, `qc_summary_by_cluster.tsv` and thresholds retain stored metric units. `cluster_sample_composition.tsv` includes sample IDs and display labels; fractions sum to one within each cluster. Missing cluster/annotation labels are retained as a distinct missing category. `_provenance/figure_methods.json` records units, overlays, sample palette, pagination and diagnostic availability.

The supplemental PDF preserves all eligible existing diagnostic families: distributions, annotation/cluster violins, nuclear-versus-UMI, annotation mitochondrial candidate context, annotation UMAP, every available metric UMAP and hypothetical candidate retention. Supplemental metric axes use stored units; candidate displays remain exploratory. Full mode exports those same plots individually in details/. The plotting helper is bundled as scripts/review_qc_diagnostics.R.
