# Fixed figure catalog — baseline inventory and rollout

Rollout: skills 06/07/08 now use `paper_v1`, default PNG 300 dpi, with fixed pagination; their baseline PDF filenames below identify the retained families, not the new default extension. Skill 08 additionally supports explicitly configured target-gene FeaturePlots. All other families retain their existing renderer pending the next rollout. See figure-style.md in the three updated skills for their implemented contract.

This catalog records existing plot families before visual redesign. A family may produce multiple files or pages by sample, scenario, comparison, database, or metric. All listed families are retained. Conditions below describe current behavior, not permission to run additional analyses. No-figure stages are intentional.

| Skill | Figure family / current filename | Trigger and multiplicity | Disposition |
|---|---|---|---|
| 01 | No figure | Standardization writes objects and tables | Preserve |
| 02 | qc_diagnosis.png | One counts-versus-detected-genes scatter per Seurat run or STARsolo sample | Preserve |
| 03 | qc_distribution_by_sample.png | QC metrics available; one faceted family | Preserve in atlas |
| 03 | qc_distributions.png | QC metrics available; histograms by metric | Preserve in atlas |
| 03 | qc_metrics_per_cluster.png | Cluster metadata available | Preserve in atlas |
| 03 | qc_distribution_by_celltype.png | Annotation metadata available | Preserve in atlas |
| 03 | n_UMIs_vs_n_genes.png | Both count metrics available; sample facets | Preserve in atlas |
| 03 | nuclear_frac_vs_n_UMIs.png | Both metrics available; sample facets | Preserve in atlas |
| 03 | filter_mito_by_annotation.png | Mitochondrial metric and annotation available | Preserve in atlas |
| 03 | umap_annotation_qc_context.png | UMAP and annotation available | Preserve in atlas |
| 03 | umap_qc_<metric>.png | One panel per available metric when UMAP exists | Preserve in atlas |
| 03 | candidate_retention_by_sample.png | Candidate thresholds available; never formal filtering | Preserve in atlas |
| 03 | qc_atlas.pdf | Ordered container for all eligible families above; individual PNGs currently require full mode | Preserve every eligible panel |
| 04 | No figure | Approved filter produces retention/decision tables | Preserve; do not invent plots silently |
| 05 | umap_by_batch__<scenario>.png | Requested; one panel per available batch field and completed scenario | Preserve |
| 05 | umap_by_sample__<scenario>.png | Requested; sample field and completed scenario | Preserve |
| 05 | umap_by_condition__<scenario>.png | Requested; condition field and completed scenario | Preserve |
| 05 | umap_by_label__<scenario>.png | Requested; biological label fields and completed scenario | Preserve |
| 05 | metric_tradeoff.png | Requested; nonempty method summary | Preserve |
| 05 | score_heatmap.png | Requested; completed metric rows | Preserve |
| 05 | score_barplot.png | Requested; nonempty method summary | Preserve |
| 05 | ranking_plot.png | Requested; method summary; ranking interpretation still requires valid scientific evidence | Preserve, review interpretation |
| 05 | marker_dotplot__<scenario>.png | Requested gene programs with matching genes; supported non-graph scenario | Preserve |
| 05 | program_retention.png | Requested gene programs with computable retention | Preserve |
| 06 | <scenario>_umap_diagnostics.pdf | Finalized/fixed clustering; selected metadata panels; optional duplicate PNG | Preserve |
| 06 | <scenario>_resolution_stability.png | Resolution scan per scenario | Preserve |
| 06 | <scenario>_umap_clusters_by_resolution.png | Resolution scan; panels share the same UMAP | Preserve all resolutions |
| 06 | <scenario>_clustree_resolution.png | Resolution scan and clustree available; otherwise record failure | Preserve |
| 06 | <first_scenario>_elbow.png | PCA diagnostic currently for first scenario | Preserve; comparison expansion requires explicit change |
| 07 | top_marker_dotplot.pdf | Configured top markers across existing clusters | Preserve |
| 08 | cluster_umap.pdf | prepare_review; one cluster view | Preserve |
| 08 | cluster_sample_umap.pdf | prepare_review; cluster and sample panels | Preserve |
| 08 | canonical_marker_dotplot.pdf | prepare_review; configured canonical genes with matches | Preserve |
| 08 | annotated_umap.pdf | apply_confirmed; broad and fine annotation panels | Preserve |
| 08 | cluster_sample_condition_umap.pdf | apply_confirmed; available audit grouping fields | Preserve |
| 09 | No figure | Subset export produces object/matrix/tables | Preserve |
| 10 | <task>_group_mean_heatmap.png | Configured summary groups; at least two signatures and groups | Preserve |
| 11 | volcano.pdf | Each successful population-by-comparison task | Preserve |
| 11 | MA_plot.pdf | Result includes baseMean | Preserve |
| 11 | pseudobulk_PCA.pdf | Sample-level normalized pseudobulk data | Preserve |
| 11 | top_DE_heatmap.pdf | At least two selected genes in pseudobulk data | Preserve |
| 11 | DEG_count_summary.pdf | Significant DE genes across tasks | Preserve |
| 11/12 | enrichment_dotplot_overview.<format> | GSEA terms pass display selection; overview by database and direction | Preserve |
| 11/12 | enrichment_ora_overview.<format> | ORA terms pass display selection | Preserve |
| 11/12 | enrichment_dotplot_<database>_ora[_pageN].<format> | Eligible ORA terms; configured terms-per-page | Preserve every page |
| 11/12 | gsea_nes_<database>[_pageN].<format> | Eligible GSEA terms; configured terms-per-page | Preserve every page |
| 13 | sample_composition.pdf | Sample counts; currently 24 samples per page | Preserve every sample |
| 13 | cell_type_proportions_by_condition.pdf | Sample proportions; currently 12 cell types per page | Preserve every cell type |
| 13 | sample_proportion_heatmap.pdf | Sample proportions; currently 25 cell types per page | Preserve every cell type |
| 13 | effect_summary.pdf | Each nonempty method result; currently 30 features per page | Preserve every feature |
| 13 | milo_da_beeswarm.pdf | Milo method succeeds and plotting is available | Preserve; failure must be surfaced |
| 13 | milo_da_graph.pdf | Milo neighborhood graph can be built and plotted | Preserve; failure must be surfaced |

## Project supplements requiring a home

The PAH project contains additional plots outside the bundled executors. Inventory them separately; do not claim they are already implemented by these skills. They remain intact in their source project.

| Additional family | Proposed owner | Required distinction |
|---|---|---|
| Target/canonical-gene FeaturePlots | 08, configurable gene lists | Expression assay, shared scale and split condition need explicit rules |
| Pan-EC, EC subtype, common-cell marker dotplots | 07/08, named marker panels | Do not hardcode PAH-specific genes into generic skills |
| QC FeaturePlot grids | 03 | Reuse available QC metrics; distinguish grids from the existing atlas pages |
| Cluster/cell-type composition panels | 06 descriptive / 13 inferential | Show sample denominators; pooled counts are not biological replication |
| Target-gene boxplots by sample/condition | 11 | Cell distributions and sample-level expression are different estimands |
| Target-gene violins by sample/condition | 11 | Preserve descriptive cell-level labeling |
| Target-gene sample-expression heatmaps | 11 | Distinguish row z-score from absolute expression and cross-dataset effect heatmaps |
| Candidate QC filter review panels | 03/04 boundary | Preserve candidate-versus-approved distinction; no new filtering decisions |
| Target-gene effect/expression summary | 11 | Preserve contrast direction, test method and evidence level |
| Cross-dataset concordance and gene-table figures | Cross-dataset workflow | No existing single-dataset skill owns this analysis; retain as project outputs |

## Change and verification rules

Before changing a family, record the final style, file pattern, trigger, per-page capacity, and verification status. Keep a stable order and filenames; distinguish families from physical file count. Never silently omit a family because a renderer or optional method fails. Do not alter clustering, significance thresholds, gene selection, or model fitting to make figures look better.

`scripts/audit_figures.py` inventories source call sites, all existing image files, hashes, and PNG navigation sheets. Its scanner is a coverage aid, not a parser proving completeness: manually review plotting wrappers, dynamic names, base graphics, optional branches, and every PDF page. A listed file is not a visually reviewed file. Default E2E success does not cover optional methods or every high-cardinality layout.
