# Core GRN contract

## Scope and inputs

This skill distills sc06 `04-grn` steps 01–02 into configurable input aggregation, GRNBoost2, motif pruning and AUCell. Step 01's PCA/resolution scan belongs to preprocessing skill 06; an existing reviewed membership column is required for metacells. sc06 steps 03–06 and their interpretation are deferred. No full derivative Seurat object is saved.

Seurat: `input.type=seurat`, object RDS/QS, assay RNA, counts_source raw_umi. Split Seurat v5 count layers are joined in memory. Matrix: sparse coordinate Matrix Market counts (.mtx or .mtx.gz, as exported by skill 09), headerless features/barcodes TSV, metadata TSV first column cell ID; orientation genes_by_cells or cells_by_genes; feature_column is one-based. Dense Matrix Market array format is unsupported by the R reader. Supply RNA-only unique identifiers. All cells must match metadata exactly. `metadata.sample` is required; configured condition and batch must be constant per sample.

Declare input species/genome/gene_identifier and matching resource declarations. This is a provenance declaration, not proof from filenames. TF list: one ID per line, no duplicates. Ranking databases: supported genes_vs_motifs.rankings.feather files. Motif annotations must match the database motif collection and gene identifier system. No automatic ortholog conversion, database substitution or downloads. Coverage files show TF and database matches; small overlap is a review concern even when nonzero. At least two variable TFs must match. rank_threshold must be below each database's gene count.

## Inference units and stages

`inference.mode` has no inferred default. single_cell uses original cells; metacell sums raw counts by (sample, configured cell_type, membership_column). Aggregation never crosses those strata. Original labels remain in cell_membership.tsv. Check n_cells and normalized-entropy purity for configured metadata fields. `metacell.min_cells` defaults to 1; a smaller unit stops the run without dropping cells. At least three units and three variable genes are required. Constant/zero genes are omitted only from inference, recorded in feature_status.tsv, and retained in the single-cell AUCell background.

`prepare` exports and audits, stopping at awaiting_input_confirmation. `infer` requires workflow.review_reason and reuses a SHA256-bound prepared input with unchanged input/resource configuration and Python versions. `run` performs both for already selected inputs/parameters. Config declarations are compared exactly; preserve them across stages. Source inputs/resources must be outside output_dir. No implicit resumption inside failed inference: preserve partial files and use a new output directory. TF/motif databases are external, never bundled in skills.

GRNBoost2 uses the multiprocessing Arboreto implementation, a fixed seed (default 777), and grn.workers (default 1). Context pruning uses custom_multiprocessing, positive modules by upstream default, ctx.min_genes=20, rank_threshold=1500, auc_threshold=0.05 and nes_threshold=3.0. ctx.min_genes filters pre-pruning modules; regulons.min_genes=20 filters motif-pruned gene sets and applies consistently to exports/scoring. Upstream regulon IDs such as TF(+) and contexts remain unchanged. Edge importance/target weights are not p-values or evidence of binding. One seeded inference is not a stability benchmark.

## Activity and outputs

AUCell uses raw counts for all original cells/genes; original IDs and metadata are retained. aucell.rank_fraction=0.05 defines floor(n_genes × fraction) ranked genes and must yield at least two. aucell.batch_size=500 controls memory; each batch is scored serially with seed + batch_index. Zero-expression ties are shuffled by AUCell; batch size/order/seed affect tie ranks, so record them. Normalized AUC values are rankings/enrichment scores, not fractional expression-program usage. Sparse counts are retained, while rankings and inference workers may allocate dense matrices.

| Deliverable | Contract |
|---|---|
| input_audit.json, feature_status.tsv | Source/resource hashes, dimensions and inference gene exclusions |
| cell_membership.tsv, inference_unit_audit.tsv | Every cell's inference unit; size and metadata purity |
| tf_coverage.tsv, database_coverage.tsv | Variable TF and ranking-database gene overlap |
| adjacencies.tsv | TF, target, GRNBoost2 importance; coexpression candidates |
| motif_enrichment.tsv | Complete native ctx motif enrichment table |
| motif_row_status.tsv | Per-row target size; empty leading edges excluded from regulon aggregation |
| regulon_summary.tsv | Native regulon IDs, TF, size, context and retained/excluded status |
| regulon_targets.tsv, regulons.gmt, regulon_list.rds | Retained TF-target weights and equivalent gene sets |
| regulon_activity.tsv, regulon_activity.rds | Original cells × retained regulons, AUCell AUC |
| regulon_activity_by_sample.tsv | Sample × optional cell type/condition/batch descriptive means and n_cells |
| task_status.tsv, handoff.json | Stage outcomes and core handoff; downstream not run |

Review starts at RESULTS.md. `_provenance/grn_input/` contains matrix/ID inputs and inference loom. Other provenance includes commands, versions, prepared binding, R/AUCell session information and run_manifest.json. No automatic regulon binarization, Seurat-assay attachment or advanced downstream plots.

Runtime profile `04-grn`: default R for input/AUCell, pyscenic Python for inference. Root precedence runtime.pixi_root, SCRNA_PIXI_ROOT, ~/projects/scrna_envs. Explicit runtime.pyscenic_python overrides without fallback. Analysis does not install/repair dependencies. Use repository environment deployment separately.

Sources: [pySCENIC tutorial](https://pyscenic.readthedocs.io/en/latest/tutorial.html), [multiprocessing GRN implementation](https://github.com/aertslab/pySCENIC/blob/master/src/pyscenic/cli/arboreto_with_multiprocessing.py). Follow installed CLI behavior when package versions change.
