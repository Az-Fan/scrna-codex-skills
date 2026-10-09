# cNMF discovery contract

## Inputs and runtimes

Seurat input declares `input.type=seurat`, `input.object`, `input.assay` (default RNA), and `input.counts_source=raw_umi`. Integer validation alone cannot prove counts are uncorrected. RDS/QS are supported; Seurat v5 counts layers are joined in memory. No full derivative object is saved.

For matrix input use:

```json
{"type": "matrix", "counts_source": "raw_umi", "counts": "data/subset_counts.mtx", "features": "data/features.tsv", "barcodes": "data/barcodes.tsv", "metadata": "data/subset_metadata.tsv", "orientation": "genes_by_cells", "feature_column": 1}
```

Matrix Market/gzip is supported. Orientation (`genes_by_cells` or `cells_by_genes`) is declared, not inferred. Features/barcodes are headerless TSV; feature_column is one-based (2 for 10x symbols). Metadata TSV has cell IDs in its first column, compatible with skill 09 exports. IDs must match exactly without extras or duplicates; gene IDs must be unique. No automatic symbol conversion or feature-type filtering occurs: provide RNA-only features.

`metadata.sample` is required. Optional condition, cell_type and batch columns must contain non-empty values when configured; condition/batch must be constant within a sample. Optional `metadata.reduction` plots an existing Seurat embedding without fitting new coordinates. Omit it for matrix inputs, which produce sample heatmaps without R.

Counts must be finite non-negative integers; zero-count cells stop execution. Zero-expression genes are audited and omitted. Choose `num_highvar_genes` no larger than expressed-gene count and each k smaller than both HVG and cell counts. HVG selection and variance normalization are performed by cNMF.

R export uses registered `05-pathway_program/.pixi/envs/default/bin/Rscript`; cNMF uses `.pixi/envs/cnmf/bin/python`. Root precedence: runtime.pixi_root, SCRNA_PIXI_ROOT, ~/projects/scrna_envs. An explicit `runtime.cnmf_python` overrides Python without fallback. Dependency checks use the same resolution. One BLAS thread per process; cnmf.workers controls process parallelism. Native regression/refitting may densify expression, so choose population and workers for available RAM.

## Rank and density review

- `discover` scans candidate ranks and stops for user review.
- `consensus` reuses the same discovery, requiring explicit consensus_k and workflow.selection_reason. Keep input roles and preparation parameters unchanged. Only new ranks from the original scan can be exported.
- `run` scans and exports user-specified consensus_k in one execution.

`n_iter` counts independent seeded repetitions, whereas `max_nmf_iter` limits solver iterations. Default 100 repetitions follows upstream examples. Fixture runs with fewer repeats are functional tests, not scientific stability evidence. Rank diagnostics show unfiltered replicate silhouette and reconstruction error; neither establishes a biological optimum or density-filtered consensus quality.

Smaller density_threshold means stricter replicate-spectra filtering. Default 0.02 follows sc06 and requires clustergram review; it is not universal. local_neighborhood_size sets neighbors relative to repetitions. Too-strict filtering may fail; never silently loosen it. Use a fresh output run for parameter experiments and preserve previous consensus.

## Outputs

| File in k_<k> | Meaning |
|---|---|
| usage_raw.tsv | Cells × programs, upstream refitted usage |
| usage.tsv | Cells × programs, row sum 1 |
| spectra_scores.tsv | Genes × programs, signed regression scores |
| spectra_tpm.tsv | Genes × programs, TPM-scale weights |
| top_genes.tsv | Ranked genes, scores and TPM weights per program |
| programs.gmt | Up to top_n positive-score genes per program |
| program_id_mapping.tsv | Upstream IDs to k-namespaced exported IDs |
| usage_summary_by_sample.tsv | Sample × optional population/condition/batch, mean fractional usage and cell counts |
| program_review.tsv | Pending candidate names, evidence, sample consistency and QC confounding |

cNMF 1.7 load_results actually returns spectra as genes × programs despite contradictory docstring wording. Validate axes against usage programs. Program IDs have no fixed biological correspondence across k or runs. Cross-k Spearman correlations report all pairs; constant programs can yield undefined values, and correlations do not certify one-to-one matches.

Heatmaps show descriptive per-sample/population means. Cells are not independent replicates; condition differences can reflect composition, technical effects, cell identity, stress or cell cycle. No automatic biological names, condition p-values, mixed-model variance decomposition or online database enrichment are produced.

Optional interpretation.gmt supplies a reviewed local species/identifier-matched GMT. ORA uses positive-score top genes and all scored genes as the eligible universe; all matched sets are tested for all programs. BH correction covers all program × gene-set tests within k. Retain complete results, coverage and universe files. This annotates genes functionally; it does not test condition differences.

## Preservation and failures

Discovery refuses existing prepared/native artifacts; consensus refuses existing selected-k directories. Failed ranks retain task_status.tsv errors and partial files, successful requested ranks remain accessible, and execution exits nonzero. Use a fresh run to repair failed output; no implicit retries or replacement of prior consensus occur.

_provenance/discovery_record.json binds input SHA256, preparation parameters, metadata roles, package versions and reusable artifact hashes. Consensus verifies these before reuse. This establishes integrity, not human authorization. Retain native factorizations and prepared input for reuse. RESULTS.md reports awaiting review, completion or failure and links scientific deliverables.

API sources: [upstream cNMF workflow](https://github.com/dylkot/cNMF), [stepwise guide](https://github.com/dylkot/cNMF/blob/main/Stepwise_Guide.md). Match installed behavior when upstream changes.
