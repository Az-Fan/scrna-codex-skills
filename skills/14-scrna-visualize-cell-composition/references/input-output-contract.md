# Input and output contract

## Input

Use exactly one of `input.object` (Seurat RDS/QS) or `input.counts_table` (TSV). An object requires `metadata.sample`; each configured composition variable must be a categorical metadata column. A counts table requires `sample`, every configured category column, and `count_column`.

`composition.denominator.mode=all_input_cells` uses all retained cells. `selected_parent` or `parent_column` recomputes proportions within each parent label. A denominator description is required and is copied into the output audit.

## Outputs

Root outputs include `composition_counts.tsv`, `composition_proportions.tsv`, `sample_coverage.tsv`, `composition_audit.tsv`, and `plot_status.tsv`. Each configured variable gets `composition_overview_<name>`, `composition_dotplot_<name>`, `composition_heatmap_<name>`, and `composition_counts_<name>` in the requested format. UMAP diagnostics are written when a requested reduction and Seurat object are available.

Technical records are under `_provenance/`, including the resolved fields, session information, figure colours, and run manifest.

## Interpretation

Each sample-level proportion is category cells divided by the declared sample denominator. Condition bars are summaries of sample proportions, not pooled-cell estimates. A changed relative proportion does not identify proliferation, death, migration, recruitment, or an absolute tissue change.
# Parent-specific figures

When `composition.parent_column` is configured, every parent has its own overview, dotplot, heatmap, and count figure. Filenames include `_parent_<index>_<label>`; within-parent proportions are never stacked or averaged across parents. `selected_parent` requires the parent column and explicit parent labels in `denominator.include`. `selected_cell_types` selects the configured variable's categories. Sample IDs must map to one condition and batch.
