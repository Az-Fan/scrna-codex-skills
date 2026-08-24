# Input and output contract

## Input modes

### STARsolo

For each declared sample, require:

```text
<starsolo_dir>/<sample_id>/Solo.out/Gene/filtered/
  matrix.mtx.gz
  features.tsv.gz
  barcodes.tsv.gz
```

Use `Velocyto/filtered/{spliced,unspliced}.mtx.gz` when present; otherwise skip `nuclear_frac`. Use a matching GTF when provided; otherwise skip `chrY_frac`.

### Seurat

Require an RDS or QS Seurat object with raw counts. A sample metadata column is optional; without one, treat the object as one sample. A batch column is optional. Preserve the input and write a derivative object.

Set `parameters.species` to `mouse` or `human` for mitochondrial, cell-cycle, and hemoglobin symbols.

Require an existing pixi executable and pixi project. Never create or change the environment.

## Metrics

- `sample_id`, `batch_id`
- `n_genes`, `n_UMIs`
- `mito_frac`, `chrY_frac`, `nuclear_frac`
- `ambient_frac_decontx`, `doublet_score`
- `phase`, `s_score`, `g2m_score`
- `hbb_score`, `is_HQ`

Write `metric_status.tsv` with `computed` or `skipped` and a reason for every optional metric.

## Output

For STARsolo, write under `<output_dir>/<sample_id>/`:

- `counts.mtx.gz`
- `features.tsv.gz`
- `metadata.tsv.gz`
- `qc_diagnosis.png`
- `metric_status.tsv`

For Seurat, write under `<output_dir>/`:

- `qc_metrics_object.rds`
- `metadata.tsv.gz`
- `qc_diagnosis.png`
- `metric_status.tsv`

Write `run_manifest.json` at the output root. Never remove cells.
