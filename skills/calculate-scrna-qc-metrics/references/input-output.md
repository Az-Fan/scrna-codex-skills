# Input and output contract

## Input

For each declared `sample_id`, require:

```text
<starsolo_dir>/<sample_id>/Solo.out/Gene/filtered/
  matrix.mtx.gz
  features.tsv.gz
  barcodes.tsv.gz
<starsolo_dir>/<sample_id>/Solo.out/Velocyto/filtered/
  spliced.mtx.gz
  unspliced.mtx.gz
  barcodes.tsv.gz
```

Require a matching gene-annotation GTF. The current executor uses mouse gene symbols, `^mt-`, mouse cell-cycle symbol conversion, and a mouse hemoglobin list.

Require an existing pixi executable and the existing QC pixi project. The skill never creates or changes the environment.

## Per-cell metadata

- `sample_id`, `batch_id`
- `n_genes`, `n_UMIs`
- `mito_frac`, `chrY_frac`, `nuclear_frac`
- `ambient_frac_decontx`, `doublet_score`
- `phase`, `s_score`, `g2m_score`
- `hbb_score`, `is_HQ`

## Output

For every sample, write under `<output_dir>/<sample_id>/`:

- `counts.mtx.gz`
- `features.tsv.gz`
- `metadata.tsv.gz`
- `qc_diagnosis.png`

Write `run_manifest.json` at the output root. No cells are removed.
