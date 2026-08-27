# Input contract

Required deliverables:

- Raw counts with genes as rows and globally unique cells as columns.
- Cell metadata containing a stable sample identifier.
- `samples.tsv` containing one row per biological sample and at least `sample_id` and `condition`.
- Provenance recording source paths, input format, species, gene identifier type, and transformations.
- A standardized Seurat object. With `output.object_format: "auto"`, the executor writes `standardized_object.qs` when `qs` is available and otherwise writes `standardized_object.rds`. Set the value explicitly to `qs` or `rds` only when a fixed serialization format is required.

Recommended canonical roles are `sample_id`, `condition`, `batch`, `cell_type`, and `cluster`. Store the original source column names in a mapping rather than forcing all project objects to use identical physical names.
