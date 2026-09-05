---
name: 01-scrna-standardize-input
description: Standardize scRNA-seq inputs and metadata from 10x matrix directories or H5 files and Seurat RDS/QS objects. Use when starting or migrating a public or local single-cell RNA-seq project, diagnosing incompatible inputs, or preparing a stable object contract for downstream skills.
---

# Standardize scRNA Input

Create an auditable input layer without changing biological values unnecessarily.

## Workflow

1. Inspect files read-only and identify the matrix orientation, format, assay, species, gene identifiers, barcodes, and available metadata.
2. Create or validate `samples.tsv`. Require unique sample identifiers and explicit condition labels; retain batch when available.
3. Map source columns to the canonical roles `sample`, `condition`, `batch`, and `cell_type`. Do not silently rename ambiguous fields.
4. Preserve raw integer counts. Keep normalized values in a separate assay or layer.
5. Make cell identifiers globally unique while retaining original barcodes and provenance.
6. Validate matrix/metadata alignment, duplicated genes, missing samples, and condition-batch confounding.
7. Write the standardized object, `samples.tsv`, field mapping, and provenance report to a new output directory.

## Guardrails

- Prefer a processed count matrix for exploration when it is documented and adequate; recommend FASTQ requantification only for a concrete reason.
- Never infer experimental groups from filenames without reporting the inference.
- Never overwrite raw inputs.
- Stop when species, sample identity, or matrix orientation cannot be established safely.

Read [references/input-contract.md](references/input-contract.md) for the canonical contract. Use the repository tool `toolkit/python/validate_project.py` to validate sample tables.

## Execution

Create a project-local JSON config from [references/config.example.json](references/config.example.json). Run `scripts/run.py --config <config>` in the selected compute context before changing data. Treat a nonzero exit as blocking. Rerun with `--execute` only after the dry-run manifest is correct. Use the bundled default R executor; set a safe argv array under `executor.argv` only to override it. Never use a shell command string.

## Result organization

Keep primary figures, complete scientific tables, and review decisions directly accessible. Store execution manifests, session information, logs, and workflow state under `_provenance/`; do not list them as primary results. Read [references/output-layout.md](references/output-layout.md) when configuring outputs, locating legacy records, or adding custom plots and diagnostics.
