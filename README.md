# scRNA-seq Codex Skills

Nine Codex/Wisp-compatible skills for auditable single-cell RNA-seq analysis. The authoritative source is this repository; server-side analysis uses the existing pixi environments.

## Released skills

- `scrna-calculate-qc-metrics`
- `scrna-review-qc`
- `scrna-standardize-input`
- `scrna-analyze-subset`
- `scrna-annotate-cells`
- `scrna-find-cluster-markers`
- `scrna-benchmark-integration`
- `scrna-score-programs`
- `scrna-run-differential-analysis`

## Install the fixed release

```bash
git clone --branch v1.1.0 --depth 1 git@github.com:Az-Fan/scrna-codex-skills.git
python scrna-codex-skills/scripts/install_skills.py
```

On Windows PowerShell, use `py` when that is the configured launcher. Pass `--target .wisp/skills` for a project-local Wisp installation. Existing installations are never overwritten unless `--force` is supplied.

## Source architecture

`toolkit/` is the only source of shared Python and R runtime code. `release/runtime-manifest.json` declares the minimal runtime required by each skill. The installer and packager inject those files when building a self-contained installation; generated copies are not committed under `skills/`.

Validate and smoke-test assembly:

```bash
python scripts/validate_runtime_manifest.py
python scripts/smoke_install.py
```

Build `.skill` archives when a release is ready:

```bash
python scripts/package_skills.py --output dist
```

## Reproducible fixture

Generate the small four-sample Seurat fixture inside a compatible server pixi environment:

```bash
Rscript tests/fixtures/create_fixture.R tests/fixtures/tiny_scrna.rds
```

The generated RDS is intentionally ignored by Git. The generator is deterministic and provides 80 cells, two conditions, two batches, two cell types, two clusters, integer counts, QC metadata, and a UMAP embedding for future end-to-end tests.

The client needs Codex and/or WispScience. The registered analysis server must have pixi and the project environments under `~/projects/scrna_envs`. See each skill's `references/compatibility.md`.
