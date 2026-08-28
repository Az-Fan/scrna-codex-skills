# scRNA-seq Codex Skills

Ten Codex/Wisp-compatible skills for auditable single-cell RNA-seq analysis. The canonical development working tree is `/home/faz_laptop/projects/scrna-codex-skills` on `ssh:xiyouyun`; GitHub is its synchronized distribution remote. Server-side analysis uses the existing pixi environments.

## Released skills

- `01-scrna-standardize-input`
- `02-scrna-calculate-qc-metrics`
- `03-scrna-review-qc`
- `04-scrna-preprocess-and-cluster`
- `05-scrna-benchmark-integration`
- `06-scrna-find-cluster-markers`
- `07-scrna-annotate-cells`
- `08-scrna-analyze-subset`
- `09-scrna-score-programs`
- `10-scrna-run-differential-analysis`

The numeric prefixes follow the primary execution flow. Integration benchmarking
(`05`) is conditional when batch correction needs comparison, and focused subset
reanalysis (`08`) is a downstream branch after a usable broad annotation exists.

## Install the fixed release

Clone the fixed tag from GitHub, then assemble all ten self-contained skill directories with the bundled installer:

```bash
git clone --branch v2.0.0 --depth 1 git@github.com:Az-Fan/scrna-codex-skills.git
python3 scrna-codex-skills/scripts/install_skills.py --target .wisp/skills
```

For a personal Codex installation, use `--target ~/.codex/skills`. On Windows PowerShell, use `py` when that is the configured launcher. Existing installations are never overwritten unless `--force` is supplied.

To update an existing clone and replace an installed Wisp copy with this fixed release:

```bash
git -C scrna-codex-skills fetch --tags
git -C scrna-codex-skills checkout v2.0.0
python3 scrna-codex-skills/scripts/install_skills.py --target .wisp/skills --force
```

Restart the agent session after installation or update so skill discovery is refreshed. Installation changes agent instructions and bundled executors only; it does not install R/Python packages or alter the server pixi environments.

## Source architecture

`toolkit/` is the only source of shared Python and R runtime code. `release/runtime-manifest.json` declares the minimal runtime required by each skill. The installer and packager inject those files when building a self-contained installation; generated copies are not committed under `skills/`.

Validate and smoke-test assembly:

```bash
python3 scripts/validate_runtime_manifest.py
python3 scripts/smoke_install.py
```

Run the real ten-skill deterministic fixture test in the registered server environment:

```bash
python3 tests/e2e/run_fixture_e2e.py
```

Build `.skill` archives when a release is ready:

```bash
python3 scripts/package_skills.py --output dist
```

## Reproducible fixture

Generate the small four-sample Seurat fixture inside a compatible server pixi environment:

```bash
Rscript tests/fixtures/create_fixture.R tests/fixtures/tiny_scrna.rds
```

The generated RDS is intentionally ignored by Git. The generator is deterministic and provides 80 cells, two conditions, two batches, two cell types, two clusters, integer counts, QC metadata, and a UMAP embedding for end-to-end tests.

The client needs Codex and/or WispScience. The registered analysis server must have pixi and the project environments under `~/projects/scrna_envs`. See each skill's `references/compatibility.md`.

Repository synchronization and release invariants are enforced by [AGENTS.md](AGENTS.md).
