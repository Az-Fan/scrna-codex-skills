# scRNA-seq Codex Skills

Nine self-contained Codex skills for auditable single-cell RNA-seq analysis. The authoritative source is this repository; server-side analysis uses the existing pixi environments.

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

## Install

```bash
git clone --branch v1.1.0 --depth 1 git@github.com:Az-Fan/scrna-codex-skills.git
python scrna-codex-skills/scripts/install_skills.py
```

On Windows PowerShell, use `py` when that is the configured launcher. Pass `--target PATH` for another Codex-compatible directory. Existing installations are never overwritten unless `--force` is supplied.

The client needs Codex and/or WispScience. The registered analysis server must have pixi and the project environments under `~/projects/scrna_envs`. See each skill's `references/compatibility.md`.

## Verify

```bash
python scripts/smoke_install.py
```
