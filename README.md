# scRNA-seq Codex Skills

Seven self-contained Codex skills for auditable single-cell RNA-seq analysis. The authoritative source is this repository; server-side analysis uses the existing locked pixi environments.

## Released skills

- `scrna-standardize-input`
- `scrna-analyze-subset`
- `scrna-annotate-cells`
- `scrna-find-cluster-markers`
- `scrna-benchmark-integration`
- `scrna-score-programs`
- `scrna-run-differential-analysis`

## Install

Clone a fixed release and install all skills:

```bash
git clone --branch v1.0.0 --depth 1 git@github.com:Az-Fan/scrna-codex-skills.git
python scrna-codex-skills/scripts/install_skills.py
```

On Windows PowerShell, use `py` instead of `python` when that is the configured launcher. To install into another Codex-compatible directory, pass `--target PATH`. Existing installations are never overwritten unless `--force` is supplied.

The client only needs Codex and/or WispScience. The registered analysis server must have pixi and the project environments under `~/projects/scrna_envs`; see each skill's `references/compatibility.md` and run its bundled dependency checker there.

## Verify

```bash
python scripts/smoke_install.py
```

The smoke suite copies every skill into an isolated temporary directory, compiles Python sources, and starts both the executor and dependency-check CLI without relying on the repository-level toolkit.
