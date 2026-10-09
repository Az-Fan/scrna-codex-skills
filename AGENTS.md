# Repository operating rules

## Canonical source and synchronization

- Treat this Git checkout as the source of truth for the task. Confirm its remote and branch before publishing; repository use does not require a particular host, username, or absolute path.
- Keep machine-specific checkout, runtime and host settings in ignored `local-development.json` (see `local-development.example.json`). These settings document the local workspace and do not change skill runtime contracts or select publishing destinations automatically.
- Maintain the seven skill environment profiles, locks and supplemental sources under `environments/`; deploy them to the shared runtime directory with `scripts/manage_environments.py`. Do not independently edit installed public environment configuration as source. After changes, update bundle checksums and verify environment audit plus relevant selected-method checks.
- Do not edit installed skills or build directories as independent source copies. Reconcile deliberate development branches through Git before distributing their outputs.
- Make source changes in the verified working tree. Build installed skills and `.skill` packages from that working tree through the repository scripts; never edit generated or installed copies as source.
- Before changing files, check the working tree, branch, remote, and upstream divergence. Preserve unrelated changes.
- After changing files, run the relevant checks and end-to-end tests, commit the verified source, push the commit to GitHub, then verify the remote branch resolves to the same commit.
- Create a release tag only from a clean, verified commit. Push the tag and verify the GitHub tag resolves to the exact same commit. Never move or overwrite an existing release tag.
- Install and update other computers from a fixed GitHub release tag when reproducibility matters. Use the default branch only when intentionally testing unreleased development.
- If the checkout and intended distribution remote differ unexpectedly, stop publishing or installing, diagnose the divergence, and reconcile it before continuing.

## Release gate

- Keep every released skill in `release/runtime-manifest.json` and build its self-contained runtime files with `scripts/install_skills.py` or `scripts/package_skills.py`.
- Run manifest validation, install smoke tests, package checks, and `tests/e2e/run_fixture_e2e.py` before a release.
- Do not claim install readiness unless every released skill passes and the source fixture remains unchanged.
