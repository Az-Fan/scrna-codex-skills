# Repository operating rules

## Canonical source and synchronization

- Treat `/home/faz_laptop/projects/scrna-codex-skills` on `ssh:xiyouyun` as the only development working tree and source of truth.
- Treat `git@github.com:Az-Fan/scrna-codex-skills.git` as the synchronized distribution remote, not as a second independently edited copy.
- Do not develop or preserve divergent copies in the Windows workspace, Codex/WispScience installed skill directories, build directories, or another computer.
- Make source changes in the canonical server working tree. Build installed skills and `.skill` packages from that working tree through the repository scripts; never edit generated or installed copies as source.
- Before changing files, check the canonical working tree, branch, remote, and upstream divergence. Preserve unrelated changes.
- After changing files, run the relevant checks and end-to-end tests, commit the verified source, push the commit to GitHub, then verify the remote branch resolves to the same commit.
- Create a release tag only from a clean, verified commit. Push the tag and verify the GitHub tag resolves to the exact same commit. Never move or overwrite an existing release tag.
- Install and update other computers from a fixed GitHub release tag when reproducibility matters. Use the default branch only when intentionally testing unreleased development.
- If the server commit and GitHub commit differ unexpectedly, stop publishing or installing, diagnose the divergence, and reconcile it in the canonical server working tree before continuing.

## Release gate

- Keep every released skill in `release/runtime-manifest.json` and build its self-contained runtime files with `scripts/install_skills.py` or `scripts/package_skills.py`.
- Run manifest validation, install smoke tests, package checks, and `tests/e2e/run_fixture_e2e.py` before a release.
- Do not claim install readiness unless every released skill passes and the source fixture remains unchanged.
