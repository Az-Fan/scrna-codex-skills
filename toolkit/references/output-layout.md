# Result layout

The stage output directory is the browsing surface for scientific results. Keep primary figures, complete scientific tables, summaries, and human review decisions accessible there. Do not create empty category directories or duplicate every table and figure in a second format.

## Technical records

New runs write execution manifests, session information, workflow state, task manifests, and runner logs under `<output_dir>/_provenance/`. Session files use `session_info.txt`. The shared runner writes default validation plans under `<config-directory>/_provenance/`; an explicit `--manifest` path is respected. Place tmux supervisor logs/status there too using the launcher's explicit `--log` and `--status` arguments.

Custom follow-up plots follow the same rule: the figure and any useful scientific summary belong in the result directory, while the plot's manifest/session record belongs in `_provenance/`. Avoid treating all JSON or all TSV files alike: recommendation status, design audits, coverage, and missing/failed task reports can be essential to interpreting results. Surface material failures or limitations in the handoff even when their technical details are stored below `_provenance/`.

Read `_provenance/workflow_state.json` for guided clustering state. For existing runs created by older versions, fall back to `workflow_state.json` at the result root. The same new-path-first lookup applies to old manifests and session files (`session_info.txt` or `sessionInfo.txt`). The executors do not migrate old files; rerunning into an old output directory does not remove its existing records.

## Objects, supplements, and project rules

Keep optional expanded diagnostics in `details/` only when requested or needed to explain a result. Preserve complete scientific tables; top-N tables and plots are views, not replacements. Do not suppress design failures, untestable populations, missing metrics, or database failures to make the directory look clean.

Honor project-specific object and figure conventions when supported by the executor. If a project requires objects in `data/interim/` but the executor fixes object paths under `output_dir`, identify that limitation before execution; do not silently move the object afterward or break finalization and downstream references. Object-path and plot-format changes require their own contract verification. This layout update does not change scientific computations, object destinations, or figure formats.

Do not delete manifests, sessions, workflow state, or large objects as cosmetic cleanup. Existing-result migration is a separate operation: inventory exact paths and downstream references, validate destinations, preserve original manifests, and record moves and checksums. Do not infer that similarly sized objects are duplicates. Historical provenance remains immutable; follow the project's history-retention location when one is specified.

## Handoff

Link the main figures, complete result tables, necessary decisions, and downstream object. Mention `_provenance/` once as the troubleshooting location; do not enumerate its files as deliverables. Report validation and any material failures honestly.
