# Result layout

For figure families and rollout status, read [figure-catalog.md](figure-catalog.md). Its baseline list includes conditional and optional plots; only explicitly marked skills currently use the new paper style.

The stage output directory is the browsing surface for scientific results. Keep primary figures, complete scientific tables, summaries, and human review decisions accessible there. Do not create empty category directories or duplicate every table and figure in a second format.

## Result entry point

Executing `scripts/run.py --execute` writes `RESULTS.md` in every stage output root, with relative links grouped as main figures, complete tables and summaries, review decisions, and downstream objects/matrices. Comparison directories receive their own `RESULTS.md`, linked from the root. Copying the output directory preserves these links. Dry runs only validate and plan; they do not create a result entry point.

The entry point reports completed, partial, failed, or awaiting-confirmation state from the actual executor outcome and current status audits. Optional metric/plot skips remain visible as limitations without turning successful execution into failure. Empty enrichment, missing genes, low replication, and unresolved recommendations must remain visible. A successful process does not imply that every scientific task succeeded.

Each execution records the files written or changed during that run in `_provenance/result_delivery.json`. Unchanged files retained on disk from earlier runs are excluded from current deliverables and counted in the entry point. The guided clustering finalization and annotation apply stages explicitly inherit the preceding stage's registered review material, label those links as retained, and update the root entry point. Legacy two-stage outputs without a delivery registry identify pre-existing files as preceding-stage material. Scientific paths and files are preserved; the index does not certify historical results.

## Technical records

New runs write execution manifests, session information, workflow state, task manifests, and runner logs under `<output_dir>/_provenance/`. Session files use `session_info.txt`. Integration exchange matrices live in `_provenance/exchange/`; scoring hash caches live in `_provenance/score_cache/` while complete score TSV matrices remain in `scores/`. Default resource downloads use `_provenance/resource_cache/`, created only when a resource requires caching; explicit external cache directories are respected. Figures directories are created when figures are saved. The shared runner writes default validation plans under `<config-directory>/_provenance/`; an explicit `--manifest` path is respected. Place tmux supervisor logs/status there too using the launcher's explicit `--log` and `--status` arguments.

Custom follow-up plots follow the same rule: the figure and any useful scientific summary belong in the result directory, while the plot's manifest/session record belongs in `_provenance/`. Avoid treating all JSON or all TSV files alike: recommendation status, design audits, coverage, and missing/failed task reports can be essential to interpreting results. Surface material failures or limitations in the handoff even when their technical details are stored below `_provenance/`.

Read `_provenance/workflow_state.json` for guided clustering state. For existing runs created by older versions, fall back to `workflow_state.json` at the result root. The same new-path-first lookup applies to old manifests and session files (`session_info.txt` or `sessionInfo.txt`).

The shared execution runner records a unique run ID and the actual child exit status, including failures before an R manifest is written. For skills 11–13, execution into a nonempty output directory first preserves the whole prior directory under the sibling `.<output-name>-previous-runs/<run-id>/`, then creates a fresh current output directory. The current manifest records `previous_output`. Keep inputs and configs outside this output directory; overlapping input/config paths are rejected before any archive move. Use the documented `scripts/run.py --execute` entry point to get this supervision and rerun isolation. Other skills retain their existing cache/finalization/overwrite rules.

## Objects, supplements, and project rules

Keep optional expanded diagnostics in `details/` only when requested or needed to explain a result. Preserve complete scientific tables; top-N tables and plots are views, not replacements. Do not suppress design failures, untestable populations, missing metrics, or database failures to make the directory look clean.

Honor project-specific object and figure conventions when supported by the executor. If a project requires objects in `data/interim/` but the executor fixes object paths under `output_dir`, identify that limitation before execution; do not silently move the object afterward or break finalization and downstream references. Object-path and plot-format changes require their own contract verification. This layout update does not change scientific computations, object destinations, or figure formats.

Do not delete manifests, sessions, workflow state, or large objects as cosmetic cleanup. Existing-result migration is a separate operation: inventory exact paths and downstream references, validate destinations, preserve original manifests, and record moves and checksums. Do not infer that similarly sized objects are duplicates. Historical provenance remains immutable; follow the project's history-retention location when one is specified.

## Handoff

Lead the handoff with `RESULTS.md`. Link the main figures, complete result tables, necessary decisions, and downstream object. Mention `_provenance/` once as the troubleshooting location; do not enumerate its files as deliverables. Report validation and any material failures honestly.
