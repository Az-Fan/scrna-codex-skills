# Approved-filter contract

The decision table is the scientific source of truth. It must contain one row per input cell, a unique cell identifier, and every configured inclusion or exclusion field.

The executor calculates:

`keep = all(include_all_true) AND NOT any(exclude_any_true)`

Every configured decision field must contain only unambiguous boolean values (`TRUE/FALSE`, `1/0`, or `yes/no`). The executor rejects missing or unknown values. Set `exclude_any_true` to an empty array when no exclusion flag is needed.

`expected_retained_cells` is mandatory and must match the approved review. The filter does not estimate or adjust it.

The output object retains raw counts and input metadata for kept cells, plus configured decision-table columns. It intentionally drops normalized layers, reductions, graphs, neighbors, and exploratory cluster state so downstream preprocessing starts from the approved raw-count handoff.

Required outputs are the filtered object, complete compressed cell-decision table, sample and optional condition retention tables, decision-reason counts, approval record, session information, and run manifest.

## Output format

`output.object_format=auto` chooses QS when `qs` is available and otherwise RDS. Set `qs` or `rds` for a fixed format, or provide `output.object_name` with the explicit extension. An explicit QS request requires `qs`; run the dependency checker with `--config` before execution. Multiple raw-count layers are joined in memory and verified against the full cell set before applying decisions.
