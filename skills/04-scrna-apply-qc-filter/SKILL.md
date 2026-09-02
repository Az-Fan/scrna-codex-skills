---
name: 04-scrna-apply-qc-filter
description: Apply an explicitly approved cell-level QC decision table to a Seurat RDS/QS object, preserving raw counts and producing an auditable filtered handoff object, per-cell decisions, sample retention summaries, and provenance. Use after QC review and human approval; do not use to invent thresholds or make filtering decisions.
---

# Apply Approved scRNA QC Filter

Apply an existing, reviewed decision rather than deriving a new threshold.

## Required inputs

- A Seurat RDS/QS object containing raw counts.
- A cell-level TSV or TSV.GZ decision table covering exactly the same cell IDs.
- The cell-ID column, boolean columns that must all be true, and optional boolean columns for which any true value excludes a cell.
- The sample column and optional condition column used for retention audits.
- An explicit approval record and expected retained-cell count.
- A new output directory.

The decision table may be produced by project code or a reviewed QC workflow. Reuse that table; do not recreate its scientific rules in this skill.

## Workflow

1. Read [references/filter-contract.md](references/filter-contract.md) and create a config from [references/config.example.json](references/config.example.json).
2. Run `scripts/run.py --config <config>` and review the resolved object, decision table, approval record, expected cell count, output directory, and command.
3. Execute with `scripts/run.py --config <config> --execute` only after those fields match the approved decision.
4. Review per-sample and per-condition retention before handing the filtered object to preprocessing.

If writing the filtered object may exceed 10 minutes or outlive the remote session, read [references/long-running-execution.md](references/long-running-execution.md) and launch the confirmed command with `scripts/run_in_tmux.py`.

## Guarantees

- Require `approval.status=approved`; never infer approval from a populated threshold or decision table.
- Require exact cell-set agreement and reject duplicated cell IDs or missing decision values.
- Refuse to overwrite an existing filtered object.
- Preserve the selected raw-count assay and rebuild a counts-only Seurat handoff object without inherited reductions or graphs.
- Atomically save and reread the object before finalizing it.
- Preserve the complete per-cell decision table and write sample, condition, and reason summaries.
- Stop when the retained count differs from the explicitly approved expected count.
