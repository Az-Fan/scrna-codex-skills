---
name: compare-cell-communication
description: Prepare, run, visualize, and compare ligand-receptor cell-communication analyses for scRNA-seq data, primarily with CellChat-compatible workflows. Use when constructing condition-specific communication networks, comparing interaction counts or strengths, focusing on selected ligand-receptor pairs, or auditing whether sample and cell-type support is adequate.
---

# Compare Cell Communication

Parameterize grouping explicitly and distinguish descriptive networks from replicated inference.

## Workflow

1. Validate cell-type, sample, condition, and optional batch columns.
2. Tabulate cells per sample, condition, and cell type; flag missing or sparse sender/receiver populations.
3. Select the species-compatible interaction database and record its version.
4. Run a pooled descriptive network only when that matches the question.
5. For condition comparison, construct comparable objects with consistent labels and preprocessing.
6. Compare interaction number and weight, pathway-level networks, centrality, and selected ligand-receptor pairs.
7. Check whether apparent differences are driven by cell abundance or a single sample.
8. Save machine-readable interaction tables alongside plots.

## Guardrails

- Never retain inherited labels such as `early/late` when the project comparison is different.
- Do not claim differential communication from pooled groups without acknowledging sample replication.
- Do not interpret database presence as evidence that both genes are detectably expressed.
- Preserve group and cell-type mappings used in every comparison.

Read [references/communication-design.md](references/communication-design.md) before comparing groups.

