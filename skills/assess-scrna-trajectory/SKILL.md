---
name: assess-scrna-trajectory
description: Assess whether a scRNA-seq population supports trajectory analysis and, after biological confirmation, run hypothesis-led lineage, pseudotime, dynamic-gene, and gene-module analyses. Use when proposing roots and endpoints, evaluating continuous state structure, migrating a lineage workflow to a new tissue, or preventing inherited lineage assumptions from being applied blindly.
---

# Assess scRNA Trajectory

Require a biological lineage hypothesis before fitting pseudotime.

## Workflow

1. Define the target population, biological process, expected direction, and alternative explanations.
2. Inspect sample, condition, batch, cell-cycle, QC, and annotation structure in candidate embeddings.
3. Test whether cells form a sufficiently connected and continuous manifold.
4. Propose candidate roots, terminals, and branches using markers, time information, or established biology.
5. Produce a trajectory configuration and pause for confirmation before committing roots or endpoints.
6. Run the selected method with sensitivity checks across dimensions, neighbors, seeds, and reasonable roots.
7. Identify dynamic genes and modules while controlling for sample and condition where possible.
8. Validate pseudotime against markers, sample coverage, known ordering, and alternative embeddings.

## Guardrails

- Never carry lineage names, roots, or endpoints from another dataset.
- Do not infer temporal progression solely from a two-dimensional UMAP.
- Do not force a trajectory across disconnected cell identities.
- Report when condition or batch is indistinguishable from pseudotime.

Read [references/trajectory-checklist.md](references/trajectory-checklist.md) before fitting a lineage.
