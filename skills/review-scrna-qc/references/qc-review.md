# QC review contract

Recognized aliases, in priority order:

| Canonical metric | Metadata aliases | Candidate direction |
|---|---|---|
| n_genes | n_genes, nFeature_RNA | lower and upper |
| n_UMIs | n_UMIs, nCount_RNA | lower and upper |
| mito_frac | mito_frac, percent.mt | upper |
| nuclear_frac | nuclear_frac | upper |
| ambient_frac | ambient_frac_decontx, ambient_frac | upper |
| doublet_score | doublet_score | upper |
| hbb_score | hbb_score | upper |
| s_score | s_score, S.Score | descriptive only |
| g2m_score | g2m_score, G2M.Score | descriptive only |

Fractions and percentages are reviewed on their stored scale; the script does not silently rescale them.
Columns containing no finite numeric values are unavailable even if the column name exists.

Candidate bounds combine global 1st/99th percentiles with median ± 3 MAD. Explicit values in `thresholds` override generated candidates. Generated thresholds are screening suggestions, not accepted biological decisions.

The generated bounds are global exploratory heuristics, not validated universal QC cutoffs. Review sample distributions, tissue biology, chemistry, sequencing depth, and expected rare populations before accepting any threshold. Treat DecontX and synthetic-doublet scores as model outputs rather than ground-truth labels.

The joint candidate scenario applies all available filtering bounds only to estimate per-sample and per-group retention. It never subsets or saves a filtered object.

When coordinates exist, generate one UMAP for every available QC metric. Auto-detect `seurat_clusters`/`cluster` and common annotation columns for grouped violin plots. Record every generated or skipped plot in `plot_status.tsv`.
