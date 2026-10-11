"""记录基因检出率，显式筛选低频基因，保留全部细胞。"""
import math
import warnings
import numpy as np
import pandas as pd
from scipy import sparse
from gene_filter_contract import validate_thresholds


def filter_counts(counts, genes, min_cells=1, min_fraction=0):
    validate_thresholds(min_cells, min_fraction)
    counts = sparse.csr_matrix(counts)
    counts.sum_duplicates()
    if counts.shape[0] == 0 or counts.shape[1] == 0 or counts.shape[1] != len(genes):
        raise ValueError("Counts must have non-empty cells and matching gene identifiers")
    if not np.isfinite(counts.data).all() or (counts.data < 0).any() or not np.allclose(counts.data, np.rint(counts.data), rtol=0, atol=1e-8):
        raise ValueError('Expected non-negative integer raw UMI counts')
    threshold = max(min_cells, math.ceil(min_fraction * counts.shape[0]))
    detected = np.asarray((counts > 0).sum(axis=0)).ravel().astype(int)
    retained = detected >= threshold
    audit = pd.DataFrame({'gene': list(genes), 'n_detected_cells': detected,
                          'detection_fraction': detected / counts.shape[0],
                          'total_counts': np.asarray(counts.sum(axis=0)).ravel(),
                          'effective_min_cells': threshold,
                          'status': np.where(detected == 0, 'zero_counts', np.where(retained, 'retained', 'low_detection'))})
    filtered = counts[:, retained].tocsr()
    if not retained.any():
        raise ValueError('Gene filter retained no genes; review thresholds')
    if (np.asarray(filtered.sum(axis=1)).ravel() <= 0).any():
        raise ValueError('Gene filtering produced zero-count cells; review thresholds, no cells were silently removed')
    return filtered, retained, audit


def annotate_top_genes(table, metrics, warning_min_cells=20):
    validate_thresholds(warning_min_cells, min_cells_field="reporting.low_detection_cells")
    metrics = metrics[['gene', 'n_detected_cells', 'detection_fraction', 'total_counts']]
    result = table.merge(metrics, on='gene', how='left', validate='many_to_one', sort=False)
    if result.n_detected_cells.isna().any():
        raise ValueError('Top gene IDs do not match feature audit')
    result['low_detection'] = result.n_detected_cells < warning_min_cells
    return result


def write_top_gene_tables(scores, tpm, metrics, destination, top_n=50, warning_min_cells=20):
    destination.mkdir(parents=True, exist_ok=True)
    rows = []
    for program in scores:
        for rank, (gene, score) in enumerate(scores[program].nlargest(top_n).items(), 1):
            rows.append({'program': program, 'rank': rank, 'gene': gene, 'score': score,
                         'tpm_weight': tpm.loc[gene, program], 'positive_score': bool(score > 0)})
    ranked = annotate_top_genes(pd.DataFrame(rows), metrics, warning_min_cells)
    ranked.to_csv(destination / 'top_genes.tsv', sep='\t', index=False)
    summary = ranked.groupby('program', sort=False).agg(top_genes=('gene', 'size'), low_detection_genes=('low_detection', 'sum')).reset_index()
    summary['warning_min_cells'] = warning_min_cells
    summary.to_csv(destination / 'top_gene_detection_summary.tsv', sep='\t', index=False)
    if ranked.low_detection.any():
        warnings.warn(f'{int(ranked.low_detection.sum())} top-gene entries detected in fewer than {warning_min_cells} cells; inspect top_gene_detection_summary.tsv', UserWarning)
    rows = []
    for program in tpm:
        for rank, (gene, weight) in enumerate(tpm[program].nlargest(top_n).items(), 1):
            rows.append({'program': program, 'rank': rank, 'gene': gene, 'tpm_weight': weight, 'score': scores.loc[gene, program]})
    annotate_top_genes(pd.DataFrame(rows), metrics, warning_min_cells).to_csv(destination / 'top_genes_by_tpm.tsv', sep='\t', index=False)
    return ranked
