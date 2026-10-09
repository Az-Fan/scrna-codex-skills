#!/usr/bin/env python3
"""Tiny artificial resources for GRN execution tests; no biological evidence."""
import csv
from pathlib import Path
import sys
import numpy as np
import pyarrow as pa
import pyarrow.feather as feather


def main():
    root = Path(sys.argv[1])
    root.mkdir(parents=True, exist_ok=True)
    genes = ["mt-Nd1", "mt-Co1", "Kdr", "Pecam1", "Cdh5", "Col1a1", "Col3a1", "Dcn", "Actb", "Gapdh"] + [f"Gene{i}" for i in range(1, 11)]
    rng = np.random.default_rng(777)
    rankings = np.array([rng.permutation(len(genes)) for _ in range(50)], dtype=np.int16)
    motifs = [f"synthetic_motif_{i}" for i in range(len(rankings))]
    table = {gene: pa.array(rankings[:, i]) for i, gene in enumerate(genes)}
    table["motifs"] = pa.array(motifs)
    feather.write_feather(pa.table(table), root / "toy.genes_vs_motifs.rankings.feather", version=2)
    # These are fixture expression columns assigned synthetic TF roles, not a TF catalogue.
    (root / "tfs.txt").write_text("Kdr\nPecam1\nCol1a1\nDcn\n")
    with (root / "annotations.tbl").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#motif_id", "gene_name", "motif_similarity_qvalue", "orthologous_identity", "description"])
        for motif in motifs:
            for tf in ["Kdr", "Pecam1", "Col1a1", "Dcn"]:
                writer.writerow([motif, tf, 0, 1, "synthetic execution fixture"])


if __name__ == "__main__":
    main()
