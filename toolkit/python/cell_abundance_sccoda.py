#!/usr/bin/env python3
"""Run one scCODA comparison from sample-by-cell-type counts."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pertpy as pt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--numerator", required=True)
    parser.add_argument("--denominator", required=True)
    parser.add_argument("--comparison-id", required=True)
    parser.add_argument("--reference", default="automatic")
    parser.add_argument("--covariates", default="")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--num-samples", type=int, default=10000)
    parser.add_argument("--num-warmup", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--absence-threshold", type=float, default=0.05)
    parser.add_argument("--save-posterior", action="store_true")
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0)
    if not counts.index.is_unique or not metadata.index.is_unique:
        raise ValueError("sample identifiers must be unique")
    missing = counts.index.difference(metadata.index)
    if len(missing):
        raise ValueError(f"metadata missing samples: {', '.join(map(str, missing))}")
    metadata = metadata.loc[counts.index].copy()
    # xarray interprets a named AnnData observation index as a dimension name.
    # scCODA expects the observation dimension itself to remain named ``obs``.
    metadata.index.name = None
    keep = metadata[args.condition].astype(str).isin([args.denominator, args.numerator])
    counts = counts.loc[keep].copy()
    metadata = metadata.loc[keep].copy()
    if counts.shape[0] < 4 or counts.shape[1] < 2:
        raise ValueError("scCODA requires at least four samples and two cell types")
    if (counts.to_numpy() < 0).any() or not np.allclose(counts, np.round(counts)):
        raise ValueError("scCODA counts must be non-negative integers")

    metadata["contrast_group"] = pd.Categorical(
        np.where(metadata[args.condition].astype(str) == args.numerator, "numerator", "denominator"),
        categories=["denominator", "numerator"],
        ordered=True,
    )
    covariates = [value for value in args.covariates.split(",") if value]
    required = ["contrast_group", *covariates]
    for column in required:
        if column not in metadata.columns:
            raise ValueError(f"metadata column not found: {column}")
    formula = " + ".join([*covariates, "contrast_group"])

    adata = ad.AnnData(
        X=counts.to_numpy(dtype=np.int64),
        obs=metadata[required].copy(),
        var=pd.DataFrame(index=counts.columns.astype(str)),
    )
    model = pt.tl.Sccoda()
    mdata = model.load(adata, type="sample_level", covariate_obs=required)
    mdata = model.prepare(
        mdata,
        formula=formula,
        reference_cell_type=args.reference,
        automatic_reference_absence_threshold=args.absence_threshold,
    )
    model.run_nuts(
        mdata,
        num_samples=args.num_samples,
        num_warmup=args.num_warmup,
        rng_key=args.seed,
    )

    effects = model.get_effect_df(mdata).reset_index()
    effects.to_csv(args.output_dir / "sccoda_effects.tsv", sep="\t", index=False)
    credible = model.credible_effects(mdata, est_fdr=args.fdr)
    credible_map = {(str(a), str(b)): bool(value) for (a, b), value in credible.items()}
    condition_rows = effects[effects["Covariate"].astype(str).str.contains("contrast_group")].copy()
    actual_reference = str(mdata["coda"].uns["scCODA_params"]["reference_cell_type"])

    rows = []
    for _, row in condition_rows.iterrows():
        covariate = str(row["Covariate"])
        cell_type = str(row["Cell Type"])
        is_credible = credible_map.get((covariate, cell_type), False)
        effect = float(row.get("log2-fold change", np.nan))
        rows.append(
            {
                "method": "sccoda",
                "comparison_id": args.comparison_id,
                "feature_type": "cell_type",
                "feature_id": cell_type,
                "annotation": cell_type,
                "numerator": args.numerator,
                "denominator": args.denominator,
                "effect": effect,
                "effect_scale": "sccoda_log2_fold_change_relative_to_reference",
                "ci_lower": float(row.get("ETI 89% lower", np.nan)),
                "ci_upper": float(row.get("ETI 89% upper", np.nan)),
                "p_value": np.nan,
                "adjusted_p_value": np.nan,
                "posterior_probability": float(row.get("Inclusion probability", np.nan)),
                "evidence_type": "posterior_inclusion_probability",
                "credible": is_credible,
                "significant": is_credible,
                "tested": True,
                "reference": actual_reference,
                "status": "completed",
                "note": "Effect is relative to the selected reference; p-values are not defined.",
            }
        )
    standardized = pd.DataFrame(rows)
    standardized.to_csv(args.output_dir / "all_results.tsv", sep="\t", index=False)
    standardized[standardized["significant"]].to_csv(
        args.output_dir / "significant_results.tsv", sep="\t", index=False
    )

    diagnostics_status = "not_available"
    try:
        import arviz as az

        inference = model.make_arviz(mdata)
        diagnostics = az.summary(inference)
        diagnostics.to_csv(args.output_dir / "posterior_diagnostics.tsv", sep="\t")
        diagnostics_status = "completed"
        if args.save_posterior:
            inference.to_netcdf(args.output_dir / "posterior.nc")
    except Exception as exc:  # diagnostics must not discard valid model effects
        (args.output_dir / "POSTERIOR_DIAGNOSTICS_WARNING.txt").write_text(str(exc) + "\n", encoding="utf-8")
        diagnostics_status = "failed"

    record = {
        "method": "scCODA via pertpy",
        "comparison_id": args.comparison_id,
        "formula": formula,
        "requested_reference": args.reference,
        "actual_reference": actual_reference,
        "estimated_fdr": args.fdr,
        "num_samples": args.num_samples,
        "num_warmup": args.num_warmup,
        "seed": args.seed,
        "diagnostics_status": diagnostics_status,
        "versions": {name: version(name) for name in ["pertpy", "anndata", "mudata", "numpyro", "jax"]},
    }
    (args.output_dir / "model_record.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
