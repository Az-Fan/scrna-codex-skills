#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


def fail(message):
    raise SystemExit(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))

    pixi_config = config.get("pixi", {})
    project_value = pixi_config.get("project")
    output_value = config.get("output_dir")
    if not project_value or not output_value:
        fail("pixi.project and output_dir are required")
    project = Path(project_value).expanduser().resolve()
    manifest = project if project.name == "pixi.toml" else project / "pixi.toml"
    if not manifest.is_file():
        fail(f"pixi.toml not found: {manifest}")

    input_config = config.get("input", {})
    input_type = input_config.get("type", "auto")
    if input_type == "auto":
        input_type = "seurat" if input_config.get("object") else "starsolo"
    if input_type not in {"seurat", "starsolo"}:
        fail("input.type must be seurat, starsolo, or auto")
    key = "object" if input_type == "seurat" else "starsolo_dir"
    if not input_config.get(key):
        fail(f"input.{key} is required for {input_type} mode")
    primary = Path(input_config[key]).expanduser().resolve()
    if input_type == "seurat" and not primary.is_file():
        fail(f"Seurat object not found: {primary}")
    if input_type == "starsolo" and not primary.is_dir():
        fail(f"STARsolo directory not found: {primary}")

    gtf_value = input_config.get("gtf_file")
    gtf = Path(gtf_value).expanduser().resolve() if gtf_value else None
    if gtf and not gtf.is_file():
        fail(f"Configured GTF not found: {gtf}")
    samples = config.get("samples", [])
    if input_type == "starsolo" and not samples:
        fail("STARsolo mode requires samples with sample_id and batch_id")
    for sample in samples:
        if not sample.get("sample_id") or not sample.get("batch_id"):
            fail("Each STARsolo sample requires sample_id and batch_id")
        base = primary / sample["sample_id"] / "Solo.out/Gene/filtered"
        missing = [str(base / name) for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz") if not (base / name).is_file()]
        if missing:
            fail("Missing primary STARsolo files:\n" + "\n".join(missing))

    configured_pixi = pixi_config.get("executable")
    candidates = [configured_pixi, shutil.which("pixi"), str(Path.home() / ".pixi/bin/pixi")]
    pixi = next((str(Path(x).expanduser().resolve()) for x in candidates if x and Path(x).expanduser().is_file()), None)
    if not pixi:
        fail("Existing pixi executable not found; set pixi.executable. The skill will not install it")
    environment = pixi_config.get("environment", "default")
    parallel_config = config.get("parallel", {})
    if not isinstance(parallel_config, dict):
        fail("parallel must be an object")
    workers = parallel_config.get("workers", 1)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        fail("parallel.workers must be a positive integer")
    ambient_config = config.get("ambient_rna", {})
    if not isinstance(ambient_config, dict):
        fail("ambient_rna must be an object")
    ambient_method = ambient_config.get("method", "decontx")
    if not isinstance(ambient_method, str):
        fail("ambient_rna.method must be decontx or skip")
    ambient_method = ambient_method.lower()
    if ambient_method not in {"decontx", "skip"}:
        fail("ambient_rna.method must be decontx or skip")
    cluster_column = ambient_config.get("cluster_column")
    if cluster_column is not None and (not isinstance(cluster_column, str) or not cluster_column.strip()):
        fail("ambient_rna.cluster_column must be null or a non-empty metadata column name")
    output = Path(output_value).expanduser().resolve()
    provenance = output / "_provenance"
    script = Path(__file__).resolve()
    driver_candidates = [script.with_name("calculate_metrics.R")]
    driver_candidates.extend(parent / "toolkit" / "R" / "calculate_metrics.R" for parent in script.parents)
    driver = next((path for path in driver_candidates if path.is_file()), None)
    if driver is None:
        fail("calculate_metrics.R not found; install a self-contained build or run from the source repository")
    command = [pixi, "run", "--manifest-path", str(manifest), "-e", environment,
               "--", "Rscript", str(driver), str(config_path)]
    plan = {
        "skill": "02-scrna-calculate-qc-metrics",
        "mode": "execute" if args.execute else "dry-run",
        "pixi_manifest": str(manifest),
        "pixi_environment": environment,
        "input_type": input_type,
        "input": str(primary),
        "gtf_file": str(gtf) if gtf else None,
        "samples": samples,
        "parallel_workers": workers,
        "ambient_rna_method": ambient_method,
        "ambient_rna_cluster_column": cluster_column,
        "output_dir": str(output),
        "provenance_dir": str(provenance),
        "command": command,
        "note": "Unavailable optional metrics are recorded as skipped; no cells or environments are changed.",
    }
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    if not args.execute:
        return
    output.mkdir(parents=True, exist_ok=True)
    provenance.mkdir(parents=True, exist_ok=True)
    execution_env = os.environ.copy()
    if workers > 1:
        for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
            execution_env[variable] = "1"
    with (provenance / "run.log").open("a", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=execution_env)
    plan["exit_status"] = result.returncode
    (provenance / "run_manifest.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
    if result.returncode:
        fail(f"QC metrics failed (exit {result.returncode}); see {provenance / 'run.log'}")


if __name__ == "__main__":
    main()
