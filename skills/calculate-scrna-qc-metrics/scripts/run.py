#!/usr/bin/env python3
import argparse
import json
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
    output = Path(output_value).expanduser().resolve()
    driver = Path(__file__).with_name("calculate_metrics.R").resolve()
    command = [pixi, "run", "--manifest-path", str(manifest), "-e", environment,
               "--", "Rscript", str(driver), str(config_path)]
    plan = {
        "skill": "calculate-scrna-qc-metrics",
        "mode": "execute" if args.execute else "dry-run",
        "pixi_manifest": str(manifest),
        "pixi_environment": environment,
        "input_type": input_type,
        "input": str(primary),
        "gtf_file": str(gtf) if gtf else None,
        "samples": samples,
        "output_dir": str(output),
        "command": command,
        "note": "Unavailable optional metrics are recorded as skipped; no cells or environments are changed.",
    }
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    if not args.execute:
        return
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)
    (output / "run_manifest.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
