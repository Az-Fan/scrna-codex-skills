#!/usr/bin/env python3
import argparse
import json
import shutil
import subprocess
import sys
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
    required = [
        ("pixi.project", config.get("pixi", {}).get("project")),
        ("input.starsolo_dir", config.get("input", {}).get("starsolo_dir")),
        ("input.gtf_file", config.get("input", {}).get("gtf_file")),
        ("samples", config.get("samples")),
        ("output_dir", config.get("output_dir")),
    ]
    missing = [name for name, value in required if not value]
    if missing:
        fail("Missing required configuration: " + ", ".join(missing))

    pixi_project = Path(config["pixi"]["project"]).expanduser().resolve()
    manifest = pixi_project if pixi_project.name == "pixi.toml" else pixi_project / "pixi.toml"
    starsolo = Path(config["input"]["starsolo_dir"]).expanduser().resolve()
    gtf = Path(config["input"]["gtf_file"]).expanduser().resolve()
    output = Path(config["output_dir"]).expanduser().resolve()
    environment = config.get("pixi", {}).get("environment", "default")

    if not manifest.is_file():
        fail(f"pixi.toml not found: {manifest}")
    if not starsolo.is_dir():
        fail(f"STARsolo directory not found: {starsolo}")
    if not gtf.is_file():
        fail(f"GTF not found: {gtf}")
    if not isinstance(config["samples"], list) or not config["samples"]:
        fail("samples must be a non-empty list")

    sample_ids = []
    for sample in config["samples"]:
        if not sample.get("sample_id") or not sample.get("batch_id"):
            fail("Each sample requires sample_id and batch_id")
        sample_ids.append(sample["sample_id"])
        base = starsolo / sample["sample_id"] / "Solo.out"
        files = [
            base / "Gene/filtered/matrix.mtx.gz",
            base / "Gene/filtered/features.tsv.gz",
            base / "Gene/filtered/barcodes.tsv.gz",
            base / "Velocyto/filtered/spliced.mtx.gz",
            base / "Velocyto/filtered/unspliced.mtx.gz",
            base / "Velocyto/filtered/barcodes.tsv.gz",
        ]
        absent = [str(path) for path in files if not path.is_file()]
        if absent:
            fail("Missing STARsolo files:\n" + "\n".join(absent))
    if len(sample_ids) != len(set(sample_ids)):
        fail("sample_id values must be unique")

    configured_pixi = config.get("pixi", {}).get("executable")
    candidates = [configured_pixi, shutil.which("pixi"), str(Path.home() / ".pixi/bin/pixi")]
    pixi = next((str(Path(item).expanduser().resolve()) for item in candidates if item and Path(item).expanduser().is_file()), None)
    if not pixi:
        fail("Existing pixi executable not found; set pixi.executable. The skill will not install it")
    driver = Path(__file__).with_name("calculate_metrics.R").resolve()
    command = [pixi, "run", "--manifest-path", str(manifest), "-e", environment,
               "--", "Rscript", str(driver), str(config_path)]
    plan = {
        "skill": "calculate-scrna-qc-metrics",
        "mode": "execute" if args.execute else "dry-run",
        "pixi_manifest": str(manifest),
        "pixi_environment": environment,
        "starsolo_dir": str(starsolo),
        "gtf_file": str(gtf),
        "samples": config["samples"],
        "output_dir": str(output),
        "command": command,
        "note": "No cells will be filtered and no environment will be modified.",
    }
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    if not args.execute:
        return
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)
    (output / "run_manifest.json").write_text(
        json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
