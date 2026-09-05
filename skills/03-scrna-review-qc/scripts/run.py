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
    config_path = Path(args.config).expanduser().resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    pixi_cfg = config.get("pixi", {})
    project_value = pixi_cfg.get("project")
    object_value = config.get("input", {}).get("object")
    sample_column = config.get("metadata", {}).get("sample")
    output_value = config.get("output_dir")
    if not all((project_value, object_value, sample_column, output_value)):
        fail("pixi.project, input.object, metadata.sample, and output_dir are required")
    project = Path(project_value).expanduser().resolve()
    manifest = project if project.name == "pixi.toml" else project / "pixi.toml"
    object_path = Path(object_value).expanduser().resolve()
    if not manifest.is_file():
        fail(f"pixi.toml not found: {manifest}")
    if not object_path.is_file() or object_path.suffix.lower() not in {".rds", ".qs"}:
        fail(f"Seurat RDS/QS object not found: {object_path}")
    configured = pixi_cfg.get("executable")
    candidates = [configured, shutil.which("pixi"), str(Path.home() / ".pixi/bin/pixi")]
    pixi = next((str(Path(x).expanduser().resolve()) for x in candidates if x and Path(x).expanduser().is_file()), None)
    if not pixi:
        fail("Existing pixi executable not found; set pixi.executable. This skill will not install it")
    output_config = config.get("output", {})
    if not isinstance(output_config, dict):
        fail("output must be an object")
    detail_level = output_config.get("detail_level", "compact")
    if not isinstance(detail_level, str) or detail_level.lower() not in {"compact", "full"}:
        fail("output.detail_level must be compact or full")
    detail_level = detail_level.lower()
    output = Path(output_value).expanduser().resolve()
    script = Path(__file__).resolve()
    driver_candidates = [script.with_name("review_qc.R")]
    driver_candidates.extend(parent / "toolkit" / "R" / "review_qc.R" for parent in script.parents)
    driver = next((path for path in driver_candidates if path.is_file()), None)
    if driver is None:
        fail("review_qc.R not found; install a self-contained build or run from the source repository")
    environment = pixi_cfg.get("environment", "default")
    command = [pixi, "run", "--manifest-path", str(manifest), "-e", environment,
               "--", "Rscript", str(driver), str(config_path)]
    plan = {"skill": "03-scrna-review-qc", "mode": "execute" if args.execute else "dry-run",
            "pixi_manifest": str(manifest), "pixi_environment": environment,
            "input_object": str(object_path), "sample_column": sample_column,
            "output_detail_level": detail_level,
            "output_dir": str(output), "command": command,
            "note": "Review only: no filtering, object mutation, package installation, or environment changes."}
    print(json.dumps(plan, ensure_ascii=False, indent=2))
    if not args.execute:
        return
    output.mkdir(parents=True, exist_ok=True)
    provenance = output / "_provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    with (provenance / "run.log").open("a", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    plan["exit_status"] = result.returncode
    (provenance / "run_manifest.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
    if result.returncode:
        fail(f"QC review failed (exit {result.returncode}); see {provenance / 'run.log'}")


if __name__ == "__main__":
    main()
