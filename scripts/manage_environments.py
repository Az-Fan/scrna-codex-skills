#!/usr/bin/env python3
"""Audit or explicitly deploy the versioned shared Pixi environment bundle."""
import argparse
import datetime as dt
import hashlib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BUNDLE = REPO / "environments"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def audit(bundle, target, selected):
    report = {}
    for name in selected:
        spec = bundle["profiles"][name]
        files = {}
        for filename, record in spec["files"].items():
            source = BUNDLE / record["source"]
            if digest(source) != record["sha256"]:
                raise ValueError("Bundle checksum mismatch: " + str(source))
            current = digest(target / name / filename)
            files[filename] = {"expected_sha256": record["sha256"], "actual_sha256": current, "matches": current == record["sha256"]}
        interpreters = {environment: (target / name / ".pixi/envs" / environment / "bin" / ("Rscript" if environment == "default" else "python")).is_file() for environment in spec["environments"]}
        report[name] = {"files": files, "installed_interpreters": interpreters}
    return report


def deploy(bundle, target, selected, force=False):
    # 校验全部源文件后再复制；备份仅包含将被替换的配置，不搬动二进制环境。
    audit(bundle, target, selected)
    conflicts = [target / name / filename for name in selected for filename in bundle["profiles"][name]["files"] if (target / name / filename).exists()]
    if conflicts and not force:
        raise ValueError("Existing configuration requires --force; it will be backed up: " + str(conflicts[0]))
    target.mkdir(parents=True, exist_ok=True)
    backup = None
    if conflicts:
        backup_root = target / ".environment-backups"
        backup_root.mkdir(exist_ok=True)
        backup = Path(tempfile.mkdtemp(prefix=dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ-"), dir=backup_root))
        for path in conflicts:
            destination = backup / path.relative_to(target)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, destination)
    published = []
    try:
        for name in selected:
            for filename, record in bundle["profiles"][name]["files"].items():
                destination = target / name / filename
                destination.parent.mkdir(parents=True, exist_ok=True)
                published.append(destination)
                shutil.copy2(BUNDLE / record["source"], destination)
    except BaseException:
        for path in reversed(published):
            prior = backup / path.relative_to(target) if backup else None
            if prior and prior.is_file():
                shutil.copy2(prior, path)
            else:
                path.unlink(missing_ok=True)
        raise
    return str(backup) if backup else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=Path.home() / "projects/scrna_envs")
    parser.add_argument("--profile", action="append", help="limit to named profiles")
    parser.add_argument("--apply", action="store_true", help="copy versioned configuration; default is read-only audit")
    parser.add_argument("--force", action="store_true", help="back up and replace existing configuration")
    parser.add_argument("--install", action="store_true", help="with --apply, explicitly install locked runtimes and supplemental packages")
    parser.add_argument("--pixi", default=str(Path.home() / ".pixi/bin/pixi"))
    args = parser.parse_args()
    if args.install and not args.apply:
        parser.error("--install requires --apply")
    target = args.target.expanduser().resolve()
    if target == BUNDLE or BUNDLE in target.parents:
        parser.error("Deploy outside the canonical environment bundle")
    bundle = json.loads((BUNDLE / "bundle.json").read_text())
    selected = args.profile or list(bundle["profiles"])
    unknown = set(selected) - set(bundle["profiles"])
    if unknown:
        parser.error("Unknown profiles: " + ", ".join(sorted(unknown)))
    backup = deploy(bundle, target, selected, args.force) if args.apply else None
    if args.install:
        for name in selected:
            project = target / name
            for environment in bundle["profiles"][name]["environments"]:
                subprocess.run([args.pixi, "install", "--frozen", "--run-post-link-scripts", "--manifest-path", str(project / "pixi.toml"), "-e", environment], check=True)
            helper = project / "install_supplemental.R"
            if helper.is_file():
                subprocess.run([args.pixi, "run", "--frozen", "--no-install", "--manifest-path", str(project / "pixi.toml"), "-e", "default", "--", "Rscript", str(helper), "--install"], check=True)
    report = audit(bundle, target, selected)
    for name, record in report.items():
        helper = target / name / "install_supplemental.R"
        if "install_supplemental.R" in bundle["profiles"][name]["files"] and helper.is_file() and record["installed_interpreters"]["default"]:
            if not all(record["files"][filename]["matches"] for filename in ("install_supplemental.R", "supplemental-r.json")):
                record["supplemental"] = {"passed": False, "error": "Supplemental check skipped: helper or package manifest differs from the versioned bundle"}
                continue
            result = subprocess.run([str(target / name / ".pixi/envs/default/bin/Rscript"), str(helper), "--check"], text=True, capture_output=True)
            record["supplemental"] = {"passed": result.returncode == 0, "details": result.stdout.strip(), "error": result.stderr.strip()}
    passed = all(all(x["matches"] for x in r["files"].values()) and all(r["installed_interpreters"].values()) and r.get("supplemental", {}).get("passed", True) for r in report.values())
    print(json.dumps({"mode": "deploy" if args.apply else "audit", "target": str(target), "backup": backup, "passed": passed, "profiles": report}, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
