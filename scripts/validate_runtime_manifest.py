#!/usr/bin/env python3
"""Validate canonical runtime declarations and source/build separation."""
import json
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
data=json.loads((REPO/"release/runtime-manifest.json").read_text(encoding="utf-8"))
skills=data["skills"]
actual={
    p.name
    for p in (REPO/"skills").iterdir()
    if p.is_dir() and p.name[:2].isdigit() and p.name[2:].startswith("-scrna-")
}
errors=[]
if data.get("schema_version") != 2: errors.append(f"unsupported runtime manifest schema: {data.get('schema_version')}")
if set(skills)!=actual: errors.append(f"skill set mismatch: manifest={sorted(skills)} source={sorted(actual)}")
if list(skills) != sorted(skills): errors.append("skills must be declared in numeric execution order")
for name,spec in skills.items():
    scripts=REPO/"skills"/name/"scripts"
    if not (scripts/"run.py").is_file(): errors.append(f"{name}: scripts/run.py missing")
    for language,folder,suffix in (("python",REPO/"toolkit/python",".py"),("r",REPO/"toolkit/R",".R")):
        for filename in spec.get(language,[]):
            if not filename.endswith(suffix): errors.append(f"{name}: invalid {language} runtime name {filename}")
            if not (folder/filename).is_file(): errors.append(f"{name}: canonical runtime missing: {folder/filename}")
            if (scripts/filename).exists(): errors.append(f"{name}: generated runtime committed in source: {filename}")
    references=REPO/"skills"/name/"references"
    for filename in spec.get("references",[]):
        source=REPO/"toolkit/references"/filename
        if not filename.endswith(".md"): errors.append(f"{name}: invalid shared reference name {filename}")
        if not source.is_file(): errors.append(f"{name}: canonical reference missing: {source}")
        if (references/filename).exists(): errors.append(f"{name}: generated reference committed in source: {filename}")
if errors: raise SystemExit("\n".join(errors))
print(f"VALID runtime manifest: {len(skills)} skills")
