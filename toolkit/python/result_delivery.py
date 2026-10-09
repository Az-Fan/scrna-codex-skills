#!/usr/bin/env python3
"""Build portable result navigation without changing scientific artifacts.

Snapshots use filesystem generation metadata, never read large count/object files.
Only changed files and explicitly inherited review stages are current deliverables.
"""
import csv
import datetime as dt
import json
import re
from pathlib import Path
from urllib.parse import quote

INDEX = "RESULTS.md"
REGISTRY = "result_delivery.json"
TECHNICAL = {"_provenance", "exchange", "resource_cache"}
FIGURES = {".png", ".pdf", ".svg", ".jpg", ".jpeg", ".html"}
OBJECTS = {".rds", ".qs", ".h5ad", ".h5", ".mtx"}
SUCCESS = {"completed", "complete", "computed", "generated", "ok", "available", "success", "empty", "not_requested"}
FAILURES = {"failed", "error", "invalid_design", "missing_dependency", "skipped_low_replicates", "partial"}


def snapshot(output):
    root = Path(output)
    result = {}
    if root.exists():
        for path in root.rglob("*"):
            if path.is_file() and not path.is_symlink() and path.name != INDEX:
                stat = path.stat()
                result[path.relative_to(root).as_posix()] = [stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]
    return result


def scientific(name):
    return not any(part in TECHNICAL for part in Path(name).parts) and Path(name).name != INDEX


def legacy_review_files(skill, names):
    """Recognize only documented review artifacts; filenames do not prove lineage."""
    if skill.startswith("08-"):
        patterns = (r"annotation_review\.tsv", r"cluster_markers\.tsv",
                    r"(?:cluster_umap|cluster_sample_umap|canonical_marker_dotplot)\.(?:png|pdf)",
                    r"clustered_object\.(?:qs|rds)")
    else:
        patterns = (r"scenario_(?:summary|cluster_similarity)\.tsv",
                    r"preprocessed_clustered_object\.(?:qs|rds)",
                    r"[A-Za-z0-9_.-]+_(?:resolution_stability|umap_clusters_by_resolution|clustree_resolution|elbow)\.(?:tsv|png|pdf)")
    return {name for name in names if len(Path(name).parts) == 1
            and any(re.fullmatch(pattern, name) for pattern in patterns)}


def label(value):
    return str(value).replace("\n", " ").replace("|", "\\|").replace("[", "\\[").replace("]", "\\]")


def link(name, prefix=""):
    return f"[{label(name)}]({quote(prefix + name, safe='/')})"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_evidence(root, names):
    evidence = []
    broken = False
    # Interpret small status/audit tables, never complete expression/result tables.
    for name in sorted(names):
        path = root / name
        if path.suffix != ".tsv" or not any(key in path.name for key in ("status", "audit", "coverage", "summary", "recommendation", "method_runs", "skipped_metrics")):
            if path.name == "recommendation_status.json":
                try:
                    recommendation = json.loads(path.read_text(encoding="utf-8"))
                    if recommendation.get("status") != "resolved":
                        reasons = recommendation.get("reasons", [])
                        evidence.append((name, "recommendation=" + str(recommendation.get("status")) + "; " + str(reasons), True, ""))
                except (OSError, ValueError) as exc:
                    evidence.append((name, f"Could not read recommendation: {exc}", True, ""))
                    broken = True
            continue
        if path.stat().st_size > 2_000_000:
            continue
        try:
            table = rows(path)
        except (OSError, UnicodeError, csv.Error) as exc:
            evidence.append((name, f"Could not read status evidence: {exc}", False, ""))
            broken = True
            continue
        for row in table:
            identity = next((row.get(key) for key in ("task_id", "comparison_id", "scenario", "metric", "family", "signature", "symbol", "sample") if row.get(key)), "")
            scope = row.get("task_id") or row.get("comparison_id") or ""
            messages = []
            bad = False
            for key in ("status", "enrichment_status"):
                value = str(row.get(key, "")).lower()
                if value and value not in SUCCESS:
                    messages.append(f"{key}={value}")
                    bad |= value in FAILURES or value.startswith(("failed", "invalid", "missing_dependency"))
                elif value == "empty":
                    messages.append(f"{key}=empty (no eligible results)")
            for key in ("low_replication_warning", "low_coverage"):
                if str(row.get(key, "")).lower() in {"true", "1"}:
                    messages.append(f"{key}=true")
            for key in ("found_in_object", "evidence_complete", "passes_min_cells"):
                if str(row.get(key, "")).lower() in {"false", "0"}:
                    messages.append(f"{key}=false")
            if row.get("missing_genes") and row["missing_genes"] != "NA":
                messages.append("missing_genes=" + row["missing_genes"])
            if messages:
                detail = next((row.get(key) for key in ("reason", "message", "interpretation", "notes") if row.get(key) and row[key] != "NA"), "")
                text = "; ".join(filter(None, [identity, *messages, detail]))
                evidence.append((name, text, bad, scope))
    return evidence, broken


def category(name):
    path = Path(name)
    if path.suffix.lower() in FIGURES:
        return "Main figures"
    if path.suffix.lower() in OBJECTS or path.name.endswith(".mtx.gz"):
        return "Downstream objects and matrices"
    if any(word in path.name.lower() for word in ("decision", "review", "threshold", "recommendation")):
        return "Review and decisions"
    return "Complete tables and summaries"


PURPOSES = {
    "all_comparisons.tsv": "Complete differential results across comparisons",
    "all_genes.tsv": "Complete eligible and ineligible gene results",
    "task_status.tsv": "Task outcomes and reasons for skipped or failed analyses",
    "design_audit.tsv": "Sample design and inference audit",
    "gene_status.tsv": "Requested gene availability",
    "plot_status.tsv": "Generated and skipped plot families",
    "signature_coverage.tsv": "Gene signature coverage and missing genes",
    "threshold_review.tsv": "QC thresholds for human review",
    "annotation_review.tsv": "Cluster annotation decisions for confirmation",
    "recommendation_status.json": "Integration recommendation evidence and review requirements",
    "recommendation.md": "Integration recommendation and interpretation",
    "cell_filter_decisions.tsv.gz": "Complete per-cell inclusion decisions",
    "all_method_results.tsv": "Complete abundance results across methods",
    "sample_cell_counts.tsv": "Per-sample cell counts used by the analysis",
    "sample_cell_proportions.tsv": "Per-sample composition within the stated denominator",
}


def purpose(name):
    base = Path(name).name
    if base in PURPOSES:
        return PURPOSES[base]
    title = category(name)
    if title == "Main figures":
        return "Figure: " + Path(name).stem.replace("_", " ")
    if title == "Downstream objects and matrices":
        return "Analysis object or matrix for downstream use"
    if title == "Review and decisions":
        return "Human review or confirmed decision record"
    if "significant" in base or "top_" in base:
        return "Selected results; see complete tables for the full tested set"
    if "audit" in base or "status" in base or "coverage" in base:
        return "Interpretation and completeness audit"
    if "summary" in base:
        return "Scientific summary"
    return "Complete scientific table or supporting result"


def write_index(root, files, retained, status, evidence, run_id, skill, tasks=(), old_count=0, provenance="_provenance/", project=None, legacy_retained=()):
    lines = ["# Results", "", f"Status: **{status}**", "", f"Project: **{label(project or skill)}**", "", f"Analysis: `{skill}`", ""]
    if status == "awaiting confirmation":
        lines += ["Review the decision material before confirming the next stage.", ""]
    if old_count:
        lines += [f"{old_count} unchanged files from earlier runs remain on disk and are excluded from this run's deliverables.", ""]
    if evidence:
        lines += ["## Limitations and skipped outputs", ""]
        grouped = {}
        for name, message, _, _ in evidence:
            grouped.setdefault(name, set()).add(message)
        for name, messages in grouped.items():
            lines.append(f"- {link(name)}: " + label("; ".join(sorted(messages))))
        lines.append("")
    if tasks:
        lines += ["## Tasks", ""] + [f"- {link(name + '/' + INDEX)} — {task_status}" for name, task_status in tasks] + [""]
    for title in ("Main figures", "Complete tables and summaries", "Review and decisions", "Downstream objects and matrices"):
        names = sorted(name for name in files if category(name) == title)
        if names:
            lines += [f"## {title}", ""]
            for name in names:
                suffix = (" — legacy review material; stage provenance unverified" if name in legacy_retained
                          else " — retained from the preceding review stage" if name in retained else "")
                lines.append(f"- {link(name)} — {purpose(name)}{suffix}")
            lines.append("")
    if not files and not tasks:
        lines += ["No scientific deliverables were written by this run.", ""]
    lines += [f"Execution records and troubleshooting: [_provenance/]({provenance}).", ""]
    (root / INDEX).write_text("\n".join(lines), encoding="utf-8")


def finish(output, before, skill, config, returncode, run_id=None):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    technical = root / "_provenance"
    technical.mkdir(exist_ok=True)
    run_id = run_id or dt.datetime.now(dt.timezone.utc).isoformat()
    after = snapshot(root)
    changed = {name for name, signature in after.items() if before.get(name) != signature}
    current = {name for name in changed if scientific(name)}
    action = config.get("workflow", {}).get("action", "run")
    inherit = (skill.startswith("06-") and action == "finalize_resolution") or (skill.startswith("08-") and action == "apply_confirmed") or (skill.startswith("16-") and action == "consensus") or (skill.startswith("17-") and action == "infer")
    retained = set()
    legacy_retained = set()
    if inherit:
        try:
            prior = json.loads((technical / REGISTRY).read_text(encoding="utf-8"))
            if not isinstance(prior, dict) or prior.get("skill") != skill:
                raise ValueError("Delivery registry is not for this skill")
            for field in ("current_files", "retained_stage_files", "legacy_unverified_files"):
                values = prior.get(field, [])
                if not isinstance(values, list) or any(not isinstance(name, str) or Path(name).is_absolute()
                                                     or ".." in Path(name).parts for name in values):
                    raise ValueError("Delivery registry has invalid artifact paths")
            candidates = set(prior.get("current_files", [])) | set(prior.get("retained_stage_files", []))
            legacy_retained = set(prior.get("legacy_unverified_files", []))
        except (OSError, ValueError, TypeError):
            # Only known review filenames are eligible; never certify their history.
            candidates = legacy_review_files(skill, before)
            legacy_retained = candidates
        retained = {name for name in candidates if name in after and name not in current and scientific(name)}
    files = current | retained
    evidence, broken = read_evidence(root, changed)
    workflow = technical / "workflow_state.json"
    pending = skill.startswith("08-") and action == "prepare_review"
    if skill.startswith(("06-", "16-", "17-")) and "_provenance/workflow_state.json" in changed:
        try:
            pending |= "awaiting" in str(json.loads(workflow.read_text()).get("status", ""))
        except (OSError, ValueError):
            broken = True
    partial = broken or any(item[2] for item in evidence) or any(Path(name).name == "ERROR.txt" for name in current)
    status = "failed" if returncode else "partial" if partial else "awaiting confirmation" if pending else "completed"
    old = {name for name in before if name in after and scientific(name)} - files
    task_entries = []
    task_roots = sorted({str(Path(name).parts[0] + "/" + Path(name).parts[1]) for name in files if len(Path(name).parts) > 2 and Path(name).parts[0] == "comparisons"})
    for name in task_roots:
        prefix = name + "/"
        task_files = {f[len(prefix):] for f in files if f.startswith(prefix)}
        task_retained = {f[len(prefix):] for f in retained if f.startswith(prefix)}
        task_evidence = []
        for source, message, bad, scope in evidence:
            if source.startswith(prefix):
                task_evidence.append((source[len(prefix):], message, bad, scope))
            elif scope and (name.split("/")[-1] == scope or name.split("/")[-1].startswith(scope + "__")):
                task_evidence.append(("../../" + source, message, bad, scope))
        task_status = "partial" if any(e[2] for e in task_evidence) else "completed"
        if "ERROR.txt" in task_files:
            task_status = "failed"
        elif returncode:
            task_status = "incomplete (run failed)"
        child = root / name
        write_index(child, task_files, task_retained, task_status, task_evidence, run_id, skill, provenance="../../_provenance/", project=config.get("project", {}).get("id"))
        task_entries.append((name, task_status))
    root_files = {name for name in files if not name.startswith("comparisons/")}
    write_index(root, root_files, retained, status, evidence, run_id, skill, task_entries, len(old), project=config.get("project", {}).get("id"), legacy_retained=legacy_retained)
    registry = {"schema_version": 1, "skill": skill, "run_id": run_id, "status": status,
                "exit_status": returncode, "action": action,
                "entrypoints": [INDEX] + [name + "/" + INDEX for name, _ in task_entries], "current_files": sorted(current),
                "retained_stage_files": sorted(retained), "legacy_unverified_files": sorted(legacy_retained & retained), "unclaimed_previous_files": sorted(old)}
    (technical / REGISTRY).write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    return registry


def supervise(output, before, skill, config, executor_exit_status, record, run_id=None):
    """Preserve the scientific exit status and report any failed delivery separately."""
    record["executor_exit_status"] = executor_exit_status
    try:
        delivery = finish(output, before, skill, config, executor_exit_status, run_id)
        record["delivery_status"] = "completed"
        record["result_status"] = delivery["status"]
        record["result_entrypoint"] = INDEX
        artifacts = record.setdefault("artifacts", [])
        existing = {str(item.get("path")) for item in artifacts if isinstance(item, dict)}
        for name in delivery["entrypoints"]:
            path = Path(output) / name
            if str(path) not in existing:
                artifacts.append({"path": str(path), "bytes": None, "sha256": None, "mutable": True})
        return executor_exit_status
    except Exception as exc:
        record["delivery_status"] = "failed"
        record["delivery_error"] = f"{type(exc).__name__}: {exc}"
        record["status"] = "failed"
        record["exit_status"] = executor_exit_status or 1
        return executor_exit_status or 1
