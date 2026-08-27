#!/usr/bin/env python3
"""Validate the minimal YAML/TSV contract shared by the scRNA skills."""

import argparse
import csv
import sys
from pathlib import Path


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("samples", type=Path)
    parser.add_argument("--required", nargs="+", default=["sample_id", "condition"])
    args = parser.parse_args()
    if not args.samples.is_file():
        return fail(f"sample table not found: {args.samples}")
    with args.samples.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        missing = [name for name in args.required if name not in fields]
        if missing:
            return fail("missing columns: " + ", ".join(missing))
        rows = list(reader)
    if not rows:
        return fail("sample table contains no samples")
    ids = [row[args.required[0]].strip() for row in rows]
    if any(not value for value in ids):
        return fail(f"blank values in {args.required[0]}")
    if len(ids) != len(set(ids)):
        return fail(f"duplicate values in {args.required[0]}")
    print(f"OK: {len(rows)} samples; columns={','.join(fields)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
