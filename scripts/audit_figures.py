#!/usr/bin/env python3
"""Inventory plot call sites and existing artifacts without changing analysis files."""
import argparse
import csv
import hashlib
import re
from pathlib import Path

PLOT_CALL = re.compile(r"(?:ggsave|savefig|(?:grDevices::)?pdf|(?:grDevices::)?png|save_plot|save_enrichment_plot|DimPlot|FeaturePlot|DotPlot|VlnPlot|heatmap|plotDAbeeswarm|plotNhoodGraphDA)\s*\(")


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def callsites(root):
    rows = []
    for path in sorted(root.rglob("*")):
        if path.suffix not in {".R", ".r", ".py", ".Rmd", ".ipynb"} or "__pycache__" in path.parts:
            continue
        for number, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if PLOT_CALL.search(line):
                rows.append({"source": str(path), "line": number, "code": line.strip()})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = callsites(repo / "toolkit") + callsites(repo / "skills")
    write_tsv(out / "skill_plot_calls.tsv", rows, ["source", "line", "code"])
    project_rows = callsites(args.project / "scripts")
    write_tsv(out / "project_plot_calls.tsv", project_rows, ["source", "line", "code"])
    figures = []
    for path in sorted((args.project / "results").rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".png", ".pdf", ".svg", ".jpg", ".jpeg", ".tiff"}:
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        figures.append({"id": f"F{len(figures)+1:03d}", "path": str(path), "name": path.name,
                        "format": path.suffix[1:], "bytes": path.stat().st_size,
                        "sha256": digest.hexdigest(), "review_status": "inventoried_not_visually_reviewed"})
    write_tsv(out / "existing_figures.tsv", figures, ["id", "path", "name", "format", "bytes", "sha256", "review_status"])
    with (out / "figure_index.md").open("w", encoding="utf-8") as handle:
        handle.write("# Existing figure index\n\nEvery file is inventoried; listing is not visual verification. PDF links refer to complete documents, not just first pages.\n\n")
        for row in figures:
            handle.write(f"- {row['id']} [{Path(row['path']).relative_to(args.project)}](<{row['path']}>)\n")
    # Contact sheets are navigation previews of existing PNGs, not replacement results.
    from PIL import Image, ImageDraw
    pngs = [row for row in figures if row["format"] == "png"]
    for start in range(0, len(pngs), 12):
        sheet = Image.new("RGB", (1500, 1200), "white")
        draw = ImageDraw.Draw(sheet)
        for index, row in enumerate(pngs[start:start+12]):
            x, y = (index % 3) * 500, (index // 3) * 300
            with Image.open(row["path"]) as original:
                preview = original.convert("RGB")
                preview.thumbnail((480, 255))
                sheet.paste(preview, (x + (500-preview.width)//2, y))
            draw.text((x+10, y+260), row["id"] + " " + row["name"][:60], fill="black")
        sheet.save(out / f"contact_sheet_{start//12+1:02d}.png")
    print(f"{len(rows)} skill call sites; {len(project_rows)} project call sites; {len(figures)} figures; {len(pngs)} PNG previews")


if __name__ == "__main__":
    main()
