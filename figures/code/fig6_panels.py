#!/usr/bin/env python3
"""Rebuild or publish approved three-strain Fig6 from packaged inputs."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parents[1]
FINAL = PACKAGE / "fig6_final"
TABLES = [
    "three_strains_nine_substrates.csv",
    "three_strain_summary.csv",
    "three_strain_change_counts.csv",
    "experimental_fit_flags.csv",
    "MAG_utilization_52_substrates.csv",
    "panel_a_top12_substrates.csv",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update_csv(path, new_rows, key, excluded):
    rows = []
    if path.exists():
        with path.open(encoding="utf-8-sig", newline="") as f:
            rows = [r for r in csv.DictReader(f) if not excluded(r[key])]
    rows += new_rows
    rows.sort(key=lambda r: r[key])
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def publish():
    panel_dir = PACKAGE / "panels/Fig6"
    panel_dir.mkdir(parents=True, exist_ok=True)
    for f in (FINAL / "panels").glob("Fig6*.*"):
        shutil.copy2(f, panel_dir / f.name)
    for ext in ["png", "pdf", "svg"]:
        shutil.copy2(
            FINAL / f"Fig6_recommended_three_strains.{ext}",
            PACKAGE / f"full_figures/Fig6_combined.{ext}",
        )
    shutil.copy2(
        FINAL / "figure_legend.txt", PACKAGE / "full_figures/Fig6_figure_legend.txt"
    )
    rows = []
    for name in TABLES:
        origin = FINAL / "source_data" / name
        target = PACKAGE / "source_data/fig6" / name
        shutil.copy2(origin, target)
        with target.open(encoding="utf-8-sig", newline="") as f:
            n = sum(1 for _ in csv.reader(f)) - 1
        rows.append(
            dict(
                package_file=str(target.relative_to(PACKAGE)),
                upstream_source=str(origin.relative_to(PACKAGE)),
                upstream_element="whole_file",
                rows=n,
                sha256=digest(target),
                purpose="Approved 20261012 three-strain Fig6",
            )
        )
    update_csv(
        PACKAGE / "source_data_manifest.csv",
        rows,
        "package_file",
        lambda value: value.startswith("source_data/fig6/"),
    )
    pdf = PACKAGE / "full_figures/Fig6_combined.pdf"
    png = pdf.with_suffix(".png")
    row = dict(
        figure="Fig6",
        component_panel_pdfs=";".join(
            str(f.relative_to(PACKAGE)) for f in sorted(panel_dir.glob("*.pdf"))
        ),
        combined_pdf=str(pdf.relative_to(PACKAGE)),
        combined_pdf_sha256=digest(pdf),
        combined_png=str(png.relative_to(PACKAGE)),
        combined_png_sha256=digest(png),
        png_dpi=600,
    )
    update_csv(
        PACKAGE / "full_figures/full_figure_manifest.csv",
        [row],
        "figure",
        lambda value: value == "Fig6",
    )
    print("Fig6: published five approved panels, full figure and six source tables.")


def main():
    subprocess.run(
        [sys.executable, str(FINAL / "code/draw_fig6_three_strains.py")], check=True
    )
    publish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--publish-only", action="store_true")
    args = parser.parse_args()
    publish() if args.publish_only else main()
