#!/usr/bin/env python3
"""Rebuild or publish the approved condition-wise waterfall Fig5."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parents[1]
FINAL = PACKAGE / "fig5_final"
TABLES = [
    "bacillus_substrate_growth.csv",
    "ecoli_exchange_fluxes.csv",
    "ecoli_substrate_growth.csv",
    "condition_values_and_changes.csv",
    "waterfall_ranked_coordinates.csv",
    "mae_summary.csv",
]
DBS = ["CarveMe", "BiGG", "ModelSEED", "VMH"]


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
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def publish():
    panel_dir = PACKAGE / "panels/Fig5"
    full = PACKAGE / "full_figures"
    data = PACKAGE / "source_data/fig5"
    for folder in [panel_dir, full, data]:
        folder.mkdir(parents=True, exist_ok=True)
    for group in "abc":
        for db in DBS:
            for ext in ["png", "pdf", "svg"]:
                name = f"Fig5{group}_{db}.{ext}"
                shutil.copy2(FINAL / "panels" / name, panel_dir / name)
    for ext in ["png", "pdf", "svg"]:
        shutil.copy2(
            FINAL / f"Fig5_final_waterfall.{ext}", full / f"Fig5_combined.{ext}"
        )
    shutil.copy2(FINAL / "figure_legend.txt", full / "Fig5_figure_legend.txt")
    rows = []
    for name in TABLES:
        origin = FINAL / "source_data" / name
        target = data / name
        shutil.copy2(origin, target)
        with target.open(encoding="utf-8-sig", newline="") as f:
            count = sum(1 for _ in csv.reader(f)) - 1
        rows.append(
            dict(
                package_file=str(target.relative_to(PACKAGE)),
                upstream_source=str(origin.relative_to(PACKAGE)),
                upstream_element="whole_file",
                rows=count,
                sha256=digest(target),
                purpose="Approved 20261012 Fig5 waterfall: archived paired values and derived changes",
            )
        )
    update_csv(
        PACKAGE / "source_data_manifest.csv",
        rows,
        "package_file",
        lambda value: value.startswith("source_data/fig5/"),
    )
    pdf = full / "Fig5_combined.pdf"
    png = full / "Fig5_combined.png"
    row = dict(
        figure="Fig5",
        component_panel_pdfs=";".join(
            f"panels/Fig5/Fig5{g}_{db}.pdf" for g in "abc" for db in DBS
        ),
        combined_pdf=str(pdf.relative_to(PACKAGE)),
        combined_pdf_sha256=digest(pdf),
        combined_png=str(png.relative_to(PACKAGE)),
        combined_png_sha256=digest(png),
        png_dpi=600,
    )
    update_csv(
        full / "full_figure_manifest.csv",
        [row],
        "figure",
        lambda value: value == "Fig5",
    )
    import build_panel_manifest

    build_panel_manifest.main()
    print(
        "Fig5: published 12 waterfall panels, full figure, legend and six source tables."
    )


def main():
    subprocess.run(
        [sys.executable, str(FINAL / "code/draw_fig5_waterfall.py")], check=True
    )
    subprocess.run([sys.executable, str(FINAL / "code/validate_final.py")], check=True)
    publish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--publish-only", action="store_true")
    args = parser.parse_args()
    publish() if args.publish_only else main()
