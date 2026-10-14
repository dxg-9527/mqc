#!/usr/bin/env python3
"""Publish/rebuild the author-approved Fig. 4 (20260910).
The self-contained native drawing code and unchanged source data are in fig4_final/.
No bitmap cropping or model optimization is performed.
"""
from pathlib import Path
import csv
import hashlib
import shutil
import subprocess
import sys
import fitz

PACKAGE = Path(__file__).resolve().parents[1]
FINAL = PACKAGE / "fig4_final"
STEM = "Fig4_final"
CODES = ["a", "b1", "b2", "b3", "b4", "c"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def publish():
    dest = PACKAGE / "panels" / "Fig4"
    dest.mkdir(parents=True, exist_ok=True)
    full = PACKAGE / "full_figures"
    full.mkdir(exist_ok=True)
    for code in CODES:
        stem = "Fig4_" + code
        src = FINAL / "panels"
        shutil.copy2(src / (stem + ".pdf"), dest / (stem + ".pdf"))
        if code == "c":
            for ext in ["svg", "png"]:
                shutil.copy2(src / (stem + "." + ext), dest / (stem + "." + ext))
        else:
            shutil.copy2(src / (stem + ".pptx"), dest / (stem + ".pptx"))
            with fitz.open(src / (stem + ".pdf")) as doc:
                doc[0].get_pixmap(dpi=600, alpha=False).save(dest / (stem + ".png"))
                (dest / (stem + ".svg")).write_text(
                    doc[0].get_svg_image(text_as_path=False), encoding="utf-8"
                )
    for source, target in [
        (STEM + ".pdf", "Fig4_combined.pdf"),
        (STEM + "_600dpi.png", "Fig4_combined.png"),
        (STEM + ".pptx", "Fig4_combined.pptx"),
    ]:
        shutil.copy2(FINAL / source, full / target)
    shutil.copy2(FINAL / "figure_legend.md", full / "Fig4_figure_legend.md")
    pairs = [
        (p, PACKAGE / "source_data" / "fig4" / p.name)
        for p in sorted((FINAL / "source_data" / "fig4").glob("*.csv"))
    ]
    pairs += [
        (p, PACKAGE / "source_data" / "fig4" / p.name)
        for p in sorted((FINAL / "Fig4_panel_c" / "data").glob("*.csv"))
    ]
    sources = [
        r
        for r in read(PACKAGE / "source_data_manifest.csv")
        if not r["package_file"].startswith("source_data/fig4/")
    ]
    for src, target in pairs:
        shutil.copy2(src, target)
        sources.append(
            dict(
                package_file=str(target.relative_to(PACKAGE)),
                upstream_source=str(src.relative_to(PACKAGE)),
                upstream_element="author-approved Fig4 20260910; unchanged archived values",
                rows=len(read(target)),
                sha256=sha(target),
                purpose="Fig4 source data",
            )
        )
    write(PACKAGE / "source_data_manifest.csv", sources)
    rows = [r for r in read(full / "full_figure_manifest.csv") if r["figure"] != "Fig4"]
    rows.append(
        dict(
            figure="Fig4",
            component_panel_pdfs=";".join(
                "panels/Fig4/Fig4_" + c + ".pdf" for c in CODES
            ),
            combined_pdf="full_figures/Fig4_combined.pdf",
            combined_pdf_sha256=sha(full / "Fig4_combined.pdf"),
            combined_png="full_figures/Fig4_combined.png",
            combined_png_sha256=sha(full / "Fig4_combined.png"),
            png_dpi=600,
        )
    )
    write(full / "full_figure_manifest.csv", rows)
    import build_panel_manifest

    build_panel_manifest.main()
    print(
        "Fig4 20260910 published: six code-generated panels, approved combined artwork, eight source tables."
    )


def main():
    subprocess.run([sys.executable, str(FINAL / "code" / "rebuild_all.py")], check=True)
    publish()


if __name__ == "__main__":
    if "--publish-only" in sys.argv:
        publish()
    else:
        main()
