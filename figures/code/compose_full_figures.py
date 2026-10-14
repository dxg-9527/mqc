#!/usr/bin/env python3
"""Compose manuscript-level figures from the code-generated panel PDFs.

The standalone panels intentionally omit the a/b/c labels.  This script adds
those labels only to the combined figures and keeps all panel artwork vector in
the PDF output.  A 600 dpi PNG is rendered from each combined PDF for convenient
preview and manuscript insertion.
"""

from __future__ import annotations

import csv
import hashlib
import tempfile
import shutil
from dataclasses import dataclass
from pathlib import Path

import fitz
import matplotlib.pyplot as plt
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject
from panel_style import apply_style


PACKAGE = Path(__file__).resolve().parents[1]
PANELS = PACKAGE / "panels"
OUTPUT = PACKAGE / "full_figures"
POINTS_PER_INCH = 72.0
PNG_DPI = 600


@dataclass(frozen=True)
class Placement:
    panel: str
    box: tuple[float, float, float, float]
    """Panel path and normalized (left, top, width, height) placement box."""


@dataclass(frozen=True)
class Label:
    text: str
    x: float
    y: float
    """Label and normalized top-left page coordinates."""


@dataclass(frozen=True)
class FigureLayout:
    name: str
    size_inches: tuple[float, float]
    placements: tuple[Placement, ...]
    labels: tuple[Label, ...]


def grid_boxes(
    *,
    rows: int,
    columns: int,
    left: float,
    top: float,
    right: float,
    bottom: float,
    horizontal_gap: float,
    vertical_gap: float,
) -> list[tuple[float, float, float, float]]:
    width = (right - left - horizontal_gap * (columns - 1)) / columns
    height = (bottom - top - vertical_gap * (rows - 1)) / rows
    return [
        (
            left + column * (width + horizontal_gap),
            top + row * (height + vertical_gap),
            width,
            height,
        )
        for row in range(rows)
        for column in range(columns)
    ]


def layouts() -> tuple[FigureLayout, ...]:
    # Fig2 presents the cohort-level overview; Fig3 groups quantitative outcomes.
    fig2_boxes = grid_boxes(
        rows=2,
        columns=2,
        left=0.040,
        top=0.028,
        right=0.985,
        bottom=0.985,
        horizontal_gap=0.042,
        vertical_gap=0.055,
    )
    fig2_files = (
        "Fig2/Fig2a_overall_error_prevalence.pdf",
        "Fig2/Fig2b_module_error_prevalence.pdf",
        "Fig2/Fig2c_biomass_issue_prevalence.pdf",
        "Fig2/Fig2d_module_burden_distribution.pdf",
    )
    fig3_files = (
        *(f"Fig3/Fig3a_ATP_{db}.pdf" for db in ("CarveMe", "BiGG", "ModelSEED", "VMH")),
        *(
            f"Fig3/Fig3b_NADH_{db}.pdf"
            for db in ("CarveMe", "BiGG", "ModelSEED", "VMH")
        ),
        *(
            f"Fig3/Fig3c_{stage}_growth_distribution.pdf"
            for stage in ("Initial", "Middle", "Final")
        ),
        "Fig3/Fig3d_CarveMe_growth_ECDF.pdf",
        "Fig3/Fig3d_VMH_growth_ECDF.pdf",
        "Fig3/Fig3d_final_range_summary.pdf",
    )
    fig3_boxes = []
    row_tops = [0.025, 0.267, 0.509, 0.751]
    for top, columns in zip(row_tops, (4, 4, 3, 3)):
        fig3_boxes.extend(
            grid_boxes(
                rows=1,
                columns=columns,
                left=0.030,
                top=top,
                right=0.992,
                bottom=top + 0.227,
                horizontal_gap=0.009,
                vertical_gap=0,
            )
        )
    return (
        FigureLayout(
            "Fig2",
            (13.6, 8.6),
            tuple(Placement(path, box) for path, box in zip(fig2_files, fig2_boxes)),
            tuple(
                Label(letter, box[0] - 0.026, box[1] - 0.012)
                for letter, box in zip("abcd", fig2_boxes)
            ),
        ),
        FigureLayout(
            "Fig3",
            (15.8, 14.8),
            tuple(Placement(path, box) for path, box in zip(fig3_files, fig3_boxes)),
            tuple(
                Label(letter, 0.007, top - 0.010)
                for letter, top in zip("abcd", row_tops)
            ),
        ),
    )


def add_panel(
    page: PageObject,
    panel_path: Path,
    box: tuple[float, float, float, float],
    width: float,
    height: float,
) -> None:
    panel = PdfReader(str(panel_path)).pages[0]
    panel_width = float(panel.mediabox.width)
    panel_height = float(panel.mediabox.height)
    left, top, box_width, box_height = box
    target_width = box_width * width
    target_height = box_height * height
    scale = min(target_width / panel_width, target_height / panel_height)
    x = left * width + (target_width - panel_width * scale) / 2
    y = (
        height
        - (top * height + target_height)
        + (target_height - panel_height * scale) / 2
    )
    transformation = Transformation().scale(scale).translate(x, y)
    page.merge_transformed_page(panel, transformation, over=True)


def label_overlay(layout: FigureLayout, path: Path) -> None:
    # Use the same font embedding in standalone and generate_all invocation.
    apply_style()
    width, height = layout.size_inches
    figure = plt.figure(figsize=(width, height), facecolor="none")
    for label in layout.labels:
        figure.text(
            label.x,
            1 - label.y,
            label.text,
            ha="left",
            va="top",
            fontsize=16,
            fontfamily="Liberation Sans",
            fontweight="bold",
            color="black",
        )
    figure.savefig(path, format="pdf", transparent=True)
    plt.close(figure)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compose(layout: FigureLayout) -> tuple[Path, Path]:
    width = layout.size_inches[0] * POINTS_PER_INCH
    height = layout.size_inches[1] * POINTS_PER_INCH
    page = PageObject.create_blank_page(width=width, height=height)
    for placement in layout.placements:
        add_panel(page, PANELS / placement.panel, placement.box, width, height)

    with tempfile.TemporaryDirectory(prefix="mqc_full_figure_") as temp_directory:
        overlay_path = Path(temp_directory) / "labels.pdf"
        label_overlay(layout, overlay_path)
        overlay = PdfReader(str(overlay_path)).pages[0]
        page.merge_page(overlay, over=True)

    pdf_path = OUTPUT / f"{layout.name}_combined.pdf"
    writer = PdfWriter()
    writer.add_page(page)
    with pdf_path.open("wb") as handle:
        writer.write(handle)

    png_path = OUTPUT / f"{layout.name}_combined.png"
    document = fitz.open(pdf_path)
    try:
        pixmap = document[0].get_pixmap(
            matrix=fitz.Matrix(PNG_DPI / 72, PNG_DPI / 72), alpha=False
        )
        pixmap.set_dpi(PNG_DPI, PNG_DPI)
        pixmap.save(png_path)
        (OUTPUT / f"{layout.name}_combined.svg").write_text(
            document[0].get_svg_image(text_as_path=True), encoding="utf-8"
        )
    finally:
        document.close()
    return pdf_path, png_path


def write_manifest(
    layout_list: tuple[FigureLayout, ...], outputs: dict[str, tuple[Path, Path]]
) -> None:
    manifest = OUTPUT / "full_figure_manifest.csv"
    with manifest.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "figure",
                "component_panel_pdfs",
                "combined_pdf",
                "combined_pdf_sha256",
                "combined_png",
                "combined_png_sha256",
                "png_dpi",
            ),
        )
        writer.writeheader()
        for layout in layout_list:
            pdf_path, png_path = outputs[layout.name]
            writer.writerow(
                {
                    "figure": layout.name,
                    "component_panel_pdfs": ";".join(
                        f"panels/{placement.panel}" for placement in layout.placements
                    ),
                    "combined_pdf": str(pdf_path.relative_to(PACKAGE)),
                    "combined_pdf_sha256": sha256(pdf_path),
                    "combined_png": str(png_path.relative_to(PACKAGE)),
                    "combined_png_sha256": sha256(png_path),
                    "png_dpi": PNG_DPI,
                }
            )


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    # Fig4–6 use the author-approved self-contained figure layouts.
    layout_list = tuple(
        layout for layout in layouts() if layout.name not in ("Fig5", "Fig6", "Fig4")
    )
    outputs = {layout.name: compose(layout) for layout in layout_list}
    write_manifest(layout_list, outputs)
    import fig4_panels

    fig4_panels.publish()
    import fig5_panels

    fig5_panels.publish()
    import fig6_panels

    fig6_panels.publish()
    print(
        f"Complete: wrote {len(layout_list)+3} combined figures, including approved Fig4–6"
    )


if __name__ == "__main__":
    main()
