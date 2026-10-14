"""Render PPT, replace only the data-panel bitmap with its original vector PDF,
and export the assembled figure at 183 mm width. Does not crop raster panels.
"""

from pathlib import Path
import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
import fitz

ROOT = Path(__file__).resolve().parent.parent
STEM = "Fig4_final"
AUDIT = ROOT / "review"
RAW = AUDIT / "office_pdf"


def render_office():
    RAW.mkdir(parents=True, exist_ok=True)
    office = shutil.which("soffice")
    if not office:
        raise RuntimeError(
            "LibreOffice (soffice) is required to render editable PPT files."
        )
    inputs = [ROOT / (STEM + ".pptx")] + sorted((ROOT / "panels").glob("*.pptx"))
    with tempfile.TemporaryDirectory(prefix="mqc_fig4_office_") as tmp:
        profile = (Path(tmp) / "profile").as_uri()
        subprocess.run(
            [
                office,
                "-env:UserInstallation=" + profile,
                "--headless",
                "--convert-to",
                "pdf",
                "--outdir",
                str(RAW),
                *map(str, inputs),
            ],
            check=True,
            timeout=180,
        )
    for src in inputs:
        assert (RAW / (src.stem + ".pdf")).exists(), src


def png(doc, path, dpi=600):
    pix = doc[0].get_pixmap(dpi=dpi, alpha=False)
    pix.save(path)
    return [pix.width, pix.height]


def main(skip_office=False):
    if not skip_office:
        render_office()
    doc = fitz.open(RAW / (STEM + ".pdf"))
    page = doc[0]
    images = page.get_images(full=True)
    assert len(images) == 1, "Expected only the Panel c PNG in the PPT render"
    im = images[0]
    rects = page.get_image_rects(im[0])
    assert len(rects) == 1
    rect = rects[0]
    assert abs(rect.x0 - 12.65 * 72) < 0.1 and abs(rect.y0 - 0.41 * 72) < 0.1
    assert abs(rect.height - 9.18 * 72) < 0.1
    # Remove the raster drawing command, keeping all native text and paths.
    pattern = rb"/" + re.escape(im[7].encode()) + rb"\s+Do\b"
    removed = 0
    for xref in page.get_contents():
        stream = doc.xref_stream(xref)
        replacement, n = re.subn(pattern, b"", stream)
        if n:
            doc.update_stream(xref, replacement)
            removed += n
    assert removed == 1
    panel = fitz.open(
        ROOT / "Fig4_panel_c" / "output" / "Fig4_panel_c_post_QC_validation.pdf"
    )
    page.show_pdf_page(rect, panel, 0, overlay=True)
    page.clean_contents(sanitize=True)
    vector_source = AUDIT / "Fig4_design_size_vector.pdf"
    doc.save(vector_source, garbage=4, deflate=True)
    doc.close()
    panel.close()
    vector = fitz.open(vector_source)
    assert (
        vector[0].get_image_info() == []
    ), "Final plot must contain no rendered bitmap objects"
    final = fitz.open()
    width = 183 / 25.4 * 72
    height = width * vector[0].rect.height / vector[0].rect.width
    final.new_page(width=width, height=height).show_pdf_page(
        fitz.Rect(0, 0, width, height), vector, 0
    )
    out = ROOT / (STEM + ".pdf")
    final.save(out, garbage=4, deflate=True)
    size = png(final, ROOT / (STEM + "_600dpi.png"))
    final[0].get_pixmap(
        matrix=fitz.Matrix(2400 / width, 2400 / width), alpha=False
    ).save(ROOT / "Fig4_preview.png")
    gray = final[0].get_pixmap(dpi=300, colorspace=fitz.csGRAY, alpha=False)
    gray.save(AUDIT / "Fig4_183mm_grayscale.png")
    shutil.copy2(out, AUDIT / "Fig4_183mm_double_column.pdf")
    final.close()
    vector.close()
    for raw in sorted(RAW.glob("Fig4_*.pdf")):
        target = ROOT / "panels" / raw.name
        shutil.copy2(raw, target)
        with fitz.open(target) as single:
            png(single, target.with_suffix(".png"), 300)
    for ext in ["pdf", "svg", "png"]:
        shutil.copy2(
            ROOT
            / "Fig4_panel_c"
            / "output"
            / ("Fig4_panel_c_post_QC_validation." + ext),
            ROOT / "panels" / ("Fig4_c." + ext),
        )
    with fitz.open(out) as check:
        fonts = sorted({f[3] for f in check.get_page_fonts(0, full=True)})
        spans = [
            s
            for b in check[0].get_text("dict")["blocks"]
            if "lines" in b
            for l in b["lines"]
            for s in l["spans"]
            if s["text"].strip()
        ]
        assert all("LiberationSans" in s["font"].replace(" ", "") for s in spans)
        assert not check[0].get_image_info()
        chars = check[0].get_text()
        assert (
            "Boundary normalization" not in chars
            and "Functional re-testing" not in chars
        )
        assert "Reducing equivalents" in chars and "NADHs" not in chars
        for required in [
            "665.802",
            "520.556",
            "18.333",
            "14.167",
            "20.833",
            "0.169",
            "0.252",
        ]:
            assert required in chars, required
        checks = {
            "final_pdf_width_mm": round(check[0].rect.width / 72 * 25.4, 3),
            "final_pdf_height_mm": round(check[0].rect.height / 72 * 25.4, 3),
            "final_png_dpi": 600,
            "final_png_pixels": size,
            "rendered_bitmap_objects": len(check[0].get_image_info()),
            "pdf_fonts": fonts,
            "min_base_font_pt_at_183mm": 13.5 * 183 / 406.4,
            "min_span_pt_including_superscripts": min(s["size"] for s in spans),
            "office_rendered_panels": len(list(RAW.glob("*.pdf"))),
            "scientific_data_panel": "Original Matplotlib vector PDF; no raster crop or manual point movement",
        }
    (AUDIT / "export_check.json").write_text(
        json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(checks, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-office", action="store_true")
    main(parser.parse_args().skip_office)
