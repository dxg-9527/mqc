"""Build editable artwork from archived inputs. No metabolic simulation."""

import json
from drawing import Drawing, ROOT, E
from panels_ab import frame, panel_a, panel_b1, panel_b2
from panels_lower import panel_b3, panel_b4, panel_c, legend
from embed_fonts import embed

STEM = "Fig4_final"


def standalone_panels():
    """Rebuild individual panels from drawing functions, never crop a bitmap."""
    folder = ROOT / "panels"
    folder.mkdir(exist_ok=True)
    items = [
        ("a", panel_a, "No-carbon phenotype screen"),
        ("b1", panel_b1, "FDR-dependent redox cycle"),
        ("b2", panel_b2, "ABAH-dependent redox cycle"),
        ("b3", panel_b3, "Proton-supported ATP synthesis"),
        ("b4", panel_b4, "Unsupported metabolite production"),
    ]
    for code, fn, title in items:
        d = Drawing()
        fn(d)
        shapes = list(d.slide.shapes)
        x0 = min(s.left for s in shapes) / 914400
        y0 = min(s.top for s in shapes) / 914400
        x1 = max(s.left + s.width for s in shapes) / 914400
        y1 = max(s.top + s.height for s in shapes) / 914400
        d.text(d.slide.shapes, x0, y0 - 0.62, x1 - x0, 0.50, title, 16, bold=True)
        pad = 0.18
        dx = pad - x0
        dy = pad - (y0 - 0.62)
        for shape in d.slide.shapes:
            shape.left += E(dx)
            shape.top += E(dy)
        d.prs.slide_width = E(x1 - x0 + 2 * pad)
        d.prs.slide_height = E(y1 - y0 + 0.62 + 2 * pad)
        out = folder / f"Fig4_{code}.pptx"
        d.prs.save(out)
        embed(out)
    return folder


def main():
    d = Drawing()
    for fn in [frame, panel_a, panel_b1, panel_b2, panel_b3, panel_b4, panel_c, legend]:
        fn(d)
    d.prs.core_properties.title = "MQC diagnosis and targeted refinement of iCN718"
    d.slide.notes_slide.notes_text_frame.text = "Refined from the accepted 6.1 layout. Panel c is generated from archived CSVs. No model optimization was run. See figure_legend.md and README.md for conditions, units, historical imbalance evidence, independent objective tests and provenance."
    out = ROOT / (STEM + ".pptx")
    d.prs.save(out)
    embed(out)
    (ROOT / "code" / "text_inventory.json").write_text(
        json.dumps(d.texts, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    issues = d.check_text()
    (ROOT / "code" / "text_fit_check.json").write_text(
        json.dumps(issues, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    if issues:
        raise RuntimeError(
            "Text dimensions failed: " + json.dumps(issues, ensure_ascii=False)
        )
    standalone_panels()
    print(json.dumps(issues, ensure_ascii=False, indent=2))
    print(out)


if __name__ == "__main__":
    main()
