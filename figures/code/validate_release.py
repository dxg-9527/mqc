"""Validate packaged figure data, output completeness and recorded hashes."""

from pathlib import Path
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    data = read(ROOT / "source_data_manifest.csv")
    for row in data:
        path = ROOT / row["package_file"]
        assert digest(path) == row["sha256"], path
        assert len(read(path)) == int(row["rows"]), path
    panels = read(ROOT / "panel_code_data_manifest.csv")
    assert len(panels) == 41
    assert len({r["standalone_panel"] for r in panels}) == 41
    for row in panels:
        assert (ROOT / row["code_file"]).is_file()
        for name in row["source_data_files"].split(";"):
            assert (ROOT / name).is_file(), name
        for field in ["png_file", "pdf_file", "svg_file"]:
            assert (ROOT / row[field]).is_file(), row[field]
        assert digest(ROOT / row["png_file"]) == row["png_sha256"]
    figures = read(ROOT / "full_figures/full_figure_manifest.csv")
    assert {r["figure"] for r in figures} == {"Fig2", "Fig3", "Fig4", "Fig5", "Fig6"}
    for row in figures:
        for ext in ["png", "pdf"]:
            assert (
                digest(ROOT / row["combined_" + ext])
                == row["combined_" + ext + "_sha256"]
            )
    assert (ROOT / "full_figures/Fig1.png").is_file()
    import fig2_module_burden

    burden = fig2_module_burden.derive()
    assert len(burden) == 24
    result = {
        "status": "passed",
        "source_tables": len(data),
        "standalone_panels": len(panels),
        "complete_figures": len(figures) + 1,
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
