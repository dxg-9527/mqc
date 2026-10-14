#!/usr/bin/env python3
"""Build the upload-oriented panel/code/data mapping manifest."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
PANELS = PACKAGE / "panels"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_files(stem: str) -> list[str]:
    if stem.startswith("Fig4_"):
        if stem == "Fig4_a":
            return ["source_data/fig4/diagnostic_quantitative_data.csv"]
        if stem == "Fig4_c":
            return [
                "source_data/fig4/no_carbon_objectives.csv",
                "source_data/fig4/matched_glucose_retest.csv",
            ]
        return [
            "source_data/fig4/causal_path_fluxes.csv",
            "source_data/fig4/diagnostic_quantitative_data.csv",
            "source_data/fig4/reaction_revisions.csv",
            "source_data/fig4/reaction_stoichiometry.csv",
        ]
    if stem.startswith("Fig2a"):
        return ["source_data/fig2/overall_error_prevalence.csv"]
    if stem.startswith("Fig2b"):
        return ["source_data/fig2/module_error_prevalence.csv"]
    if stem.startswith("Fig2c"):
        return ["source_data/fig2/biomass_issue_prevalence.csv"]
    if stem.startswith("Fig2d"):
        return [
            "source_data/fig2/core_model_qc_flags.csv",
            "source_data/fig2/module_burden_distribution.csv",
        ]
    if stem.startswith("Fig3a"):
        return ["source_data/fig3/atp_profiles.csv"]
    if stem.startswith("Fig3b"):
        return ["source_data/fig3/nadh_profiles.csv"]
    if stem.startswith("Fig3c") or stem.startswith("Fig3d_final"):
        return ["source_data/fig3/growth_distributions.csv"]
    if stem.startswith("Fig3d_CarveMe"):
        return ["source_data/fig3/carveme_growth_ecdf.csv"]
    if stem.startswith("Fig3d_VMH"):
        return ["source_data/fig3/vmh_growth_ecdf.csv"]
    if stem.startswith("Fig5a"):
        return [
            "source_data/fig5/bacillus_substrate_growth.csv",
            "source_data/fig5/waterfall_ranked_coordinates.csv",
        ]
    if stem.startswith("Fig5b"):
        return [
            "source_data/fig5/ecoli_exchange_fluxes.csv",
            "source_data/fig5/waterfall_ranked_coordinates.csv",
        ]
    if stem.startswith("Fig5c"):
        return [
            "source_data/fig5/ecoli_substrate_growth.csv",
            "source_data/fig5/waterfall_ranked_coordinates.csv",
        ]
    if stem.startswith("Fig6a"):
        return ["source_data/fig6/MAG_utilization_52_substrates.csv"]
    if stem.startswith(("Fig6b", "Fig6c", "Fig6d", "Fig6e")):
        return ["source_data/fig6/three_strains_nine_substrates.csv"]
    if stem.startswith("Fig4a"):
        return ["source_data/fig4/matched_glucose_comparison.csv"]
    if stem.startswith("Fig4b"):
        return [
            "source_data/fig4/diagnostic_quantitative_data.csv",
            "source_data/fig4/reaction_revisions.csv",
            "source_data/fig4/reaction_stoichiometry.csv",
        ]
    if stem.startswith("Fig4c"):
        return [
            "source_data/fig4/causal_path_fluxes.csv",
            "source_data/fig4/reducing_power_paths.csv",
            "source_data/fig4/reaction_stoichiometry.csv",
        ]
    raise ValueError(f"No source mapping for {stem}")


def main() -> None:
    pngs = sorted(PANELS.glob("Fig*/*.png"))
    rows = []
    for png in pngs:
        stem = png.stem
        figure = png.parent.name
        pdf = png.with_suffix(".pdf")
        svg = png.with_suffix(".svg")
        if not pdf.exists() or not svg.exists():
            raise FileNotFoundError(f"Missing PDF/SVG sibling for {png}")
        code = f"code/{figure.lower()}_panels.py"
        sources = source_files(stem)
        for relative in [code, *sources]:
            if not (PACKAGE / relative).exists():
                raise FileNotFoundError(relative)
        rows.append(
            {
                "manuscript_figure": figure,
                "standalone_panel": stem,
                "code_file": code,
                "source_data_files": ";".join(sources),
                "png_file": str(png.relative_to(PACKAGE)),
                "pdf_file": str(pdf.relative_to(PACKAGE)),
                "svg_file": str(svg.relative_to(PACKAGE)),
                "png_sha256": sha256(png),
            }
        )
    if len(rows) != 41:
        raise ValueError(f"Expected 41 standalone panels, found {len(rows)}")
    output = PACKAGE / "panel_code_data_manifest.csv"
    with output.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} panel mappings to {output}")


if __name__ == "__main__":
    main()
