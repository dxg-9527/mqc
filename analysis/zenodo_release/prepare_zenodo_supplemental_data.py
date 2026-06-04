#!/usr/bin/env python3
"""Prepare supplemental MQC manuscript/source-data archives for Zenodo.

This complements ``zenodo_upload_mqc_models.py``. It packages small source-data
tables, figure-generation scripts, experimental validation inputs, and large MQC
reference resources that are unsuitable for a normal GitHub repository.
"""

from __future__ import annotations

import csv
import hashlib
import os
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("MQC_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
OUT_DIR = PROJECT_ROOT / "zenodo_upload"
PAPER_DIR = PROJECT_ROOT / "others/定量分析/result/论文"


@dataclass(frozen=True)
class PackageSpec:
    filename: str
    description: str
    members: tuple[Path, ...]
    arc_prefix: str


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel_project(path: Path) -> str:
    return path.relative_to(PROJECT_ROOT).as_posix()


def paper_source_files() -> tuple[Path, ...]:
    suffixes = {".csv", ".xlsx", ".md", ".py", ".tsv"}
    return tuple(sorted(p for p in PAPER_DIR.iterdir() if p.is_file() and p.suffix in suffixes))


def existing(paths: list[Path]) -> tuple[Path, ...]:
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing supplemental files: " + ", ".join(p.as_posix() for p in missing))
    return tuple(paths)


def package_specs() -> tuple[PackageSpec, ...]:
    analysis_inputs = existing(
        [
            PROJECT_ROOT / "tmp/bigg_new.xlsx",
            PROJECT_ROOT / "tmp/virtual_new.xlsx",
            PROJECT_ROOT / "tmp/seed_new.xlsx",
            PROJECT_ROOT / "tmp/carveme_new.xlsx",
        ]
    )

    experimental_inputs = existing(
        [
            PROJECT_ROOT / "others/定量分析/定量生长数据文章 (version 1).xlsx",
            PROJECT_ROOT / "others/定量分析/37碳源-29菌生长情况25.1.13.xlsx",
            PROJECT_ROOT / "others/定量分析/37碳源-30菌情况表.xlsx",
            PROJECT_ROOT / "others/定量分析/37碳源-31菌情况表.xlsx",
            PROJECT_ROOT / "others/定量分析/改_37碳源-29菌生长情况25.1.13.xlsx",
            PROJECT_ROOT / "others/定量分析/文献中的不同生长速率下的C13通量分布数据，可以进行拟合msb201352-sup-0009 (version 1).xlsx",
            PROJECT_ROOT / "others/定量分析/new文献中的不同生长速率下的C13通量分布数据，可以进行拟合msb201352-sup-0009 (version 1).xlsx",
            PROJECT_ROOT / "mqc/summary/biolog_bsub.tsv",
            PROJECT_ROOT / "mqc/summary/biolog_ecol.tsv",
            PROJECT_ROOT / "mqc/summary/biolog_paer.tsv",
            PROJECT_ROOT / "mqc/summary/biolog_rsol.tsv",
            PROJECT_ROOT / "mqc/summary/biolog_sone.tsv",
            PROJECT_ROOT / "mqc/summary/essential_MOPS_succ.tsv",
            PROJECT_ROOT / "mqc/summary/essential_glc_mm.csv",
            PROJECT_ROOT / "mqc/summary/essential_rich_medium.tsv",
            PROJECT_ROOT / "mqc/summary/essentiality.tsv",
            PROJECT_ROOT / "mqc/summary/essentiality_price2016.tsv",
            PROJECT_ROOT / "mqc/summary/essentiality_rich_medium.tsv",
            PROJECT_ROOT / "mqc/summary/文献模型与自动构建模型作图版xlsx.xlsx",
        ]
    )

    reference_resources = existing(
        [
            PROJECT_ROOT / "mqc/summary/merged_meta_model.xml",
            PROJECT_ROOT / "mqc/summary/general_library.xml",
            PROJECT_ROOT / "mqc/summary/meta_general_library.xml",
            PROJECT_ROOT / "mqc/summary/general_library.json",
            PROJECT_ROOT / "mqc/summary/modelseed_metabolites.xlsx",
            PROJECT_ROOT / "mqc/summary/BiGGdata_map.xlsx",
            PROJECT_ROOT / "mqc/summary/meta_general_library.xlsx",
            PROJECT_ROOT / "mqc/summary/MetacycRule.xlsx",
            PROJECT_ROOT / "mqc/summary/bigg_metacyc.xlsx",
            PROJECT_ROOT / "mqc/summary/modelseed_met.xlsx",
            PROJECT_ROOT / "mqc/summary/modelseed_reactions.xlsx",
            PROJECT_ROOT / "mqc/summary/kegg_met.xlsx",
            PROJECT_ROOT / "mqc/summary/meta_met.xlsx",
            PROJECT_ROOT / "mqc/summary/virtual_met.xlsx",
            PROJECT_ROOT / "mqc/summary/met_id_shortlist.json",
        ]
    )

    return (
        PackageSpec(
            filename="mqc_manuscript_figure_source_data.zip",
            description="Manuscript figure source tables, supplementary tables, source-data notes, and figure-generation scripts.",
            members=paper_source_files(),
            arc_prefix="manuscript_figure_source_data",
        ),
        PackageSpec(
            filename="mqc_analysis_input_tables.zip",
            description="Primary tmp/*_new.xlsx analysis input tables used for Figure 2 mechanistic clustering and model-list derivation.",
            members=analysis_inputs,
            arc_prefix="analysis_input_tables",
        ),
        PackageSpec(
            filename="mqc_experimental_validation_raw_data.zip",
            description="Raw experimental validation and literature-derived growth/flux/essentiality inputs used in the manuscript analyses.",
            members=experimental_inputs,
            arc_prefix="experimental_validation_raw_data",
        ),
        PackageSpec(
            filename="mqc_reference_resources.zip",
            description="MQC reference models and annotation resources required by model standardization and gap-filling modules.",
            members=reference_resources,
            arc_prefix="mqc_reference_resources",
        ),
    )


def write_zip(spec: PackageSpec, force: bool = False) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = OUT_DIR / spec.filename
    if zip_path.exists() and not force:
        print(f"SKIP existing archive: {zip_path.name}")
        return zip_path
    tmp_path = zip_path.with_suffix(zip_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    print(f"Creating {zip_path.name} with {len(spec.members)} files")
    with zipfile.ZipFile(tmp_path, "w", allowZip64=True) as zf:
        for src in spec.members:
            compress_type = zipfile.ZIP_STORED if src.name.endswith(".gz") else zipfile.ZIP_DEFLATED
            arcname = f"{spec.arc_prefix}/{rel_project(src)}"
            zf.write(src, arcname=arcname, compress_type=compress_type, compresslevel=6)
    tmp_path.replace(zip_path)
    return zip_path


def write_manifest(specs: tuple[PackageSpec, ...], archive_paths: dict[str, Path]) -> tuple[Path, Path, Path]:
    manifest_path = OUT_DIR / "mqc_supplemental_source_data_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "package",
                "description",
                "archive_member",
                "source_relative_path",
                "size_bytes",
                "sha256",
            ],
        )
        writer.writeheader()
        for spec in specs:
            for src in spec.members:
                writer.writerow(
                    {
                        "package": spec.filename,
                        "description": spec.description,
                        "archive_member": f"{spec.arc_prefix}/{rel_project(src)}",
                        "source_relative_path": rel_project(src),
                        "size_bytes": src.stat().st_size,
                        "sha256": sha256_file(src),
                    }
                )

    readme_path = OUT_DIR / "README_MQC_supplemental_source_data.md"
    readme_path.write_text(
        "\n".join(
            [
                "# MQC supplemental source-data archives",
                "",
                f"Prepared on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}.",
                "",
                "These archives complement the paired before/after GEM model archives in the same Zenodo draft.",
                "They contain manuscript source data, figure-generation scripts, analysis input tables,",
                "experimental validation inputs, and MQC reference resources that are unsuitable for a normal GitHub repository.",
                "",
                "## Packages",
                "",
                *[f"- `{spec.filename}`: {spec.description}" for spec in specs],
                "",
                "## Metadata",
                "",
                "- `mqc_supplemental_source_data_manifest.csv`: file-level source path, archive member path, size and SHA-256 checksum.",
                "- `checksums_sha256_supplemental.tsv`: archive-level SHA-256 checksums for supplemental files.",
                "",
                "The external `others/compounds.sqlite` and `others/chem_prop.tsv` resources were not included here because",
                "they were not directly referenced by the current MQC source package during this audit.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    upload_files = [
        *[archive_paths[spec.filename] for spec in specs],
        manifest_path,
        readme_path,
    ]
    checksum_path = OUT_DIR / "checksums_sha256_supplemental.tsv"
    with checksum_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["filename", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        for path in upload_files:
            writer.writerow(
                {
                    "filename": path.name,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    upload_files.append(checksum_path)

    upload_manifest_path = OUT_DIR / "zenodo_supplemental_upload_files.csv"
    with upload_manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["filename", "path", "size_bytes"])
        writer.writeheader()
        for path in upload_files:
            writer.writerow(
                {
                    "filename": path.name,
                    "path": path.as_posix(),
                    "size_bytes": path.stat().st_size,
                }
            )
    return manifest_path, checksum_path, upload_manifest_path


def main() -> int:
    specs = package_specs()
    archive_paths = {spec.filename: write_zip(spec) for spec in specs}
    manifest_path, checksum_path, upload_manifest_path = write_manifest(specs, archive_paths)
    print("Prepared supplemental Zenodo package:")
    print(f"  manifest: {manifest_path}")
    print(f"  checksums: {checksum_path}")
    print(f"  upload_manifest: {upload_manifest_path}")
    for spec in specs:
        archive = archive_paths[spec.filename]
        print(f"  {archive.name}: {len(spec.members)} files, {archive.stat().st_size / 1_000_000:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
