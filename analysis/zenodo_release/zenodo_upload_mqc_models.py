#!/usr/bin/env python3
"""Prepare and upload paired MQC before/after GEM model archives to Zenodo.

The script creates database-level ZIP archives, public manifests, and optional
Zenodo draft uploads. It intentionally uploads to a draft deposition only; it
does not publish the record.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import zipfile
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.parse import quote

import requests


PROJECT_ROOT = Path(os.environ.get("MQC_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
LOCAL_TEST_DATA = PROJECT_ROOT / "mqc/local_test_data"
TMP_DIR = PROJECT_ROOT / "tmp"
DEFAULT_OUT_DIR = PROJECT_ROOT / "zenodo_upload"
ZENODO_API = "https://zenodo.org/api"


@dataclass(frozen=True)
class ModelPair:
    database: str
    model_id: str
    before_path: Path
    after_path: Path
    before_member: str
    after_member: str
    source_note: str = ""


@dataclass(frozen=True)
class ManifestRow:
    database: str
    status: str
    model_id: str
    paired_model_filename: str
    archive_name: str
    archive_member: str
    source_relative_path: str
    size_bytes: int
    sha256: str
    source_note: str


def rel_project(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def modelseed_candidates(stem: str) -> set[str]:
    """Return plausible normalized ModelSEED identifiers.

    The original files contain mixed suffixes and occasional typos, e.g.
    ``.1_modelseed``, ``.1_modelfeed``, ``.1_moelseed`` and ``.1_Model``.
    The QC outputs usually remove the version and many of these suffixes.
    """
    candidates = {stem}
    no_version = re.sub(r"\.\d+", "", stem)
    candidates.add(no_version)
    suffix_patterns = [
        r"_?modelseed$",
        r"_?modelfeed$",
        r"_?moelseed$",
        r"_?modeseed$",
        r"_?model$",
    ]
    for candidate in list(candidates):
        lowered = candidate
        for pattern in suffix_patterns:
            lowered = re.sub(pattern, "", lowered, flags=re.IGNORECASE)
        candidates.add(lowered)
    return {c for c in candidates if c}


def pair_direct_xml(database: str, before_dir: Path, after_dir: Path) -> list[ModelPair]:
    before_by_name = {p.name: p for p in before_dir.glob("*.xml")}
    pairs: list[ModelPair] = []
    missing: list[str] = []
    for after_path in sorted(after_dir.glob("*.xml")):
        before_path = before_by_name.get(after_path.name)
        if before_path is None:
            missing.append(after_path.name)
            continue
        pairs.append(
            ModelPair(
                database=database,
                model_id=after_path.stem,
                before_path=before_path,
                after_path=after_path,
                before_member=f"before/{database}/{before_path.name}",
                after_member=f"after/{database}/{after_path.name}",
            )
        )
    if missing:
        raise RuntimeError(f"{database}: missing {len(missing)} before files, e.g. {missing[:5]}")
    return pairs


def pair_modelseed() -> list[ModelPair]:
    before_dir = LOCAL_TEST_DATA / "GCF_modelseed"
    after_dir = TMP_DIR / "seed_model"
    before_by_name = {p.name: p for p in before_dir.glob("*.xml")}
    candidate_map: dict[str, list[Path]] = {}
    for before_path in before_dir.glob("*.xml"):
        for key in modelseed_candidates(before_path.stem):
            candidate_map.setdefault(key, []).append(before_path)

    pairs: list[ModelPair] = []
    missing: list[str] = []
    ambiguous: list[tuple[str, list[str]]] = []
    for after_path in sorted(after_dir.glob("*.xml")):
        before_path = before_by_name.get(after_path.name)
        if before_path is None:
            matches: list[Path] = []
            for key in modelseed_candidates(after_path.stem):
                matches.extend(candidate_map.get(key, []))
            unique_matches = sorted(set(matches), key=lambda p: p.name)
            if len(unique_matches) == 1:
                before_path = unique_matches[0]
            elif len(unique_matches) > 1:
                ambiguous.append((after_path.name, [p.name for p in unique_matches]))
            else:
                missing.append(after_path.name)
        if before_path is not None:
            pairs.append(
                ModelPair(
                    database="modelseed",
                    model_id=after_path.stem,
                    before_path=before_path,
                    after_path=after_path,
                    before_member=f"before/modelseed/{before_path.name}",
                    after_member=f"after/modelseed/{after_path.name}",
                    source_note="ModelSEED before filename normalized from version/suffix where needed.",
                )
            )
    if missing or ambiguous:
        raise RuntimeError(
            "modelseed pairing failed: "
            f"missing={missing[:10]}, ambiguous={ambiguous[:3]}"
        )
    return pairs


def pair_carveme() -> tuple[list[ModelPair], list[dict[str, str]]]:
    model_list = LOCAL_TEST_DATA / "embl_gems/model_list.tsv"
    after_dir = TMP_DIR / "carveme_model"
    rows_by_key: dict[str, list[dict[str, str]]] = {}
    all_rows: list[dict[str, str]] = []
    with model_list.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            file_path = row.get("file_path", "")
            if not file_path:
                continue
            base = Path(file_path).name
            if base.endswith(".xml.gz"):
                key = base[:-7]
            elif base.endswith(".xml"):
                key = base[:-4]
            else:
                key = Path(base).stem
            record = dict(row)
            record["_model_key"] = key
            record["_before_path"] = str(LOCAL_TEST_DATA / "embl_gems" / file_path)
            rows_by_key.setdefault(key, []).append(record)
            all_rows.append(record)

    pairs: list[ModelPair] = []
    missing: list[str] = []
    for after_path in sorted(after_dir.glob("*.xml")):
        if after_path.name.endswith("_xml.xml"):
            key = after_path.name[:-8]
        else:
            key = after_path.stem
        matches = rows_by_key.get(key, [])
        if not matches:
            missing.append(after_path.name)
            continue
        row = matches[0]
        before_path = Path(row["_before_path"])
        if not before_path.exists():
            missing.append(after_path.name)
            continue
        pairs.append(
            ModelPair(
                database="carveme",
                model_id=key,
                before_path=before_path,
                after_path=after_path,
                before_member=f"before/carveme/{row['file_path']}",
                after_member=f"after/carveme/{after_path.name}",
                source_note="CarveMe before file selected from embl_gems/model_list.tsv.",
            )
        )
    if missing:
        raise RuntimeError(f"carveme: missing {len(missing)} before files, e.g. {missing[:5]}")
    paired_keys = {pair.model_id for pair in pairs}
    unpaired_rows = [row for row in all_rows if row["_model_key"] not in paired_keys]
    return pairs, unpaired_rows


def build_pairs() -> tuple[list[ModelPair], list[dict[str, str]]]:
    pairs: list[ModelPair] = []
    pairs.extend(pair_direct_xml("bigg", LOCAL_TEST_DATA / "bigg_data", TMP_DIR / "bigg_model"))
    pairs.extend(pair_direct_xml("vmh", LOCAL_TEST_DATA / "virtual_metabolic_human", TMP_DIR / "virtual_model"))
    pairs.extend(pair_modelseed())
    carveme_pairs, carveme_unpaired = pair_carveme()
    pairs.extend(carveme_pairs)
    return pairs, carveme_unpaired


def iter_archive_plan(pairs: list[ModelPair]) -> dict[str, list[tuple[Path, str]]]:
    plan: dict[str, list[tuple[Path, str]]] = {}
    for pair in pairs:
        before_archive = f"mqc_before_{pair.database}_models.zip"
        after_archive = f"mqc_after_{pair.database}_models.zip"
        plan.setdefault(before_archive, []).append((pair.before_path, pair.before_member))
        plan.setdefault(after_archive, []).append((pair.after_path, pair.after_member))
    return plan


def write_zip(archive_path: Path, members: list[tuple[Path, str]], force: bool) -> None:
    if archive_path.exists() and not force:
        print(f"SKIP existing archive: {archive_path.name}")
        return
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = archive_path.with_suffix(archive_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    print(f"Creating {archive_path.name} with {len(members)} files")
    with zipfile.ZipFile(tmp_path, "w", allowZip64=True) as zf:
        for idx, (src, arcname) in enumerate(members, 1):
            compress_type = zipfile.ZIP_STORED if src.name.endswith(".gz") else zipfile.ZIP_DEFLATED
            zf.write(src, arcname=arcname, compress_type=compress_type, compresslevel=6)
            if idx % 500 == 0 or idx == len(members):
                print(f"  {archive_path.name}: {idx}/{len(members)} files")
    tmp_path.replace(archive_path)


def write_csv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def prepare(out_dir: Path, force: bool) -> None:
    pairs, carveme_unpaired = build_pairs()
    out_dir.mkdir(parents=True, exist_ok=True)
    archive_plan = iter_archive_plan(pairs)

    expected_counts = {
        ("before", "bigg"): 73,
        ("after", "bigg"): 73,
        ("before", "vmh"): 806,
        ("after", "vmh"): 806,
        ("before", "modelseed"): 79,
        ("after", "modelseed"): 79,
        ("before", "carveme"): 5585,
        ("after", "carveme"): 5585,
    }
    observed_counts: dict[tuple[str, str], int] = {}
    for pair in pairs:
        observed_counts[("before", pair.database)] = observed_counts.get(("before", pair.database), 0) + 1
        observed_counts[("after", pair.database)] = observed_counts.get(("after", pair.database), 0) + 1
    if observed_counts != expected_counts:
        raise RuntimeError(f"Unexpected before/after counts: {observed_counts}")

    for archive_name, members in sorted(archive_plan.items()):
        write_zip(out_dir / archive_name, members, force=force)

    manifest_rows: list[ManifestRow] = []
    for pair in pairs:
        before_archive = f"mqc_before_{pair.database}_models.zip"
        after_archive = f"mqc_after_{pair.database}_models.zip"
        before_sha = sha256_file(pair.before_path)
        after_sha = sha256_file(pair.after_path)
        manifest_rows.append(
            ManifestRow(
                database=pair.database,
                status="before_mqc",
                model_id=pair.model_id,
                paired_model_filename=pair.after_path.name,
                archive_name=before_archive,
                archive_member=pair.before_member,
                source_relative_path=rel_project(pair.before_path),
                size_bytes=pair.before_path.stat().st_size,
                sha256=before_sha,
                source_note=pair.source_note,
            )
        )
        manifest_rows.append(
            ManifestRow(
                database=pair.database,
                status="after_mqc",
                model_id=pair.model_id,
                paired_model_filename=pair.before_path.name,
                archive_name=after_archive,
                archive_member=pair.after_member,
                source_relative_path=rel_project(pair.after_path),
                size_bytes=pair.after_path.stat().st_size,
                sha256=after_sha,
                source_note="MQC quality-controlled model from tmp database-specific output directory.",
            )
        )

    manifest_path = out_dir / "mqc_before_after_model_manifest.csv"
    write_csv(manifest_path, [asdict(row) for row in manifest_rows], list(asdict(manifest_rows[0]).keys()))

    summary_rows: list[dict[str, object]] = []
    for (status, database), count in sorted(observed_counts.items()):
        total_size = sum(
            row.size_bytes
            for row in manifest_rows
            if row.status == ("before_mqc" if status == "before" else "after_mqc")
            and row.database == database
        )
        summary_rows.append(
            {
                "status": status,
                "database": database,
                "model_count": count,
                "total_size_bytes": total_size,
                "total_size_gb": f"{total_size / 1_000_000_000:.3f}",
            }
        )
    write_csv(
        out_dir / "mqc_before_after_model_summary.csv",
        summary_rows,
        ["status", "database", "model_count", "total_size_bytes", "total_size_gb"],
    )

    if carveme_unpaired:
        fieldnames = [
            "assembly_accession",
            "taxid",
            "organism_name",
            "infraspecific_name",
            "file_path",
            "_model_key",
        ]
        write_csv(out_dir / "carveme_model_list_unpaired_records.csv", carveme_unpaired, fieldnames)

    readme = out_dir / "README_MQC_Zenodo_upload.md"
    readme.write_text(
        "\n".join(
            [
                "# MQC before/after GEM model archives",
                "",
                f"Prepared on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}.",
                "",
                "This draft upload contains paired genome-scale metabolic models before and after MQC quality control.",
                "The files are grouped by model source database and by QC status to keep the Zenodo record below",
                "the 100-file per-record limit.",
                "",
                "## Included archives",
                "",
                "- `mqc_before_bigg_models.zip` and `mqc_after_bigg_models.zip`: 73 paired BiGG models.",
                "- `mqc_before_vmh_models.zip` and `mqc_after_vmh_models.zip`: 806 paired VMH models.",
                "- `mqc_before_modelseed_models.zip` and `mqc_after_modelseed_models.zip`: 79 paired ModelSEED models.",
                "- `mqc_before_carveme_models.zip` and `mqc_after_carveme_models.zip`: 5585 paired CarveMe models.",
                "",
                "## Metadata files",
                "",
                "- `mqc_before_after_model_manifest.csv`: model-level pairing, archive member paths, source paths, sizes and SHA-256 checksums.",
                "- `mqc_before_after_model_summary.csv`: counts and total sizes by source database and QC status.",
                "- `checksums_sha256.tsv`: SHA-256 checksums for all uploaded archives and metadata files.",
                "- `carveme_model_list_unpaired_records.csv`: CarveMe `model_list.tsv` records not paired with QC outputs, if present.",
                "",
                "## Notes",
                "",
                "- CarveMe before-QC files were selected by matching QC output model names to `mqc/local_test_data/embl_gems/model_list.tsv`.",
                "- BiGG and VMH before-QC files were selected by exact filename matching to the MQC output files.",
                "- ModelSEED before-QC files were matched by exact filename where possible and otherwise by normalized assembly identifiers.",
                "- This script creates a Zenodo draft only; publication and DOI release must be done after manual metadata review.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    upload_files = sorted(
        [
            *(out_dir.glob("mqc_before_*_models.zip")),
            *(out_dir.glob("mqc_after_*_models.zip")),
            manifest_path,
            out_dir / "mqc_before_after_model_summary.csv",
            readme,
            *([out_dir / "carveme_model_list_unpaired_records.csv"] if carveme_unpaired else []),
        ],
        key=lambda p: p.name,
    )
    checksum_rows = []
    for path in upload_files:
        checksum_rows.append(
            {
                "filename": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    checksums_path = out_dir / "checksums_sha256.tsv"
    with checksums_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["filename", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checksum_rows)

    upload_files.append(checksums_path)
    write_csv(
        out_dir / "zenodo_upload_files.csv",
        [{"filename": p.name, "path": p.as_posix(), "size_bytes": p.stat().st_size} for p in upload_files],
        ["filename", "path", "size_bytes"],
    )

    print("Prepared upload package:")
    print(f"  output_dir: {out_dir}")
    print(f"  paired_models: {len(pairs)} before + {len(pairs)} after")
    print(f"  upload_file_count: {len(upload_files)}")
    for row in summary_rows:
        print(f"  {row['status']} {row['database']}: {row['model_count']} models, {row['total_size_gb']} GB")


def read_token(token_file: Path) -> str:
    raw = token_file.read_text(encoding="utf-8").strip()
    if "=" in raw and "\n" not in raw:
        raw = raw.split("=", 1)[1].strip()
    token = raw.splitlines()[0].strip()
    if token.lower().startswith("bearer "):
        token = token.split(None, 1)[1].strip()
    if " " in token:
        # Allow token files such as "token value : <TOKEN>" without exposing
        # the token or forcing the user to reformat the file.
        token = max(token.split(), key=len)
    token = token.strip("'\"")
    if not token:
        raise RuntimeError(f"No token found in {token_file}")
    return token


def zenodo_request(method: str, url: str, token: str, **kwargs) -> requests.Response:
    headers = kwargs.pop("headers", {})
    headers.setdefault("Authorization", f"Bearer {token}")
    response = requests.request(method, url, headers=headers, timeout=kwargs.pop("timeout", (30, None)), **kwargs)
    if response.status_code >= 400:
        redacted_url = url.split("?access_token=", 1)[0]
        raise RuntimeError(f"Zenodo {method} failed {response.status_code} for {redacted_url}: {response.text[:1000]}")
    return response


def create_or_load_deposition(out_dir: Path, token: str, new_deposition: bool) -> dict[str, object]:
    state_path = out_dir / "zenodo_deposition.json"
    if state_path.exists() and not new_deposition:
        return json.loads(state_path.read_text(encoding="utf-8"))

    response = zenodo_request("POST", f"{ZENODO_API}/deposit/depositions", token, json={})
    deposition = response.json()
    state_path.write_text(json.dumps(deposition, indent=2, ensure_ascii=False), encoding="utf-8")
    return deposition


def upload(out_dir: Path, token_file: Path, new_deposition: bool, manifest_file: Path | None = None) -> None:
    token = read_token(token_file)
    upload_manifest = manifest_file or out_dir / "zenodo_upload_files.csv"
    if not upload_manifest.exists():
        raise RuntimeError(f"Run prepare first; missing {upload_manifest}")

    deposition = create_or_load_deposition(out_dir, token, new_deposition=new_deposition)
    deposition_id = deposition["id"]
    bucket_url = deposition["links"]["bucket"]
    print(f"Zenodo draft deposition id: {deposition_id}")

    existing_response = zenodo_request("GET", f"{ZENODO_API}/deposit/depositions/{deposition_id}/files", token)
    existing_names = {item["filename"] for item in existing_response.json()}

    with upload_manifest.open(newline="") as handle:
        upload_rows = list(csv.DictReader(handle))

    for row in upload_rows:
        path = Path(row["path"])
        filename = path.name
        if filename in existing_names:
            print(f"SKIP already uploaded: {filename}")
            continue
        print(f"Uploading {filename} ({int(row['size_bytes']) / 1_000_000_000:.3f} GB)")
        with path.open("rb") as handle:
            headers = {
                "Content-Type": "application/octet-stream",
                "Content-Length": str(path.stat().st_size),
            }
            url = f"{bucket_url}/{quote(filename)}"
            # The bucket endpoint accepts token authentication via query params.
            response = requests.put(
                url,
                params={"access_token": token},
                data=handle,
                headers=headers,
                timeout=(30, None),
            )
        if response.status_code >= 400:
            raise RuntimeError(f"Upload failed for {filename}: {response.status_code} {response.text[:1000]}")
        existing_names.add(filename)
        print(f"  uploaded: {filename}")
        time.sleep(1)

    final_response = zenodo_request("GET", f"{ZENODO_API}/deposit/depositions/{deposition_id}", token)
    final_state = final_response.json()
    (out_dir / "zenodo_deposition_after_upload.json").write_text(
        json.dumps(final_state, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    files_response = zenodo_request("GET", f"{ZENODO_API}/deposit/depositions/{deposition_id}/files", token)
    (out_dir / "zenodo_files_after_upload.json").write_text(
        json.dumps(files_response.json(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print("Upload completed.")
    print(f"Draft deposition id: {deposition_id}")
    print(f"Draft URL: https://zenodo.org/deposit/{deposition_id}")
    print(f"Uploaded file count visible via API: {len(files_response.json())}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="Create ZIP archives and manifests.")
    prepare_parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    prepare_parser.add_argument("--force", action="store_true", help="Recreate existing archives.")

    upload_parser = subparsers.add_parser("upload", help="Upload prepared files to a Zenodo draft deposition.")
    upload_parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    upload_parser.add_argument("--token-file", type=Path, default=PROJECT_ROOT / "zenodo.txt")
    upload_parser.add_argument("--manifest-file", type=Path, default=None, help="CSV listing files to upload.")
    upload_parser.add_argument("--new-deposition", action="store_true", help="Ignore saved deposition state and create a new draft.")

    args = parser.parse_args(argv)
    if args.command == "prepare":
        prepare(args.out_dir, force=args.force)
    elif args.command == "upload":
        upload(
            args.out_dir,
            token_file=args.token_file,
            new_deposition=args.new_deposition,
            manifest_file=args.manifest_file,
        )
    else:
        parser.error(f"Unknown command: {args.command}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
