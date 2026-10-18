"""Standard-library helpers shared by cluster tools (Python 3.8+)."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[1]
STAGES = ('initial_control', 'nadh_control', 'atp_control', 'net_control',
          'yield_control', 'biomass_control')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(v) for v in value]
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=str(path.parent),
                                     prefix=path.name + '.', delete=False) as stream:
        json.dump(clean(value), stream, ensure_ascii=False, indent=2,
                  allow_nan=False, default=str)
        temporary = Path(stream.name)
    os.replace(str(temporary), str(path))


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def under(root, relative):
    root = Path(root).resolve()
    path = (root / relative).resolve()
    path.relative_to(root)  # Reject absolute paths or ../ escapes outside root.
    return path


def verify_release():
    path = REPO / 'cluster/release_lock.json'
    lock = json.loads(path.read_text())
    for relative, digest in {**lock['core_hashes'], **lock['runner_hashes']}.items():
        p = under(REPO, relative)
        if not p.is_file() or sha(p) != digest:
            raise ValueError('Changed or missing release file: ' + relative)
    return sha(path)


def load_inputs(root):
    root = Path(root).resolve()
    meta = json.loads((root / 'bundle.json').read_text())
    if sha(root / 'manifest.csv') != meta['manifest_sha256']:
        raise ValueError('Input manifest checksum mismatch')
    rows = read_csv(root / 'manifest.csv')
    ids = [r['job_id'] for r in rows]
    if len(set(ids)) != len(ids) or len(rows) != meta['model_count']:
        raise ValueError('Duplicate identifiers or incorrect cohort size')
    for r in rows:
        if not r['job_id'].replace('_', '').isalnum():
            raise ValueError('Unsafe job ID')
        under(root, r['input_relative_path'])
    return rows, meta


def git_commit():
    p = subprocess.run(['git', '-C', str(REPO), 'rev-parse', 'HEAD'],
                       capture_output=True, text=True)
    return p.stdout.strip() if p.returncode == 0 else 'unavailable'


def valid_completed(attempt, release_hash, manifest_hash, row):
    try:
        record = json.loads((Path(attempt) / 'status.json').read_text())
        archive = Path(attempt) / 'evidence.tar.gz'
        return (record.get('status') == 'completed'
                and record.get('release_lock_sha256') == release_hash
                and record.get('input_manifest_sha256') == manifest_hash
                and record.get('input_xml_sha256') == row['xml_sha256']
                and archive.is_file() and sha(archive) == record.get('archive_sha256'))
    except (OSError, ValueError):
        return False
