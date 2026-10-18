"""Create a portable gzipped input set from an audited local manifest.

Only original models are copied. No original file is modified or deleted.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

from common import read_csv, sha, write_csv, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--pilot-manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--tar', action='store_true')
    args = p.parse_args()
    source = read_csv(args.manifest)
    pilot = {r['stable_model_id'] for r in read_csv(args.pilot_manifest)}
    assert len(source) == 6543 and len(pilot) == 24
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []
    for i, r in enumerate(source):
        original = Path(r['before_path'])
        if sha(original) != r['before_sha256']:
            raise ValueError('Input changed since audit: ' + str(original))
        job_id = '{}_{:06d}'.format(r['database'].lower(), i)
        relative = 'models/{}/{}.xml.gz'.format(r['database'], job_id)
        dest = args.out / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        if original.suffix == '.gz':
            shutil.copyfile(str(original), str(dest))
        else:
            with original.open('rb') as src, dest.open('wb') as raw:
                with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as gz:
                    shutil.copyfileobj(src, gz)
        h = hashlib.sha256()
        with gzip.open(dest, 'rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        if h.hexdigest() != r['decompressed_xml_sha256']:
            raise ValueError('Transferred XML changed: ' + job_id)
        rows.append(dict(job_id=job_id, database=r['database'], stable_model_id=r['stable_model_id'],
            original_model_id=r['archive_model_id'], input_relative_path=relative,
            input_filename=original.name.removesuffix('.gz') if hasattr(str, 'removesuffix')
                           else (original.name[:-3] if original.name.endswith('.gz') else original.name),
            gzip_sha256=sha(dest), xml_sha256=h.hexdigest(),
            historical_after_sha256=r['historical_after_sha256'],
            reactions=r['reactions'], metabolites=r['metabolites'],
            pilot=int(r['stable_model_id'] in pilot)))
        if (i+1) % 1000 == 0:
            print('Prepared', i+1, flush=True)
    write_csv(args.out/'manifest.csv', rows, list(rows[0]))
    write_json(args.out/'bundle.json', dict(schema=1, model_count=len(rows),
        database_counts=dict(Counter(r['database'] for r in rows)),
        manifest_sha256=sha(args.out/'manifest.csv'),
        source_audit_manifest_sha256=sha(args.manifest),
        cohort='manuscript_core_6543', contains='original_models_only', pilot_count=24))
    if args.tar:
        target = args.out.with_suffix('.tar')
        with tarfile.open(str(target), 'x') as tar:
            tar.add(str(args.out), arcname=args.out.name)
        target.with_suffix('.tar.sha256').write_text(sha(target)+'  '+target.name+'\n')
        print('Transfer archive:', target, 'bytes:', target.stat().st_size)


if __name__ == '__main__':
    main()
