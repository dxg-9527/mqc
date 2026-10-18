"""Summarize all expected models, retain failed attempts, and prepare return bundles."""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import tarfile
import time

from common import load_inputs, read_csv, sha, valid_completed, verify_release, write_csv, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--pack',action='store_true',help='After jobs stop: create report tar and per-database evidence tars')
    p.add_argument('--database',choices=['BiGG','CarveMe','ModelSEED','VMH'],help='Pack selected DB; summary still covers all 6543')
    args=p.parse_args()
    release=verify_release()
    rows,meta=load_inputs(args.inputs)
    root=args.results.resolve()
    identity=json.loads((root/'run_identity.json').read_text())
    if identity['release_lock_sha256']!=release or identity['input_manifest_sha256']!=meta['manifest_sha256']:
        raise ValueError('Source or input identity differs from this results directory')
    report=root/'reports'/time.strftime('%Y%m%dT%H%M%S')
    report.mkdir(parents=True,exist_ok=False)
    summary,attempts,selected=[],[],{}
    for row in rows:
        candidates=[]
        for attempt in sorted((root/'tasks'/row['job_id']).glob('attempt_*')):
            try: record=json.loads((attempt/'status.json').read_text())
            except (OSError,ValueError): record=dict(status='missing_or_corrupt_status')
            good=valid_completed(attempt,release,meta['manifest_sha256'],row)
            item=dict(job_id=row['job_id'],database=row['database'],attempt=str(attempt.relative_to(root)),
                      status=record.get('status'),verified_completed=good,
                      seconds=record.get('seconds'),archive_sha256=record.get('archive_sha256'))
            attempts.append(item)
            if good: candidates.append((attempt,record))
        result=dict(job_id=row['job_id'],database=row['database'],stable_model_id=row['stable_model_id'],
                    status='missing_or_failed',input_xml_sha256=row['xml_sha256'])
        if candidates:
            # Keep every attempt. First valid completion is the predeclared representative.
            attempt,record=candidates[0]
            outcome=record.get('model_outcome',{})
            result.update(status='completed',attempt=str(attempt.relative_to(root)),
                seconds=record.get('seconds'),final_model_sha256=outcome.get('final_model_sha256'),
                final_objective_value=outcome.get('final_objective_value'),
                postgap_status=outcome.get('postgap_status'),
                max_abs_Sv=outcome.get('final_pfba_max_abs_Sv'),
                completed_attempts=len(candidates))
            selected[row['job_id']]=attempt
        summary.append(result)
    fields=['job_id','database','stable_model_id','status','input_xml_sha256','attempt','seconds',
            'final_model_sha256','final_objective_value','postgap_status','max_abs_Sv','completed_attempts']
    write_csv(report/'summary.csv',summary,fields)
    write_csv(report/'rerun.csv',[r for r in summary if r['status']!='completed'],fields)
    write_csv(report/'attempts.csv',attempts,
              ['job_id','database','attempt','status','verified_completed','seconds','archive_sha256'])
    plans=[]
    for path in sorted((root/'plans').glob('*.json')):
        plans.append(json.loads(path.read_text()))
    job_ids=sorted({x['slurm_job_id'] for x in plans if x.get('slurm_job_id')})
    if job_ids:
        try:
            response=subprocess.run(['sacct','-j',','.join(job_ids),'-P',
                '--format=JobID,State,ExitCode,Elapsed,AllocCPUS,MaxRSS,ReqMem'],
                capture_output=True,text=True,timeout=120)
            (report/'sacct.txt').write_text(response.stdout)
            (report/'sacct_error.txt').write_text(response.stderr)
        except (OSError,subprocess.TimeoutExpired) as exc:
            (report/'sacct_error.txt').write_text(repr(exc))
    # Resource accounting may omit killed processes; retain Slurm logs and missing jobs.
    totals=dict(expected=len(rows),completed=len(selected),not_completed=len(rows)-len(selected),
        by_database={db:dict(Counter(r['status'] for r in summary if r['database']==db))
                     for db in sorted({r['database'] for r in rows})},
        submitted_job_ids=job_ids,release_lock_sha256=release,input_manifest_sha256=meta['manifest_sha256'],
        interpretation='completed means workflow/export/archive verification, not demonstrated experimental accuracy')
    write_json(report/'counts.json',totals)
    print(json.dumps(totals,indent=2))
    print('Reports:',report)
    if args.pack:
        destination=root/'return_packages'/report.name
        destination.mkdir(parents=True,exist_ok=False)
        archives=[]
        # Do not extract/recompress per-model gz archives; keep failures and prior attempts too.
        dbs=[args.database] if args.database else sorted({r['database'] for r in rows})
        for db in dbs:
            path=destination/('MQC_optimized_'+db+'_evidence.tar')
            with tarfile.open(str(path),'x') as tar:
                for row in rows:
                    task=root/'tasks'/row['job_id']
                    if row['database']==db and task.exists():
                        tar.add(str(task),arcname='tasks/'+row['job_id'])
            archives.append(path)
        path=destination/'MQC_optimized_reports.tar.gz'
        with tarfile.open(str(path),'x:gz') as tar:
            for name in ['reports','plans','environments','slurm_logs','run_identity.json']:
                f=root/name
                if f.exists(): tar.add(str(f),arcname=name)
            for name in ['manifest.csv','bundle.json']:
                tar.add(str(args.inputs/name),arcname='inputs/'+name)
        archives.append(path)
        checks=[dict(file=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in archives]
        write_csv(destination/'return_manifest.csv',checks,['file','bytes','sha256'])
        (destination/'SHA256SUMS.txt').write_text(''.join(r['sha256']+'  '+r['file']+'\n' for r in checks))
        print('Return every archive + SHA256SUMS.txt + return_manifest.csv from:',destination)
    return 0


if __name__=='__main__': raise SystemExit(main())
