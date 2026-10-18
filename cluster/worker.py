"""Execute one array slot, isolate each model, and retain failures as evidence."""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time
import traceback
import uuid

from common import REPO, load_inputs, sha, under, valid_completed, verify_release, write_json


def compress_evidence(work, attempt):
    target=attempt/'evidence.tar.gz'
    temporary=attempt/'evidence.tar.gz.partial'
    with tarfile.open(str(temporary),'w:gz',compresslevel=6) as tar:
        for child in sorted(work.iterdir()):
            if child.name!='input': tar.add(str(child),arcname=child.name)
    # Stream every member to validate decompression before reclaiming scratch.
    with tarfile.open(str(temporary),'r:gz') as tar:
        for member in tar:
            if member.isfile():
                stream=tar.extractfile(member)
                while stream.read(1024*1024): pass
    os.replace(str(temporary),str(target))
    return sha(target)


def terminate_group(process):
    try: os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError: return
    try: process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try: os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError: pass
        process.wait()


def run_one(row,plan,environment):
    root=Path(plan['results'])
    model_root=root/'tasks'/row['job_id']
    model_root.mkdir(parents=True,exist_ok=True)
    for previous in sorted(model_root.glob('attempt_*')):
        if valid_completed(previous,plan['release_lock_sha256'],plan['input_manifest_sha256'],row):
            print('SKIP verified completed',row['job_id'],flush=True)
            return True
    attempt=model_root/('attempt_'+time.strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8])
    attempt.mkdir()
    record=dict(job_id=row['job_id'],database=row['database'],stable_model_id=row['stable_model_id'],
        release_lock_sha256=plan['release_lock_sha256'], input_manifest_sha256=plan['input_manifest_sha256'],
        input_xml_sha256=row['xml_sha256'], status='started',started=time.time(),
        slurm_job_id=os.environ.get('SLURM_JOB_ID'),slurm_array_job_id=os.environ.get('SLURM_ARRAY_JOB_ID'),
        slurm_array_task_id=os.environ.get('SLURM_ARRAY_TASK_ID'))
    write_json(attempt/'status.json',record)
    scratch=Path(os.environ.get('SLURM_TMPDIR',str(root/'scratch')))
    scratch.mkdir(parents=True,exist_ok=True)
    work=Path(tempfile.mkdtemp(prefix='mqc_'+row['job_id']+'_',dir=str(scratch)))
    record['scratch']=str(work)
    write_json(attempt/'status.json',record)
    write_json(work/'environment.json',environment)
    write_json(work/'input_record.json',row)
    try:
        (work/'input').mkdir()
        filename=Path(row['input_filename']).name
        if filename!=row['input_filename'] or not filename.endswith('.xml'):
            raise ValueError('Unsafe/non-XML input filename')
        original=under(plan['inputs'],row['input_relative_path'])
        if sha(original)!=row['gzip_sha256']: raise ValueError('Gzip input checksum mismatch')
        model=work/'input'/filename
        with gzip.open(original,'rb') as src,model.open('wb') as dst: shutil.copyfileobj(src,dst)
        if sha(model)!=row['xml_sha256']: raise ValueError('Uncompressed input checksum mismatch')
        env=dict(os.environ,PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',PYTHONUNBUFFERED='1',
            OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg',
            MPLCONFIGDIR=str(work/'mpl_cache'))
        command=[sys.executable,str(REPO/'cluster/run_model.py'),'--input',str(model),'--work',str(work)]
        record['command']=command
        write_json(attempt/'status.json',record)
        with (work/'pipeline.log').open('w') as log:
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
            try: code=process.wait(timeout=plan['model_timeout_seconds'])
            except subprocess.TimeoutExpired:
                terminate_group(process)
                code=124
            except BaseException:
                terminate_group(process)
                raise
        child=json.loads((work/'outcome.json').read_text()) if (work/'outcome.json').exists() else {}
        record['returncode']=code
        record['model_outcome']=child
        record['status']='completed' if code==0 and child.get('status')=='completed' else ('timeout' if code==124 else 'failed')
    except BaseException as exc:
        record.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
        write_json(work/'worker_error.json',record)
    record['seconds']=time.time()-record['started']
    write_json(work/'worker_status.json',record)
    try:
        record['archive_sha256']=compress_evidence(work,attempt)
        record['archive_bytes']=(attempt/'evidence.tar.gz').stat().st_size
        write_json(attempt/'status.json',record)
        # Only this function's unique mkdtemp directory, after verified archival.
        if work.parent.resolve()!=scratch.resolve() or not work.name.startswith('mqc_'+row['job_id']+'_'):
            raise ValueError('Unsafe scratch cleanup target')
        shutil.rmtree(str(work))
    except Exception as exc:
        record.update(status='archive_or_cleanup_failed',archive_error=repr(exc))
        write_json(attempt/'status.json',record)
    print(row['job_id'],record['status'],round(record['seconds'],1),'seconds',flush=True)
    return record['status']=='completed'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--slot',type=int,required=True)
    args=p.parse_args()
    plan=json.loads(args.plan.read_text())
    if verify_release()!=plan['release_lock_sha256']: raise ValueError('Release changed since submission')
    rows,meta=load_inputs(plan['inputs'])
    if meta['manifest_sha256']!=plan['input_manifest_sha256']: raise ValueError('Input manifest changed')
    if not 0<=args.slot<len(plan['assignments']): raise ValueError('Invalid array index')
    rows={r['job_id']:r for r in rows}
    envpath=Path(plan['results'])/'environments'/(args.plan.stem+'_'+str(args.slot)+'_'+uuid.uuid4().hex[:8]+'.json')
    command=[sys.executable,str(REPO/'cluster/check_environment.py'),'--output',str(envpath),'--probe']
    if plan['allow_environment_drift']: command.append('--allow-environment-drift')
    result=subprocess.run(command)
    if result.returncode: return result.returncode
    environment=json.loads(envpath.read_text())
    results=[run_one(rows[key],plan,environment) for key in plan['assignments'][args.slot]]
    return int(not all(results))


if __name__=='__main__': raise SystemExit(main())
