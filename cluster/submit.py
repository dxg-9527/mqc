"""Submit one Slurm array, grouping models only when array size is limited."""
import argparse
import heapq
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import uuid

from common import REPO, git_commit, load_inputs, read_csv, verify_release, write_json


def assignments(rows, slots):
    """Greedy size balancing; no phenotype or quality outcomes used."""
    groups = [[] for _ in range(slots)]
    heap = [(0, i) for i in range(slots)]
    heapq.heapify(heap)
    for r in sorted(rows, key=lambda r: (-int(r['reactions']), r['job_id'])):
        weight, i = heapq.heappop(heap)
        groups[i].append(r['job_id'])
        heapq.heappush(heap, (weight + int(r['reactions']), i))
    return groups


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--python', default=sys.executable)
    p.add_argument('--partition', default='qcpu_23i')
    p.add_argument('--account')
    p.add_argument('--parallel', type=int, default=32)
    p.add_argument('--cpus', type=int, default=10)
    p.add_argument('--mem', default='8G')
    p.add_argument('--time', default='24:00:00')
    p.add_argument('--array-slots', type=int, default=1000)
    p.add_argument('--model-timeout', type=int, default=7200)
    p.add_argument('--pilot', action='store_true')
    p.add_argument('--database', choices=['BiGG','CarveMe','ModelSEED','VMH'])
    p.add_argument('--ids-file', type=Path, help='rerun.csv from collect.py; reads job_id')
    p.add_argument('--allow-environment-drift', action='store_true')
    p.add_argument('--submit', action='store_true', help='without this flag, only writes a plan')
    args = p.parse_args()
    if min(args.parallel, args.cpus, args.array_slots, args.model_timeout) < 1:
        p.error('Counts and timeout must be positive')
    release = verify_release()
    rows, meta = load_inputs(args.inputs)
    selected = rows
    if args.pilot:
        selected = [r for r in selected if r['pilot']=='1']
    if args.database:
        selected = [r for r in selected if r['database']==args.database]
    if args.ids_file:
        ids = {r['job_id'] for r in read_csv(args.ids_file)}
        if not ids <= {r['job_id'] for r in rows}:
            p.error('Retry list contains unknown job IDs')
        selected = [r for r in selected if r['job_id'] in ids]
    if not selected:
        p.error('No selected models')
    cap = args.array_slots
    slurm_config = 'unavailable; using explicit --array-slots'
    try:
        config = subprocess.run(['scontrol', 'show', 'config'], capture_output=True, text=True, timeout=20)
        match = re.search(r'\bMaxArraySize\s*=\s*(\d+)', config.stdout)
        if match:
            cap = min(cap, int(match.group(1)))
            slurm_config = match.group(0)
        match = re.search(r'\bmax_array_tasks=(\d+)', config.stdout)
        if match:
            cap = min(cap, int(match.group(1)))
    except (OSError, subprocess.TimeoutExpired):
        pass
    if cap < 1:
        p.error('Cluster job arrays appear disabled')
    groups = assignments(selected, min(cap, len(selected)))
    root = args.results.resolve()
    (root/'plans').mkdir(parents=True, exist_ok=True)
    (root/'slurm_logs').mkdir(exist_ok=True)
    identity = dict(release_lock_sha256=release, input_manifest_sha256=meta['manifest_sha256'],
                    expected_total=meta['model_count'])
    identity_path = root/'run_identity.json'
    if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
        p.error('Results directory belongs to another release/input cohort; use a new directory')
    write_json(identity_path, identity)
    plan_path = root/'plans'/('plan_'+uuid.uuid4().hex[:12]+'.json')
    command = ['sbatch', '--parsable', '--job-name=mqc-opt', '--partition='+args.partition,
        '--nodes=1', '--ntasks=1', '--cpus-per-task='+str(args.cpus), '--mem='+args.mem,
        '--time='+args.time, '--array=0-{}%{}'.format(len(groups)-1, min(args.parallel, len(groups))),
        '--output='+str(root/'slurm_logs/%A_%a.out'),
        '--error='+str(root/'slurm_logs/%A_%a.err')]
    if args.account:
        command.append('--account='+args.account)
    command += [str(REPO/'cluster/array.sbatch'), str(Path(args.python).resolve()), str(REPO), str(plan_path)]
    plan = dict(schema=1, **identity, git_commit=git_commit(),
        inputs=str(args.inputs.resolve()), results=str(root), assignments=groups,
        selected_count=len(selected), array_slots=len(groups),
        max_models_per_slot=max(map(len, groups)), concurrency_limit=min(args.parallel,len(groups)),
        model_timeout_seconds=args.model_timeout, allow_environment_drift=args.allow_environment_drift,
        command=command, slurm_array_limit=slurm_config, submitted=False)
    write_json(plan_path, plan)
    print('Models:',len(selected),'array tasks:',len(groups),'maximum models/task:',plan['max_models_per_slot'])
    print('Plan:',plan_path)
    print(' '.join(shlex.quote(x) for x in command))
    if args.submit:
        response = subprocess.run(command, capture_output=True, text=True)
        plan['sbatch_stdout'], plan['sbatch_stderr'] = response.stdout, response.stderr
        plan['submitted'] = response.returncode == 0
        if response.returncode == 0:
            plan['slurm_job_id'] = response.stdout.strip().split(';')[0]
        write_json(plan_path, plan)
        print(response.stdout, response.stderr)
        raise SystemExit(response.returncode)
    print('Dry run only. Repeat with --submit to submit on the cluster.')


if __name__ == '__main__':
    main()
