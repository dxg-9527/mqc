"""Record dependency versions and test that CPLEX accepts a >1000-variable LP."""
import argparse
import importlib.metadata
import json
import platform
import sys

from common import write_json

EXPECTED = {'cobra':'0.29.1', 'optlang':'1.8.2', 'cplex':'22.1.0.0',
    'python-libsbml':'5.20.1', 'numpy':'1.22.4', 'pandas':'1.5.3',
    'sympy':'1.12', 'swiglpk':'5.0.8', 'openpyxl':'3.1.2',
    'matplotlib':'3.7.5', 'memote':'0.17.0', 'd3flux':'0.2.7'}


def inspect(probe=False):
    packages = {}
    for name in EXPECTED:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    differences = {k:dict(expected=v, observed=packages[k])
                   for k,v in EXPECTED.items() if packages[k]!=v}
    if platform.python_version() != '3.8.16':
        differences['python'] = dict(expected='3.8.16', observed=platform.python_version())
    record = dict(python=platform.python_version(), packages=packages,
                  differences=differences, cplex_probe='not_run',
                  all_installed_versions={d.metadata['Name']: d.version
                    for d in importlib.metadata.distributions() if d.metadata.get('Name')})
    if probe:
        try:
            import cobra
            import cplex
            default_solver = cobra.Configuration().solver
            record['cobra_default_solver'] = str(default_solver)
            if default_solver.__name__ != 'optlang.cplex_interface':
                raise RuntimeError('The validated environment uses CPLEX as the COBRA default; '
                                   'use a dedicated environment without another preferred solver')
            problem = cplex.Cplex()
            for method in ('set_log_stream','set_results_stream','set_warning_stream','set_error_stream'):
                getattr(problem, method)(None)
            problem.parameters.threads.set(1)
            problem.variables.add(obj=[1.]*1101, lb=[0.]*1101, ub=[2.]*1101)
            problem.linear_constraints.add(
                lin_expr=[cplex.SparsePair(ind=[i],val=[1.]) for i in range(1101)],
                senses='G'*1101, rhs=[1.]*1101)
            problem.solve()
            value = problem.solution.get_objective_value()
            if abs(value - 1101.) > 1e-6:
                raise RuntimeError('Unexpected LP result')
            record.update(cplex_probe='passed', cobra_default_solver=str(cobra.Configuration().solver))
            problem.end()
        except Exception as exc:
            record.update(cplex_probe='failed', probe_error=repr(exc))
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True)
    p.add_argument('--probe', action='store_true')
    p.add_argument('--allow-environment-drift', action='store_true')
    args=p.parse_args()
    record=inspect(args.probe)
    write_json(args.output,record)
    print(json.dumps(record,indent=2))
    missing=any(v is None for v in record['packages'].values())
    bad=missing or record['cplex_probe']=='failed' or (record['differences'] and not args.allow_environment_drift)
    raise SystemExit(int(bool(bad)))


if __name__=='__main__': main()
