"""Run the unchanged optimized correction once, observing six stage snapshots."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

from common import REPO, STAGES, git_commit, sha, verify_release, write_json


def load_exported_model(path):
    import cobra
    path = Path(path)
    # MQC may export a valid SBML file named '.xml'; Path('.xml').suffix is empty.
    if path.name.lower().endswith('.xml'):
        return cobra.io.read_sbml_model(str(path))
    return cobra.io.load_json_model(str(path))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    args=p.parse_args()
    started=time.time()
    work=args.work.resolve()
    output=work/'output'
    record=dict(status='started', started=started,
                mode='model_check2_only', experimental_values_used=False,
                argv=[str(args.input),str(output)], stages=[])
    write_json(work/'outcome.json',record)
    sys.path.insert(0,str(REPO))
    try:
        record.update(git_commit=git_commit(), release_lock_sha256=verify_release(),
                      input_xml_sha256=sha(args.input))
        write_json(work/'outcome.json', record)
        import cobra
        from mqc.solver_utils import robust_pfba as pfba
        from cobra.util.solver import linear_reaction_coefficients
        from mqc.main import ModelCheck
        from mqc.control.auxiliary_supply import validate_supply_budget
        from mqc.control import initial_control,nadh_control,atp_control,net_control,yield_control,biomass_control
        classes=[initial_control.InitialPreprocess,nadh_control.Nadhs,atp_control.Atps,
                 net_control.Nets,yield_control.Yields,biomass_control.Biomasses]
        for cls,name in zip(classes,STAGES):
            original=getattr(cls,name)
            def observed(self,*a,_call=original,_name=name,**kw):
                record['stages'].append(dict(stage=_name,event='start',seconds=time.time()-started))
                write_json(work/'outcome.json',record)
                value=_call(self,*a,**kw)
                # Serialization only; no additional optimization inside the workflow.
                cobra.io.save_json_model(a[1],str(output/(_name+'_model.json')))
                record['stages'].append(dict(stage=_name,event='finished',seconds=time.time()-started))
                write_json(work/'outcome.json',record)
                return value
            setattr(cls,name,observed)
        job=ModelCheck(str(args.input),str(output))
        write_json(work/'run_config.json',dict(config=vars(job.cfg), **record))
        result=job.model_check2()  # Never call run(), which also performs detection.
        finished={r['stage'] for r in record['stages'] if r['event']=='finished'}
        if finished!=set(STAGES) or not isinstance(result,(list,tuple)) or len(result)!=3:
            raise RuntimeError('Correction returned without all six stages and final paths')
        exported=Path(result[1]).resolve()
        exported.relative_to(output)
        model=load_exported_model(exported)
        control=json.loads((output/'control_result.json').read_text())
        supply = control.get('auxiliary_supply_budget', {})
        validate_supply_budget(model, supply)
        record['auxiliary_supply_status'] = supply.get('status', 'not_reported')
        record['auxiliary_carbon_accounting'] = supply.get('carbon_accounting', 'not_reported')
        record['auxiliary_partial_formula_inputs'] = supply.get('partial_formula_inputs', {})
        if control.get('postgap_validation',{}).get('status') not in ('passed','repaired'):
            raise RuntimeError('Final post-gap validation did not pass')
        solution=model.optimize()
        record.update(final_model=str(exported.relative_to(work)),
            final_model_sha256=sha(exported), final_reactions=len(model.reactions),
            final_metabolites=len(model.metabolites), final_solver_status=solution.status,
            final_objective={r.id:float(v) for r,v in linear_reaction_coefficients(model).items()},
            final_objective_value=float(solution.objective_value) if solution.status=='optimal' else None,
            postgap_status=control['postgap_validation']['status'],
            postgap_energy_status=control['postgap_validation'].get('energy_status'),
            final_growth_status='not_positive')
        if (not record['final_objective'] or solution.status != 'optimal'
                or not math.isfinite(record['final_objective_value'])
                or record['final_objective_value'] <= 1e-7):
            raise RuntimeError('Final exported model has no positive growth in its stored QC medium')
        record['final_growth_status']='positive'
        if solution.status=='optimal':
            flux=pfba(model).fluxes
            residuals={m.id:0. for m in model.metabolites}
            violation=0.
            with (output/'final_pfba.tsv').open('w',newline='') as stream:
                w=csv.writer(stream,delimiter='\t')
                w.writerow(['reaction','flux','lower_bound','upper_bound','equation','GPR'])
                for r in model.reactions:
                    v=float(flux[r.id])
                    violation=max(violation,r.lower_bound-v,v-r.upper_bound)
                    for m,c in r.metabolites.items(): residuals[m.id]+=v*c
                    if abs(v)>1e-9: w.writerow([r.id,v,r.lower_bound,r.upper_bound,r.reaction,r.gene_reaction_rule])
            record.update(final_pfba_max_abs_Sv=max(map(abs,residuals.values()),default=0.),
                          final_pfba_max_bound_violation=violation)
            if record['final_pfba_max_abs_Sv']>1e-6 or violation>1e-6:
                raise RuntimeError('Reopened final-model flux failed numerical checks')
        else:
            raise RuntimeError('Final model reloads, but objective solve is '+solution.status)
        verify_release()
        record['status']='completed'
    except BaseException as exc:
        record.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
        traceback.print_exc()
    finally:
        record['seconds']=time.time()-started
        write_json(work/'outcome.json',record)
    return 0 if record['status']=='completed' else 1


if __name__=='__main__': raise SystemExit(main())
