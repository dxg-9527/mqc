"""Bounded native-network coverage inside the original precursor repair.

Targets come from blocked, original GPR-bearing reactions, not phenotypes.
Only structural copies open existing exchanges. The caller inserts the safe
union with its ordinary gapfill insertion routine before biomass finalization.
This is an experimental algorithm extension, not a validated enzyme assignment.
"""
import hashlib
import json
import math
import re
from pathlib import Path

import cobra
from cobra.exceptions import OptimizationError
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis
from cobra.medium import find_external_compartment
from mqc.control.postgap_control import close_inputs, diagnose, discover_targets
from mqc.solver_utils import solve_with_cold_retry


def native_core_reactions(model, original_ids, biomass):
    external = find_external_compartment(model)
    core = sorted(r.id for r in model.reactions if r.id in original_ids
        and r.id != biomass and len(r.metabolites) > 1
        and bool(r.gene_reaction_rule.strip()) and len(r.compartments) == 1
        and external not in r.compartments and r.lower_bound <= 0 <= r.upper_bound
        and r.bounds != (0, 0))
    return core, external


def is_annotated_isomerase(reaction):
    """Conservative eligibility, not proof of a strain's enzyme or direction.

    Original reversible bounds, an unambiguous EC-5 annotation and a balanced
    single-compartment 1:1 isomer conversion are all necessary. No ID/name or
    phenotype whitelist. GPR/native membership is enforced by the caller.
    """
    if not reaction.lower_bound < 0 < reaction.upper_bound or len(reaction.metabolites) != 2:
        return False
    if len(reaction.compartments) != 1:
        return False
    codes = reaction.annotation.get('ec-code', [])
    codes = [codes] if isinstance(codes, str) else codes
    if not isinstance(codes, (list, tuple)) or not codes:
        return False
    if not all(isinstance(c,str) and re.fullmatch(r'5\.\d+\.\d+\.\d+',c.strip()) for c in codes):
        return False
    (a, x), (b, y) = reaction.metabolites.items()
    if x*y >= 0 or not math.isfinite(x+y) or abs(x+y) > 1e-12:
        return False
    if not a.formula or not b.formula or not a.elements or a.elements != b.elements:
        return False
    if any(k in {'R','X'} or not math.isfinite(v) or v <= 0 for k,v in a.elements.items()):
        return False
    if a.charge is None or b.charge is None:
        return False
    return math.isfinite(a.charge) and math.isfinite(b.charge) and a.charge == b.charge


def coverage_directions(model, core, fva):
    """Preserve complete-block targets; add only missing isomerase directions."""
    targets, additional = {}, []
    for rid in sorted(core):
        low, high = float(fva.loc[rid,'minimum']), float(fva.loc[rid,'maximum'])
        if not math.isfinite(low) or not math.isfinite(high):
            raise ValueError('Nonfinite coverage capacity: '+rid)
        if abs(low) < 1e-7 and abs(high) < 1e-7:
            targets[rid] = [1,-1]
        elif is_annotated_isomerase(model.reactions.get_by_id(rid)):
            signs = ([1] if abs(high) < 1e-7 else []) + ([-1] if abs(low) < 1e-7 else [])
            if signs:
                targets[rid] = signs
                additional.append(dict(reaction=rid,signs=signs,minimum=low,maximum=high,
                    ec_codes=model.reactions.get_by_id(rid).annotation['ec-code']))
    return targets, additional


def open_existing_exchanges(model, ids):
    changes = []
    for rid in ids:
        reaction = model.reactions.get_by_id(rid)
        if not reaction.lower_bound <= 0 <= reaction.upper_bound:
            raise ValueError('Forced exchange cannot be opened for coverage: ' + rid)
        before = list(reaction.bounds)
        reaction.bounds = (-1000, 1000)
        if before != list(reaction.bounds):
            changes.append(dict(reaction=rid, before=before, after=list(reaction.bounds)))
    return changes


def chemical_gate(model, library, additions):
    if not additions:
        return 'no_new_connection'
    if len(additions) > 3:
        return 'more_than_three_new_reactions'
    external = find_external_compartment(model)
    for rid in additions:
        r = library.reactions.get_by_id(rid)
        if len(r.metabolites) <= 1 or len(r.compartments) != 1:
            return 'boundary_or_transport_addition'
        if any(m.id not in model.metabolites for m in r.metabolites):
            return 'requires_new_metabolite'
        resolved = [model.metabolites.get_by_id(m.id) for m in r.metabolites]
        if len({m.compartment for m in resolved}) != 1 or external in {m.compartment for m in resolved}:
            return 'external_or_transport_addition'
        if any(not m.formula or m.charge is None or not m.elements
               or any(k in {'R', 'X'} for k in m.elements) for m in resolved):
            return 'incomplete_metabolite_chemistry'
    with model:
        model.add_reactions([library.reactions.get_by_id(r).copy() for r in additions])
        for rid in additions:
            balance = model.reactions.get_by_id(rid).check_mass_balance()
            if any(not math.isfinite(v) or abs(v) > 1e-7 for v in balance.values()):
                return 'atom_or_charge_imbalance'
    return 'eligible'


def energy_checks(model, biomass):
    targets, coverage = discover_targets(model)
    checks = [diagnose(model, target) for target in targets]
    with model:
        close_inputs(model)
        model.objective = biomass
        model.objective_direction = 'max'
        sol = model.optimize()
        checks.append(dict(target='Growth', solver_status=sol.status,
            value=sol.objective_value, passed=sol.status == 'optimal'
            and math.isfinite(sol.objective_value) and abs(sol.objective_value) < 1e-7))
    # Absence of recognizable carrier tests must not be reported as success.
    available_groups = {t['group'] for t in targets}
    passed = {'energy', 'redox'} <= available_groups and all(c['passed'] for c in checks)
    return dict(passed=passed, checks=checks, coverage=coverage)


def target_capacities(model, targets, exchanges=None):
    rows = []
    with model:
        if exchanges is None:
            close_inputs(model)
        else:
            open_existing_exchanges(model, exchanges)
        for target in targets:
            model.objective = target['reaction']
            model.objective_direction = 'max' if target['sign'] > 0 else 'min'
            sol = model.optimize()
            rows.append(dict(**target, status=sol.status,
                capacity=float(sol.objective_value) if sol.status == 'optimal' else None))
    return rows


def safe_targets(model, targets, exchanges):
    opened = target_capacities(model, targets, exchanges)
    closed = target_capacities(model, targets)
    accepted = [t for t, a, b in zip(targets, opened, closed)
        if a['status'] == b['status'] == 'optimal'
        and math.isfinite(a['capacity']) and math.isfinite(b['capacity'])
        and a['capacity'] * t['sign'] >= .1 - 1e-7 and abs(b['capacity']) < 1e-7]
    return accepted, dict(opened=opened, source_free=closed)


def set_weighted_objective(model, costs):
    model.objective = model.problem.Objective(0, direction='min')
    coefficients = {}
    for r in model.reactions:
        cost = costs.get(r.id, 1.)
        if not math.isfinite(cost) or cost < 1:
            raise ValueError('Invalid original gapfill cost: ' + r.id)
        coefficients[r.forward_variable] = coefficients[r.reverse_variable] = cost
    model.objective.set_linear_coefficients(coefficients)


def _settings(model):
    model.solver.configuration.timeout = 90
    if 'cplex' in model.solver.interface.__name__:
        model.solver.problem.parameters.threads.set(1)


def select_native_coverage(model, big_model, info, directory):
    """Select a safe union without mutating the actual model or search object.

    Requires the original S*1000 precursor search still in memory. Native
    target flux 0.0001 there means physical flux 0.1. Original EX_gap optimum
    and candidate penalties are retained; no new boundary is constructed.
    """
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=False)
    report = dict(status='started', added=[], experimental_values_used=False,
        source='original_precursor_search_in_memory', native_search_scale=1000,
        target_physical_flux=.1, blocked_tolerance=1e-7, max_new_reactions_per_plan=3,
        screening_environment='existing_exchanges_open_in_structural_copies_only')

    def finish(status, additions=()):
        report.update(status=status, selected_union=list(additions))
        (folder/'report.json').write_text(json.dumps(report, indent=2))
        return list(additions), report

    candidates = set(info.get('gapfill_candidate_bounds', {}))
    if 'EX_gap' not in big_model.reactions or not candidates:
        return finish('skipped_no_original_precursor_search')
    biomass = info['initial_rxn']['biomass_rxn_id']
    # Retain current optlang constraints by copying, rather than round-tripping
    # a model file (JSON/SBML would omit non-stoichiometric constraints).
    native = model.copy()
    search = big_model.copy()
    _settings(native)
    _settings(search)
    native.objective = biomass
    native.objective_direction = 'max'
    baseline = native.optimize()
    report['recorded_medium_growth_before'] = baseline.objective_value
    if baseline.status != 'optimal' or not math.isfinite(baseline.objective_value) or baseline.objective_value <= 1e-7:
        return finish('skipped_original_precursor_repair_not_viable')
    report['baseline_safety'] = energy_checks(native, biomass)
    if not report['baseline_safety']['passed']:
        return finish('skipped_preexisting_source_free_or_unsupported_safety')
    original_ids = set(search.reactions.list_attr('id')) - candidates - {'EX_gap'}
    exchanges = sorted(r.id for r in native.exchanges)
    core, external = native_core_reactions(native, original_ids, biomass)
    report.update(core_count=len(core), external_compartment=external,
        original_ids=sorted(original_ids), existing_exchanges=exchanges)
    structure = native.copy()
    report['diagnostic_environment_changes'] = open_existing_exchanges(structure, exchanges)
    structure.objective = structure.problem.Objective(0)
    if not core:
        return finish('no_eligible_native_reactions')
    fva = flux_variability_analysis(structure, reaction_list=core, fraction_of_optimum=0, processes=1)
    if not all(math.isfinite(v) for v in fva[['minimum', 'maximum']].to_numpy().ravel()):
        raise RuntimeError('Nonfinite native coverage FVA')
    fva['blocked'] = (fva.minimum.abs() < 1e-7) & (fva.maximum.abs() < 1e-7)
    fva['equation'] = [native.reactions.get_by_id(r).reaction for r in fva.index]
    fva['gpr'] = [native.reactions.get_by_id(r).gene_reaction_rule for r in fva.index]
    fva.to_csv(folder/'native_core_screen.csv', index_label='reaction')
    blocked = sorted(fva.index[fva.blocked])
    report['blocked_core'] = blocked
    directions, additional = coverage_directions(native, core, fva)
    report['directionally_blocked_isomerases'] = additional
    report['target_directions'] = directions
    print('Native coverage blocked targets:', blocked, flush=True)
    print('Native coverage missing isomerase directions:', additional, flush=True)
    if not directions:
        return finish('no_blocked_native_reactions')
    search.objective = 'EX_gap'
    search.objective_direction = 'max'
    try:
        search, primary, report['precursor_solver'] = solve_with_cold_retry(search)
    except OptimizationError as error:
        report['precursor_solver'] = dict(status='unresolved', error=repr(error))
        return finish('skipped_unresolved_original_precursor_solve')
    if not math.isfinite(primary) or primary <= 1e-7:
        return finish('skipped_infeasible_original_precursor_target')
    report['precursor_optimum_preserved'] = primary
    open_existing_exchanges(search, exchanges)
    search.reactions.get_by_id('EX_gap').bounds = (primary, primary)
    existing = set(native.reactions.list_attr('id'))
    plans, proposals = {}, []
    for i, rid in enumerate(sorted(directions)):
        reaction = search.reactions.get_by_id(rid)
        for sign in directions[rid]:
            direction = 'max' if sign > 0 else 'min'
            search.objective = reaction
            search.objective_direction = direction
            capacity = search.optimize()
            proposal = dict(target=rid, sign=sign, capacity_status=capacity.status,
                physical_capacity=capacity.objective_value * 1000 if capacity.status == 'optimal' else None)
            proposals.append(proposal)
            if capacity.status != 'optimal' or capacity.objective_value * sign * 1000 < .1 - 1e-7:
                proposal['assessment'] = 'no_directional_library_capacity'
                continue
            with search:
                if sign > 0:
                    reaction.lower_bound = max(reaction.lower_bound, .1 / 1000)
                else:
                    reaction.upper_bound = min(reaction.upper_bound, -.1 / 1000)
                set_weighted_objective(search, info.get('gapfill_candidate_costs', {}))
                sol = search.optimize()
                proposal.update(status=sol.status, cost=sol.objective_value)
                if sol.status != 'optimal':
                    proposal['assessment'] = 'weighted_search_failed'
                    continue
                selected = sorted(r for r in candidates if abs(sol.fluxes[r]) > 1e-6)
                extra = tuple(r for r in selected if r not in existing)
                gate = chemical_gate(native, search, extra)
                proposal.update(selected=selected, additions=list(extra), assessment=gate)
                rows = [dict(reaction=r.id, search_flux=float(sol.fluxes[r.id]),
                    candidate=r.id in candidates, equation=r.reaction, gpr=r.gene_reaction_rule)
                    for r in search.reactions if abs(sol.fluxes[r.id]) > 1e-9]
                pd.DataFrame(rows).to_csv(folder/('selection_%03d_%+d.csv' % (i, sign)), index=False)
                if gate == 'eligible':
                    key = hashlib.sha256('|'.join(extra).encode()).hexdigest()[:10]
                    plans.setdefault(key, dict(additions=list(extra), targets=[]))['targets'].append(dict(reaction=rid, sign=sign))
    report['proposals'] = proposals
    accepted = {}
    for key, plan in sorted(plans.items()):
        trial = native.copy()
        trial.add_reactions([search.reactions.get_by_id(r).copy() for r in plan['additions']])
        plan['safety'] = energy_checks(trial, biomass)
        targets, evidence = safe_targets(trial, plan['targets'], exchanges)
        plan['target_checks'] = evidence
        plan['accepted_targets'] = targets
        plan['accepted'] = plan['safety']['passed'] and bool(targets)
        if plan['accepted']:
            accepted[key] = plan
    report['plans'] = plans
    union = sorted({r for p in accepted.values() for r in p['additions']})
    if not union:
        return finish('no_safe_native_coverage_additions')
    trial = native.copy()
    trial.add_reactions([search.reactions.get_by_id(r).copy() for r in union])
    report['union_safety'] = energy_checks(trial, biomass)
    targets = [t for p in accepted.values() for t in p['accepted_targets']]
    safe, evidence = safe_targets(trial, targets, exchanges)
    report['union_target_checks'] = evidence
    if not report['union_safety']['passed'] or len(safe) != len(targets):
        return finish('rejected_unsafe_union')
    trial.objective = biomass
    trial.objective_direction = 'max'
    growth = trial.optimize()
    report['union_recorded_medium_growth'] = growth.objective_value
    if growth.status != 'optimal' or not math.isfinite(growth.objective_value) or growth.objective_value < baseline.objective_value - 1e-7:
        return finish('rejected_recorded_medium_growth_loss')
    # Diagnostic exchange changes must never enter the actual insertion.
    report['selected_candidate_bounds'] = {r: list(search.reactions.get_by_id(r).bounds) for r in union}
    report['selected_candidate_equations'] = {r: search.reactions.get_by_id(r).reaction for r in union}
    cobra.io.save_json_model(trial, str(folder/'precommit_physical_model.json'))
    return finish('safe_union_selected_before_biomass_finalization', union)
