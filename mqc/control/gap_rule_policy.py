"""Apply a narrow subset of MQC rules throughout its original biomass repair.

No phenotype tasks, experimental rates, substrate list, or strain names enter
this policy. Native S*1000 search semantics and primary objectives are retained.
This is an isolated algorithm candidate, not a universal thermodynamic claim.
"""
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.core.solution import get_solution

from mqc.control.gapfill_chemistry import assess_gapfill_chemistry
from mqc.control.pathway_rules import disaccharide_direction, HYDROLYSIS_RULE
from mqc.control.biotin_direction import biotin_bounds, BIOTIN_RULE
from mqc.control.deacetylase_direction import deacetylase_bounds, ACETATE_RULE
from mqc.control.desulfurase_direction import desulfurase_bounds, DESULFURASE_RULE
from mqc.control.quinone_methylation import methylation_bounds, METHYL_REDOX_RULE


def enabled(cfg):
    return getattr(cfg, 'biomass_rule_consistency', True)


def hydrolysis_bounds(reaction):
    """Intersect, never enlarge, existing bounds with verified hydrolysis.

The stoichiometric signature excludes phosphorolysis, PTS, ATP-dependent
glycosylation and membrane-coupled reactions. Unknown chemistry is not a hit.
"""
    direction = disaccharide_direction(reaction)
    if direction is None:
        return reaction.bounds
    lb, ub = reaction.bounds
    if (direction > 0 and ub < 0) or (direction < 0 and lb > 0):
        raise ValueError('Forced flux conflicts with MQC hydrolysis rule: ' + reaction.id)
    return (max(lb, 0), ub) if direction > 0 else (lb, min(ub, 0))


def apply_native_rules(model, info):
    """Enforce the chemically verified subset of existing direction rules."""
    lookup = {r['id']: r for r in info['reactions']}
    changes = info.setdefault('biomass_rule_changes', [])
    for reaction in model.reactions:
        record = lookup.get(reaction.id, {})
        for rule, policy in [(HYDROLYSIS_RULE, hydrolysis_bounds),
                             (BIOTIN_RULE, biotin_bounds), (ACETATE_RULE, deacetylase_bounds),
                             (DESULFURASE_RULE, desulfurase_bounds),
                             (METHYL_REDOX_RULE, methylation_bounds)]:
            if rule not in record.get('rules', {}):
                continue
            old = reaction.bounds
            reaction.bounds = policy(reaction)
            if reaction.bounds != old:
                record['bounds'] = list(reaction.bounds)
                changes.append(dict(reaction=reaction.id, before=list(old),
                                    after=list(reaction.bounds), rule=rule,
                                    stage='before_original_biomass_precursor_search'))


def configure_candidates(model, candidate_ids, info):
    """Use resolved target chemistry, and retain policies for actual insertion.

Known imbalance costs +3, matching MQC's existing mass/charge penalty. The
base cost is 1; unknown charge is recorded but not equated to known imbalance.
This changes ranking, not a claim that a review item is chemically valid.
"""
    bounds, costs, records = {}, {}, []
    for rid in candidate_ids:
        r = model.reactions.get_by_id(rid)
        old = r.bounds
        r.bounds = hydrolysis_bounds(r)
        r.bounds = biotin_bounds(r)
        r.bounds = deacetylase_bounds(r)
        r.bounds = desulfurase_bounds(r)
        r.bounds = methylation_bounds(r)
        assessment = assess_gapfill_chemistry(r)
        known_imbalance = bool(assessment['atom_imbalance']) or (
            assessment['charge_imbalance'] is not None and assessment['charge_imbalance'] != 0)
        cost = 1.0 + (3.0 if known_imbalance else 0.0)
        bounds[rid] = list(r.bounds)
        if cost != 1:
            costs[rid] = cost
        if old != r.bounds or cost != 1:
            records.append(dict(reaction=rid, before=list(old), after=list(r.bounds),
                                search_cost=cost, chemistry=assessment))
    info['gapfill_candidate_bounds'] = bounds
    info['gapfill_candidate_costs'] = costs
    info['gapfill_rule_policy'] = dict(records=records, primary_target_unchanged=True,
        native_search_scale=1000, experimental_values_used=False,
        policy='hydrolysis_biotin_DAP_deacetylase_desulfurase_directions_methyl_redox_quarantine_known_imbalance_penalty_3')


def rule_pfba(model, costs):
    """Original pFBA with explicit per-candidate costs, same optimum first."""
    with model:
        add_pfba(model)
        coefficients = {}
        for rid, cost in costs.items():
            if rid not in model.reactions:
                continue
            if not (cost >= 1 and cost < float('inf')):
                raise ValueError('Invalid gap-filling cost for ' + rid)
            r = model.reactions.get_by_id(rid)
            coefficients[r.forward_variable] = cost
            coefficients[r.reverse_variable] = cost
        model.objective.set_linear_coefficients(coefficients)
        model.slim_optimize(error_value=None)
        return get_solution(model)


def candidate_copy(library, rid, info):
    """Never share a Reaction object or silently reopen searched directions."""
    r = library.reactions.get_by_id(rid).copy()
    if rid in info.get('gapfill_candidate_bounds', {}):
        r.bounds = tuple(info['gapfill_candidate_bounds'][rid])
    return r


def validate_final_rule_bounds(model, info, control):
    """Fail explicitly if later restoration reopened a constrained direction."""
    from mqc.control.postgap_control import PostGapValidationError
    expected = {row['reaction']: row['after'] for row in info.get('biomass_rule_changes', [])}
    expected.update(info.get('gapfill_candidate_bounds', {}))
    violations = []
    for rid, (lb, ub) in expected.items():
        if rid not in model.reactions:
            continue
        r = model.reactions.get_by_id(rid)
        if r.lower_bound < lb-1e-9 or r.upper_bound > ub+1e-9:
            violations.append(dict(reaction=rid, expected=[lb,ub], actual=list(r.bounds)))
    report = dict(status='failed' if violations else 'passed', violations=violations,
                  scope='directions_applied_by_biomass_rule_policy')
    info['biomass_rule_validation'] = control['biomass_rule_validation'] = report
    if violations:
        control['final_model_status'] = 'biomass_rule_validation_failed'
        raise PostGapValidationError('Biomass rule directions reopened after search')
