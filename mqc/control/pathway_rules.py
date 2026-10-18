"""Narrow chemistry checks for original MQC templates and growth-path rules.

No strain IDs, fitted flux bounds, phenotype targets, or new metabolic routes.
Direction interventions remain explicit trial policies, not universal biology.
"""
from mqc.defaults import C_COMPARTMENT
import math

PROTON_RULE = 'The direction of the reaction for proton-driven generation of ATP is incorrect'
HYDROLYSIS_RULE = 'Polysaccharide hydrolysis reaction is the direction of hydrolysis'
SUGAR_IDENTITY_RULE = 'Trehalose-phosphate hydrolysis has a conflicting monosaccharide identity'


def sugar_identity_conflict(reaction, resolved=None):
    """Detect one documented, atom-balanced but structurally incoherent template.

Trehalose-6P + water -> glucose-6P + fructose is neither TRE6PH nor
sucrose-6P hydrolysis. Do NOT guess whether its substrate, product, or enzyme
label was wrong. Recognize chemical identities, not strain/reaction names.
Formula-only matches, ambiguous annotations and coupled reactions are misses.
    EC 3.2.1.93: https://enzyme.expasy.org/EC/3.2.1.93
    """
    if len(reaction.metabolites) != 4:
        return None
    name = (reaction.name or '').lower()
    ec = reaction.annotation.get('ec-code', [])
    ec = {ec} if isinstance(ec, str) else set(ec) if isinstance(ec, (list, tuple, set)) else set()
    # The same NET equation could describe hydrolysis plus isomerization.
    # Only flag a record that actually claims a single supported hydrolase.
    if not (ec & {'3.2.1.26', '3.2.1.93'} or any(label in name for label in
            ['fructofuranosidase', 'trehalose-6-phosphate hydrolase', 'phosphotrehalase'])):
        return None
    if any(label in name for label in ['overall', 'lumped', 'isomerase', 'isomerization']):
        return None
    aliases = {
        'water': {'h2o', 'cpd00001', 'C00001'},
        'tre6p': {'tre6p', 'cpd00523', 'C00689'},
        'fructose': {'fru', 'cpd00082', 'C00095'},
        'g6p': {'g6p', 'cpd00079', 'C00092'},
    }
    chemistry = {'water': ('H2O', 0), 'tre6p': ('C12H21O14P', -2),
                 'fructose': ('C6H12O6', 0), 'g6p': ('C6H11O9P', -2)}
    # Recognize competing identities as conflicts rather than trusting one
    # namespace and silently discarding a contradictory explicit annotation.
    competing = {'glc__D', 'cpd00027', 'C00031', 'suc6p', 'cpd00514', 'C01083'}
    found, compartments = {}, set()
    for original, coefficient in reaction.metabolites.items():
        met = (resolved or {}).get(original.id, original)
        comp = met.compartment or ''
        if not comp or not math.isfinite(coefficient):
            return None
        compartments.add(comp)
        base = met.id
        for suffix_comp in [comp] + ([comp[2:]] if comp.startswith('C_') else []):
            for suffix in ('_' + suffix_comp, '[' + suffix_comp + ']', '@' + suffix_comp):
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
        identities = {base}
        for key in ['bigg.metabolite', 'seed.compound', 'kegg.compound']:
            value = met.annotation.get(key, [])
            if isinstance(value, str):
                identities.add(value)
            elif isinstance(value, (list, tuple, set)):
                identities.update(v for v in value if isinstance(v, str))
        roles = [role for role, names in aliases.items() if identities & names]
        if len(roles) != 1 or identities & competing:
            return None
        role = roles[0]
        if role in found or (met.formula, met.charge) != chemistry[role]:
            return None
        found[role] = coefficient
    if len(compartments) != 1 or set(found) != set(aliases):
        return None
    scale = -found['water']
    if abs(scale) < 1e-12 or any(not math.isclose(found[role] / scale, expected,
            rel_tol=1e-9, abs_tol=1e-9) for role, expected in
            {'water': -1, 'tre6p': -1, 'fructose': 1, 'g6p': 1}.items()):
        return None
    return dict(reaction=reaction.id, rule=SUGAR_IDENTITY_RULE,
        observed='trehalose-6P + water <=> fructose + glucose-6P',
        expected_hydrolysis_products='glucose + glucose-6P',
        source='https://enzyme.expasy.org/EC/3.2.1.93',
        automatic_stoichiometry_replacement=False,
        ambiguity='substrate/product/enzyme identity must be reconciled; no enzyme presence inferred')


def verified_chemistry(r):
    return all(m.formula and m.charge is not None for m in r.metabolites) and not r.check_mass_balance()


def atp_coupling(r):
    if len(r.metabolites) != 6 or not verified_chemistry(r):
        return None
    # Exact chemical species, not reaction-ID or name substring matching.
    specs = {'ATP': ({'C':10,'H':12,'N':5,'O':13,'P':3}, -4),
             'ADP': ({'C':10,'H':12,'N':5,'O':10,'P':2}, -3),
             'Pi': ({'H':1,'O':4,'P':1}, -2), 'water': ({'H':2,'O':1}, 0)}
    found = {}
    for name, (elements, charge) in specs.items():
        hits = [m for m in r.metabolites if m.elements == elements and m.charge == charge]
        if len(hits) != 1: return None
        found[name] = hits[0]
    comp = found['ATP'].compartment
    if comp not in C_COMPARTMENT or any(m.compartment != comp for m in found.values()): return None
    protons = [m for m in r.metabolites if m.elements == {'H':1} and m.charge == 1]
    inner = [m for m in protons if m.compartment == comp]
    outer = [m for m in protons if m.compartment and m.compartment != comp]
    if len(inner) != 1 or len(outer) != 1: return None
    scale = r.metabolites[found['ATP']]
    for name, coeff in {'ATP':1, 'ADP':-1, 'Pi':-1, 'water':1}.items():
        if abs(r.metabolites[found[name]] / scale - coeff) > 1e-9: return None
    n = -r.metabolites[outer[0]] / scale
    if n <= 1 or abs(r.metabolites[inner[0]] / scale - (n-1)) > 1e-9: return None
    return dict(synthesis_direction=1 if scale > 0 else -1, protons_per_ATP=n,
                inner=inner[0].id, outer=outer[0].id)


def prepare_new_atpases(model, check_model, before_ids, info, enabled=False):
    for r in model.reactions:
        if r.id in before_ids or r.id not in {'add_ATPS4_e', 'add_ATPS4rpp_p'}: continue
        signature = atp_coupling(r)
        r.notes['MQC_origin'] = 'automatic_ATP_proton_template'
        if not enabled: continue
        if signature is None or r.lower_bound > 0 or r.upper_bound < 0:
            raise ValueError('Unverified auto ATPase template: ' + r.id)
        before = r.bounds
        r.bounds = (-1000, 1000)
        check_model.reactions.get_by_id(r.id).bounds = r.bounds
        r.notes['MQC_direction_policy'] = 'trial_ATP_consuming_proton_export_allowed'
        info.setdefault('auto_atpase_direction_changes', []).append(
            dict(reaction=r.id, before=before, after=r.bounds, **signature))


def disaccharide_direction(r):
    if len(r.metabolites) not in (3, 4) or not verified_chemistry(r): return None
    if len({m.compartment for m in r.metabolites}) != 1: return None
    water = [m for m in r.metabolites if m.elements == {'H':2,'O':1} and m.charge == 0]
    di = [m for m in r.metabolites if m.elements == {'C':12,'H':22,'O':11} and m.charge == 0]
    if len(water) != 1 or len(di) != 1: return None
    scale = -r.metabolites[water[0]]
    if abs(r.metabolites[di[0]] / scale + 1) > 1e-9: return None
    mono = [m for m in r.metabolites if m not in water + di]
    if any(m.elements != {'C':6,'H':12,'O':6} or m.charge != 0 or r.metabolites[m] / scale <= 0 for m in mono): return None
    if abs(sum(r.metabolites[m] / scale for m in mono) - 2) > 1e-9: return None
    return 1 if scale > 0 else -1


def apply_disaccharide_policy(model, info, enabled=False):
    if not enabled: return
    for record in info['reactions']:
        if HYDROLYSIS_RULE not in record.get('rules', {}): continue
        r = model.reactions.get_by_id(record['id']); direction = disaccharide_direction(r)
        if direction is None: continue
        before = r.bounds
        if direction == 1 and r.upper_bound >= 0: r.bounds = (max(0, r.lower_bound), r.upper_bound)
        elif direction == -1 and r.lower_bound <= 0: r.bounds = (r.lower_bound, min(0, r.upper_bound))
        else: raise ValueError('Forced flux conflicts with hydrolysis trial: ' + r.id)
        if r.bounds != before:
            record['bounds'] = r.bounds
            info.setdefault('growth_hydrolysis_direction_changes', []).append(
                dict(reaction=r.id, before=before, after=r.bounds, rule=HYDROLYSIS_RULE,
                     scope='conditional_hydrolysis_direction_trial_before_original_biomass_gapfill'))


def growth_rule_audit(model, info, fluxes):
    rows = []
    for record in info['reactions']:
        if record['id'] not in model.reactions: continue
        r = model.reactions.get_by_id(record['id']); v = float(fluxes[r.id])
        if abs(v) < 1e-8: continue
        direction = disaccharide_direction(r)
        if direction is not None and v * direction < -1e-8:
            rows.append(dict(reaction=r.id, flux=v, equation=r.reaction,
                existing_rule_hit=HYDROLYSIS_RULE in record.get('rules', {}),
                penalty=record.get('net_penalty_points'),
                status='reverse_disaccharide_hydrolysis_requires_review'))
    return rows
