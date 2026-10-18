"""Reconcile an explicit flavodoxin annotation with an existing carrier pool.

RHEA:44140 uses protein-bound flavodoxin, not the freely exchangeable FMN pool.
Only reconcile its exact pyruvate reaction when the same model already has a
balanced, equivalent bound-carrier reaction. No new function, pool or bound.
Source: https://www.rhea-db.org/rhea/44140
"""
import json
import math

SOURCE = 'https://www.rhea-db.org/rhea/44140'
RHEA = {'44140', '44141', '44142', '44143'}
ALIASES = {
    'pyr': {'pyr', 'cpd00020', 'C00022'},
    'coa': {'coa', 'cpd00010', 'C00010'},
    'accoa': {'accoa', 'cpd00022', 'C00024'},
    'co2': {'co2', 'cpd00011', 'C00011'},
    'h': {'h', 'cpd00067', 'C00080'},
    'free_ox': {'fmn', 'cpd00050', 'C00061', 'META:FMN'},
    'free_red': {'fmnh2', 'cpd01270', 'C01847', 'META:FMNH2'},
    'bound_ox': {'fldox', 'cpd12207', 'C02869', 'META:Oxidized-flavodoxins'},
    'bound_red': {'fldrd', 'cpd12173', 'C02745', 'META:Reduced-flavodoxins'},
}
CORE = {'pyr': ({'C':3,'H':3,'O':3}, -1),
        'coa': ({'C':21,'H':32,'N':7,'O':16,'P':3,'S':1}, -4),
        'accoa': ({'C':23,'H':34,'N':7,'O':17,'P':3,'S':1}, -4),
        'co2': ({'C':1,'O':2}, 0), 'h': ({'H':1}, 1)}
EXPECTED = {'pyr':-1, 'coa':-1, 'h':-1, 'free_ox':-1,
            'accoa':1, 'co2':1, 'free_red':1}


def role(met):
    comp = met.compartment or ''
    if not comp:
        return None
    base = met.id
    for c in [comp] + ([comp[2:]] if comp.startswith('C_') else []):
        for suffix in ['_'+c, '['+c+']', '@'+c]:
            if base.endswith(suffix):
                base = base[:-len(suffix)]
                break
    ids = {base}
    for key in ['bigg.metabolite', 'seed.compound', 'kegg.compound', 'biocyc']:
        values = met.annotation.get(key, [])
        if isinstance(values, str):
            values = [values]
        if not isinstance(values, (list, tuple, set)) or not all(isinstance(v,str) for v in values):
            return None
        ids.update(values)
    hits = [key for key, aliases in ALIASES.items() if ids & aliases]
    return hits[0] if len(hits) == 1 else None


def two_hydrogen_pair(oxidized, reduced):
    if not oxidized.formula or not reduced.formula:
        return False
    if oxidized.charge is None or reduced.charge != oxidized.charge:
        return False
    if not math.isfinite(oxidized.charge):
        return False
    delta = dict(reduced.elements)
    for key,value in oxidized.elements.items():
        delta[key] = delta.get(key,0)-value
    return {k:v for k,v in delta.items() if v} == {'H':2}


def proposal(reaction, model):
    ids = reaction.annotation.get('rhea', [])
    ids = [ids] if isinstance(ids, str) else ids
    if not isinstance(ids,(list,tuple,set)) or not ids or not set(ids) <= RHEA:
        return None
    if len(reaction.metabolites) != 7 or any(word in (reaction.name or '').lower()
            for word in ['lumped', 'overall', 'coupled']):
        return None
    found = {}
    for met,value in reaction.metabolites.items():
        key = role(met)
        if key not in EXPECTED or key in found or not math.isfinite(value):
            return None
        found[key] = met
    if set(found) != set(EXPECTED):
        return None
    compartments = {m.compartment for m in found.values()}
    if len(compartments) != 1:
        return None
    comp = next(iter(compartments))
    scale = -reaction.metabolites[found['pyr']]
    if abs(scale) < 1e-12 or any(not math.isclose(reaction.metabolites[found[k]]/scale,c,
            rel_tol=1e-9,abs_tol=1e-9) for k,c in EXPECTED.items()):
        return None
    if any(found[k].elements != formula or found[k].charge != charge
           for k,(formula,charge) in CORE.items()):
        return None
    ox,red = found['free_ox'],found['free_red']
    base = dict(ox.elements);hydrogens = base.pop('H',None)
    if base != {'C':17,'N':4,'O':9,'P':1} or (hydrogens,ox.charge) not in {(21,0),(19,-2),(18,-3)}:
        return None
    if not two_hydrogen_pair(ox,red) or reaction.check_mass_balance():
        return None
    pool = {}
    for key in ['bound_ox','bound_red']:
        hits = [m for m in model.metabolites if m.compartment == comp and role(m) == key]
        if len(hits) != 1:
            return None
        pool[key] = hits[0]
    if not two_hydrogen_pair(pool['bound_ox'],pool['bound_red']):
        return None
    new = {m:c for m,c in reaction.metabolites.items() if m not in [ox,red]}
    new[pool['bound_ox']] = -scale;new[pool['bound_red']] = scale
    # Require a pre-existing, balanced counterpart: never invent a route.
    peers = pool['bound_ox'].reactions & pool['bound_red'].reactions
    counterparts = []
    for peer in peers:
        if set(peer.metabolites) != set(new) or peer.check_mass_balance():
            continue
        factor = peer.metabolites[found['pyr']]/new[found['pyr']]
        if abs(factor)>1e-12 and all(math.isclose(peer.metabolites[m],c*factor,
                rel_tol=1e-9,abs_tol=1e-9) for m,c in new.items()):
            counterparts.append(peer.id)
    if not counterparts:
        return None
    return new, dict(reaction=reaction.id,source=SOURCE,compartment=comp,
        existing_counterparts=sorted(counterparts),
        before={m.id:c for m,c in reaction.metabolites.items()},
        after={m.id:c for m,c in new.items()},
        reason='declared_protein_bound_carrier_mapped_to_free_FMN_despite_existing_separate_pool')


def repair_model(model):
    records = []
    for reaction in model.reactions:
        plan = proposal(reaction,model)
        if plan is None:
            continue
        new,record = plan
        original_annotation = dict(reaction.annotation)
        delta = {m:-c for m,c in reaction.metabolites.items()}
        for met,value in new.items():
            delta[met] = delta.get(met,0)+value
        reaction.add_metabolites(delta)
        # Retain superseded mappings in provenance, not as current identities.
        reaction.notes['MQC_ORIGINAL_CARRIER_ANNOTATION'] = json.dumps(original_annotation)
        reaction.notes['MQC_FLAVODOXIN_IDENTITY_REPAIR'] = json.dumps(record)
        reaction.annotation = {k:v for k,v in original_annotation.items()
                               if k not in ['seed.reaction','metanetx.reaction','bigg.reaction']}
        records.append(record)
    return records
