"""Physiological direction of the legacy free-sulfur desulfurase lump.

The exact cysteine -> alanine + S template abstracts sulfur mobilization;
it is not evidence for making cysteine from imported elemental sulfur.
Do not apply this policy to protein-bound persulfide transfer, sulfide-based
cysteine synthesis, desulfhydrases, or arbitrary sulfur reactions. This is a
narrow direction rule, not a thermodynamic proof of microscopic irreversibility
and not a reconstruction of the missing sulfur-carrier chemistry.

Sources: https://iubmb.qmul.ac.uk/enzyme/EC2/8/1/7.html
         https://pubmed.ncbi.nlm.nih.gov/9278392/
"""
import math
import re

DESULFURASE_RULE = 'Legacy free-sulfur cysteine desulfurase proceeds toward sulfur mobilization'
IDENTITIES = {
    'cysteine': {'cys__L', 'C00097', 'cpd00084'},
    'alanine': {'ala__L', 'C00041', 'cpd00035'},
    'sulfur': {'s', 'C00087'},
}
NAMES = {'cysteine': {'l-cysteine', 'l cysteine'},
         'alanine': {'l-alanine', 'l alanine'},
         'sulfur': {'sulfur', 'elemental sulfur'}}
CHEMISTRY = {'cysteine': ('C3H7NO2S', 0, -1),
             'alanine': ('C3H7NO2', 0, 1), 'sulfur': ('S', 0, 1)}


def desulfuration_direction(reaction):
    """Require enzyme identity, L-species identities, chemistry and exact ratio.

    Matching works with namespace IDs, annotations or exact chemical names;
    neither reaction ID nor strain identity controls the result. Opaque or
    contradictory identities and incomplete/extended equations are misses.
    """
    if len(reaction.metabolites) != 3:
        return None
    ecs = set()
    for key, value in reaction.annotation.items():
        if key.lower() not in {'ec-code', 'ec number', 'ec'}:
            continue
        values = [value] if isinstance(value, str) else value if isinstance(value, (list, tuple, set)) else []
        for entry in values:
            if isinstance(entry, str):
                ecs.update(re.findall(r'\b\d+\.\d+\.\d+\.\d+\b', entry))
    if ecs and ecs != {'2.8.1.7'}:
        return None
    if not ecs and (reaction.name or '').strip().lower() not in {
            'cysteine desulfurase', 'l-cysteine desulfurase'}:
        return None
    compartments = {m.compartment for m in reaction.metabolites}
    if len(compartments) != 1 or not next(iter(compartments)):
        return None
    found = {}
    for met, coeff in reaction.metabolites.items():
        if not math.isfinite(coeff) or coeff == 0:
            return None
        base = met.id
        comp = met.compartment
        for c in [comp] + ([comp[2:]] if comp.startswith('C_') else []):
            for suffix in ['_'+c, '['+c+']', '@'+c]:
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
        ids = {base}
        for key, value in met.annotation.items():
            if key.lower() not in {'bigg.metabolite', 'seed.compound', 'kegg.compound'}:
                continue
            values = [value] if isinstance(value, str) else value if isinstance(value, (list,tuple,set)) else []
            ids.update(v.rstrip('/').rsplit('/', 1)[-1] for v in values if isinstance(v, str))
        name = (met.name or '').strip().lower()
        # Explicit D/beta names are not rescued by a conflicting L alias.
        if name.startswith(('d-', 'd ', 'beta-', 'beta ', 'β')) or ids & {
                'ala__D', 'ala_B', 'cys__D', 'C00133', 'C00099', 'C00793'}:
            return None
        roles = [role for role in IDENTITIES if ids & IDENTITIES[role] or name in NAMES[role]]
        if len(roles) != 1 or roles[0] in found:
            return None
        role = roles[0]
        formula, charge, expected = CHEMISTRY[role]
        if (met.formula, met.charge) != (formula, charge):
            return None
        found[role] = coeff / expected
    if set(found) != set(CHEMISTRY):
        return None
    scale = found['cysteine']
    if not all(math.isclose(v, scale, rel_tol=1e-9, abs_tol=1e-12) for v in found.values()):
        return None
    return 1 if scale > 0 else -1


def desulfurase_bounds(reaction):
    direction = desulfuration_direction(reaction)
    if direction is None:
        return reaction.bounds
    lb, ub = reaction.bounds
    if (direction > 0 and ub < 0) or (direction < 0 and lb > 0):
        raise ValueError('Forced flux conflicts with legacy desulfurase direction: ' + reaction.id)
    return (max(0, lb), ub) if direction > 0 else (lb, min(0, ub))
