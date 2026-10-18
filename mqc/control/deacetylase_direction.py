"""Narrow implementation of MQC's existing acetate-production rule.

Only EC 3.5.1.47 with the complete, balanced single-compartment chemistry is
recognized. This is a physiological hydrolysis-direction policy, not a claim
that every acetate-consuming reaction is invalid or an enzyme is present.
Sources: https://enzyme.expasy.org/EC/3.5.1.47
         https://doi.org/10.1111/mmi.12467
"""
import math

ACETATE_RULE = ('The reaction containing acetate is in the direction of producing '
                'acetate except for reactions with another high energy metabolite '
                'such as ATP, Acyl-CoA')
TEMPLATE = [('C9H15N2O5', -1, -1), ('H2O', 0, -1),
            ('C7H14N2O4', 0, 1), ('C2H3O2', -1, 1)]


def deacetylase_direction(reaction):
    """Match enzyme identity AND chemistry, independently of reaction/met IDs."""
    ecs = reaction.annotation.get('ec-code', [])
    if isinstance(ecs, str):
        ecs = [ecs]
    if not isinstance(ecs, (list, tuple, set)) or set(ecs) != {'3.5.1.47'}:
        return None
    if any(word in (reaction.name or '').lower() for word in
           ['overall', 'lumped', 'coupled']):
        return None
    if len(reaction.metabolites) != len(TEMPLATE):
        return None
    compartments = {m.compartment for m in reaction.metabolites}
    if len(compartments) != 1 or not next(iter(compartments)):
        return None
    ratios = []
    for formula, charge, coefficient in TEMPLATE:
        matches = [m for m in reaction.metabolites
                   if m.elements == _elements(formula) and m.charge == charge]
        if len(matches) != 1:
            return None
        ratios.append(reaction.metabolites[matches[0]] / coefficient)
    if any(not math.isfinite(v) or v == 0 for v in ratios):
        return None
    if not all(math.isclose(v, ratios[0], rel_tol=1e-9, abs_tol=1e-12) for v in ratios):
        return None
    return 1 if ratios[0] > 0 else -1


def _elements(formula):
    # Accept formula element order changes, not missing or guessed chemistry.
    return {'C9H15N2O5': {'C':9, 'H':15, 'N':2, 'O':5},
            'H2O': {'H':2, 'O':1}, 'C7H14N2O4': {'C':7, 'H':14, 'N':2, 'O':4},
            'C2H3O2': {'C':2, 'H':3, 'O':2}}[formula]


def deacetylase_bounds(reaction):
    direction = deacetylase_direction(reaction)
    if direction is None:
        return reaction.bounds
    lb, ub = reaction.bounds
    if (direction > 0 and ub < 0) or (direction < 0 and lb > 0):
        raise ValueError('Forced flux conflicts with deacetylase hydrolysis: ' + reaction.id)
    return (max(0, lb), ub) if direction > 0 else (lb, min(0, ub))
