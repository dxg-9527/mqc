"""Synthesis direction for two explicitly identified legacy BioB templates.

This prevents using biotin biosynthesis backwards as a reductant generator.
It does NOT claim to repair the omitted radical-SAM/FeS chemistry of the old
templates, and never adds a cofactor, reaction, substrate or phenotype task.
Sources: https://iubmb.qmul.ac.uk/enzyme/EC2/8/1/6.html
         https://www.rhea-db.org/rhea/22061
"""
import math

BIOTIN_RULE = 'Legacy biotin synthase templates proceed toward biotin synthesis'
COMMON = [('C10H17N2O3', -1, -1), ('C10H15N2O3S', -1, 1)]
TEMPLATES = [COMMON + [('C3H7NO2S', 0, -1), ('C3H7NO2', 0, 1), ('H', 1, 2)],
             COMMON + [('S', 0, -2), ('H2S', 0, 1)]]


def synthesis_direction(reaction):
    """Require identity, exact ionic chemistry/ratios, and one compartment."""
    ecs = reaction.annotation.get('ec-code', [])
    ecs = [ecs] if isinstance(ecs, str) else ecs
    if ecs and set(ecs) != {'2.8.1.6'}:
        return None
    name = (reaction.name or '').strip().lower()
    if not ecs and name not in {'biotin synthase', 'biotin synthase (ala-l producing)'}:
        return None
    compartments = {m.compartment for m in reaction.metabolites}
    if len(compartments) != 1 or not next(iter(compartments)):
        return None
    for specs in TEMPLATES:
        if len(reaction.metabolites) != len(specs):
            continue
        ratios = []
        for formula, charge, coefficient in specs:
            hits = [m for m in reaction.metabolites if m.formula == formula and m.charge == charge]
            if len(hits) != 1:
                break
            ratios.append(reaction.metabolites[hits[0]] / coefficient)
        if len(ratios) != len(specs) or any(not math.isfinite(x) or x == 0 for x in ratios):
            continue
        if all(math.isclose(x, ratios[0], rel_tol=1e-9, abs_tol=1e-12) for x in ratios):
            return 1 if ratios[0] > 0 else -1
    return None


def biotin_bounds(reaction):
    direction = synthesis_direction(reaction)
    if direction is None:
        return reaction.bounds
    lb, ub = reaction.bounds
    if (direction > 0 and ub < 0) or (direction < 0 and lb > 0):
        raise ValueError('Forced flux conflicts with biotin synthesis direction: ' + reaction.id)
    return (max(0, lb), ub) if direction > 0 else (lb, min(0, ub))
