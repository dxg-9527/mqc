"""Direction of the exact quinol-driven 2 NO -> N2O reductase template.

This rule does not infer enzyme presence or correct a cytochrome/quinol EC
annotation mismatch. It only prevents using N2O oxidation to reduce quinone
in the uncoupled five-metabolite NO-reductase lump. Extended energy-coupled
reactions and N2O -> N2 reductases are deliberately outside its scope.

Sources: https://www.rhea-db.org/rhea/30215
         https://www.rhea-db.org/rhea/30211
"""
import math
import re

NO_RULE = 'Quinol-dependent nitric oxide reduction proceeds toward nitrous oxide'


def nitric_oxide_direction(reaction):
    if len(reaction.metabolites) != 5:
        return None
    codes = reaction.annotation.get('ec-code', [])
    codes = [codes] if isinstance(codes, str) else codes
    if not isinstance(codes, (list, tuple)):
        return None
    # Both annotations describe NO reduction; the cofactor mismatch, when
    # present, is not repaired or claimed to have been resolved by this rule.
    if not set(codes) & {'1.7.5.2', '1.7.2.5'} and (reaction.name or '').lower() not in {
            'nitric oxide reductase', 'nitric oxide reductase (quinol)',
            'nitrous-oxide:(acceptor) oxidoreductase (no-forming)'}:
        return None
    small = {'no': ({'no', 'C00533'}, {'nitric oxide'}, {'N': 1, 'O': 1}, -2),
             'n2o': ({'n2o', 'C00887'}, {'nitrous oxide'}, {'N': 2, 'O': 1}, 1),
             'water': ({'h2o', 'cpd00001', 'C00001'}, {'water', 'h2o'}, {'H': 2, 'O': 1}, 1)}
    found, other = {}, []
    for met, coefficient in reaction.metabolites.items():
        if not math.isfinite(coefficient) or coefficient == 0 or met.charge != 0:
            return None
        base = re.sub(r'(?:_[a-z]\d*|\[[a-z]\d*\]|@[a-z]\d*)$', '', met.id)
        identities = {base}
        for key in ('bigg.metabolite', 'seed.compound', 'kegg.compound'):
            value = met.annotation.get(key, [])
            if isinstance(value, str):
                identities.add(value)
            elif isinstance(value, (list, tuple)):
                identities.update(v for v in value if isinstance(v, str))
        roles = [key for key, (ids, names, formula, expected) in small.items()
                 if (identities & ids or (met.name or '').lower() in names) and met.elements == formula]
        if len(roles) == 1 and roles[0] not in found:
            found[roles[0]] = (met, coefficient)
        else:
            other.append((met, coefficient, base))
    if len(found) != 3 or len(other) != 2:
        return None
    if len({met.compartment for met, _ in found.values()}) != 1:
        return None
    scale = found['n2o'][1]
    if any(not math.isclose(found[key][1] / scale, data[3], abs_tol=1e-9) for key, data in small.items()):
        return None
    oxidized = [v for v in other if math.isclose(v[1]/scale, 1, abs_tol=1e-9)]
    reduced = [v for v in other if math.isclose(v[1]/scale, -1, abs_tol=1e-9)]
    if len(oxidized) != 1 or len(reduced) != 1:
        return None
    ox, _, ox_id = oxidized[0]
    red, _, red_id = reduced[0]
    pairs = {('q8','q8h2'), ('mqn8','mql8'), ('dmmq8','dmmql8')}
    if (ox_id, red_id) not in pairs or ox.compartment != red.compartment:
        return None
    expected = dict(ox.elements)
    expected['H'] = expected.get('H', 0) + 2
    if red.elements != expected or reaction.check_mass_balance():
        return None
    return 1 if scale > 0 else -1


def nitric_oxide_bounds(reaction):
    direction = nitric_oxide_direction(reaction)
    if direction is None:
        return reaction.bounds
    lb, ub = reaction.bounds
    if (direction > 0 and ub < 0) or (direction < 0 and lb > 0):
        raise ValueError('Forced flux conflicts with nitric oxide reduction: ' + reaction.id)
    return (max(0, lb), ub) if direction > 0 else (lb, min(0, ub))
