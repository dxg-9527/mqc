"""Native gap-filling may repair pathways, not invent membrane compartments.

Use the input target inventory; retain legitimate compartments already present
in eukaryotic models. C_ is the SBML namespace prefix, not an extra membrane.
This guard changes neither reaction chemistry nor experimental tasks.
"""


def compartment_key(value):
    value = (value or '').strip()
    return value[2:] if value.startswith('C_') else value


def target_inventory(model):
    return sorted({compartment_key(m.compartment) for m in model.metabolites
                   if compartment_key(m.compartment)})


def foreign_compartments(reaction, allowed, resolved=None):
    invalid = {}
    for original in reaction.metabolites:
        met = (resolved or {}).get(original.id, original)
        key = compartment_key(met.compartment)
        if not key or key not in allowed:
            invalid[met.id] = met.compartment or ''
    if not invalid:
        return None
    return dict(reaction=reaction.id, invalid_metabolite_compartments=invalid,
                reason='candidate_requires_compartment_absent_from_input_target')
