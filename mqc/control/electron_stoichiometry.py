"""Repair exact one-carrier encodings of documented two-electron reactions.

Only an existing single-electron Fe2S2 pair (+1/+2) with fully identified
substrates, full atom balance and the predicted one-charge residual qualifies.
No enzymes, donors, directions, metabolites or phenotype tasks are added.
"""
import math

CHEMISTRY = {
    '2mecdp': ('C5H10O9P2', -2), 'h2mb4p': ('C5H9O8P2', -3),
    'h': ('H', 1), 'h2o': ('H2O', 0), 'h2': ('H2', 0),
    'fd_red': ('Fe2S2', 1), 'fd_ox': ('Fe2S2', 2),
    '3mob': ('C5H7O3', -1), 'coa': ('C21H32N7O16P3S', -4),
    'ibcoa': ('C25H38N7O17P3S', -4), 'co2': ('CO2', 0),
    'nad': ('C21H26N7O14P2', -1), 'nadh': ('C21H27N7O14P2', -2),
    'nadp': ('C21H25N7O17P3', -3), 'nadph': ('C21H26N7O17P3', -4),
}
ALIASES = {'fdxrd': 'fd_red', 'fdxo_2_2': 'fd_ox', 'fdxox': 'fd_ox'}
TEMPLATES = {
    'ISPG': ({'2mecdp': -1, 'fd_red': -1, 'h': -1, 'fd_ox': 1,
              'h2mb4p': 1, 'h2o': 1}, '1/17/7/1'),
    'VOR': ({'3mob': -1, 'coa': -1, 'fd_ox': -1, 'co2': 1,
             'fd_red': 1, 'h': 1, 'ibcoa': 1}, '1/2/7/7'),
    'FNR_NADP': ({'fd_red': -1, 'h': -1, 'nadp': -1, 'fd_ox': 1,
                  'nadph': 1}, '1/18/1/2'),
    'FNR_NAD': ({'fd_red': -1, 'h': -1, 'nad': -1, 'fd_ox': 1,
                 'nadh': 1}, '1/18/1/3'),
    'HYDROGENASE': ({'h2': -1, 'fd_ox': -1, 'fd_red': 1, 'h': 2}, '1/12/7/2'),
}


def identity(met):
    """Require explicit, coherent compartments; no suffix-only guessing."""
    compartment = (met.compartment or '').lower()
    if not compartment:
        return None, None
    name = met.id.lower()
    aliases = [compartment]
    if compartment.startswith('c_'):
        aliases.append(compartment[2:])
    for alias in aliases:
        for suffix in ('_' + alias, '[' + alias + ']', '@' + alias):
            if name.endswith(suffix):
                base = name[:-len(suffix)]
                return ALIASES.get(base, base), compartment
    return None, compartment


def proposal(reaction, resolved=None, families=None):
    """Return only a fully checked correction plan; never mutate on a miss."""
    actual = {}; compartments = set(); residual = {}
    for original, value in reaction.metabolites.items():
        met = (resolved or {}).get(original.id, original)
        base, compartment = identity(met)
        if base not in CHEMISTRY or base in actual or not math.isfinite(value):
            return None
        formula, charge = CHEMISTRY[base]
        if met.formula != formula or met.charge != charge:
            return None
        actual[base] = (original, value)
        compartments.add(compartment)
        for element, count in met.elements.items():
            residual[element] = residual.get(element, 0) + count * value
        residual['charge'] = residual.get('charge', 0) + met.charge * value
    if len(compartments) != 1:
        return None
    for family, (expected, ec) in TEMPLATES.items():
        if families is not None and family not in families:
            continue
        if set(actual) != set(expected):
            continue
        anchor = next(iter(expected))
        scale = actual[anchor][1] / expected[anchor]
        if abs(scale) < 1e-12:
            continue
        if any(not math.isclose(actual[k][1] / scale, c, rel_tol=1e-9, abs_tol=1e-9)
               for k, c in expected.items()):
            continue
        tol = abs(scale) * 1e-8
        if any(abs(v) > tol for k, v in residual.items() if k != 'charge'):
            continue
        correction = actual['fd_ox'][1] * 2 + actual['fd_red'][1]
        if not math.isclose(residual.get('charge', 0), -correction,
                            rel_tol=1e-9, abs_tol=tol):
            continue
        before = {m.id: v for m, v in reaction.metabolites.items()}
        after = dict(before)
        for key in ['fd_red', 'fd_ox']:
            m, value = actual[key]
            after[m.id] = 2 * value
        return dict(reaction=reaction.id, family=family, before=before, after=after,
            charge_before=residual['charge'], charge_after=0,
            source='https://iubmb.qmul.ac.uk/enzyme/EC' + ec + '.html',
            bounds=list(reaction.bounds), gpr=reaction.gene_reaction_rule,
            enzyme_presence_confirmed=False, carrier_identity_changed=False)
    return None


def repair_model(model, resolved=None, families=None):
    plans = [p for r in model.reactions if (p := proposal(r, resolved, families))]
    for p in plans:
        r = model.reactions.get_by_id(p['reaction'])
        r.add_metabolites({m: p['after'][m.id] for m in r.metabolites}, combine=False)
        # Chemistry with a resolved mapping was checked above, not by using
        # potentially different annotations on the unmerged reference objects.
        if resolved is None and r.check_mass_balance():
            r.add_metabolites({m: p['before'][m.id] for m in r.metabolites}, combine=False)
            raise ValueError('Electron correction failed mass/charge balance: ' + r.id)
    return plans


def repair_preprocessed_pair(model, check_model):
    """Correct working chemistry; synchronize coefficients, not annotations.

The check model can retain input formula metadata. Require its coefficients
to match before mirroring, so later restoration cannot reinstate old electrons.
"""
    plans = [p for r in model.reactions if (p := proposal(r))]
    for p in plans:
        rid = p['reaction']
        if rid not in check_model.reactions:
            raise ValueError('Missing reference reaction for electron correction: ' + rid)
        original = {m.id: v for m, v in check_model.reactions.get_by_id(rid).metabolites.items()}
        if original != p['before']:
            raise ValueError('Reference stoichiometry conflict for ' + rid)
    repair_model(model)
    for p in plans:
        r = check_model.reactions.get_by_id(p['reaction'])
        r.add_metabolites({m: p['after'][m.id] for m in r.metabolites}, combine=False)
    return plans


def validate_final_templates(model, info, control):
    """Detect reintroduction; never silently repair the final exported model."""
    from mqc.control.postgap_control import PostGapValidationError
    remaining = [p for r in model.reactions if (p := proposal(r))]
    reverted = []
    for p in info.get('electron_stoichiometry_preprocessing', []):
        if p['reaction'] not in model.reactions:
            continue
        r = model.reactions.get_by_id(p['reaction'])
        actual = {m.id: v for m, v in r.metabolites.items()}
        if set(actual) != set(p['after']) or any(
            not math.isclose(actual.get(k, float('nan')), c, rel_tol=1e-9, abs_tol=1e-9)
            for k, c in p['after'].items()):
            reverted.append(p['reaction'])
    report = dict(status='failed' if remaining or reverted else 'passed',
                  remaining_legacy_templates=remaining, changed_preprocessing_templates=reverted,
                  scope='exact_supported_single_electron_templates_not_all_redox_chemistry')
    info['electron_stoichiometry_validation'] = control['electron_stoichiometry_validation'] = report
    if report['status'] == 'failed':
        control['final_model_status'] = 'electron_stoichiometry_validation_failed'
        raise PostGapValidationError('Uncorrected or reverted electron template before export')
